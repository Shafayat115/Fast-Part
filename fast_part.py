import os
import random
import time
import argparse
from Bio import SeqIO
from collections import defaultdict

# Parsing command-line arguments
def get_args():
    parser = argparse.ArgumentParser(description="A tool for efficient partitioning and clustering of sequences.")
    parser.add_argument("--fasta_file", required=True, help="Input FASTA file containing sequences")
    parser.add_argument("--output_dir", required=True, help="Directory to store output files")
    parser.add_argument("--method", required=True, choices=["mmseq", "cdhit"], help="Clustering method: 'mmseq' or 'cdhit'")
    parser.add_argument("--identity_threshold", type=float, default=0.8, help="Identity threshold for clustering")
    parser.add_argument("--cdhit_cores", type=int, default=128, help="Number of cores to use for CD-HIT (if used)")
    parser.add_argument("--train_ratio", type=float, default=0.8, help="Train set ratio for splitting clusters (5% less is used)")
    parser.add_argument("--query_cover", type=int, default=60, help="Query coverage for DIAMOND alignment")
    return parser.parse_args()

# Determine word length based on identity threshold for CD-HIT
def get_word_length(identity_threshold):
    if 0.4 <= identity_threshold < 0.5:
        return 2
    elif 0.5 <= identity_threshold < 0.6:
        return 3
    elif 0.6 <= identity_threshold < 0.7:
        return 4
    else:
        return 5

# Step 1: Division - Parsing and Separating the Fasta File by Labels
def parse_and_separate_fasta(fasta_file):
    sequences_by_label = defaultdict(list)
    for record in SeqIO.parse(fasta_file, "fasta"):
        label = record.description.split('|')[1]
        sequences_by_label[label].append(record)
    return sequences_by_label

def write_sequences_by_label(sequences_by_label, output_dir):
    label_dirs = []
    for label, sequences in list(sequences_by_label.items()):
        sanitized_label = label.replace("/", "_").replace("\\", "_")
        label_dir = os.path.join(output_dir, sanitized_label)
        os.makedirs(label_dir, exist_ok=True)
        output_file = os.path.join(label_dir, f"{sanitized_label}.fasta")
        SeqIO.write(sequences, output_file, "fasta")
        label_dirs.append(label_dir)
        yield sanitized_label, output_file
    return label_dirs

# Step 2: Clustering with CD-HIT or MMseqs2
def run_cdhit(input_file, output_file, identity_threshold, cores):
    word_length = get_word_length(identity_threshold)
    cmd = f"cd-hit -i {input_file} -o {output_file} -c {identity_threshold} -n {word_length} -T {cores}"
    os.system(cmd)

def run_mmseqs2_clustering(input_file, output_file, tmp_dir, identity_threshold):
    os.makedirs(tmp_dir, exist_ok=True)
    os.system(f"mmseqs createdb {input_file} {output_file}_db")
    os.system(f"mmseqs cluster {output_file}_db {output_file}_cluster {tmp_dir} --min-seq-id {identity_threshold} -c 0.5 --cov-mode 1")
    os.system(f"mmseqs createtsv {output_file}_db {output_file}_db {output_file}_cluster {output_file}.tsv")
    os.system(f"rm -r {tmp_dir}")

def parse_clusters(method, cluster_file):
    clusters = defaultdict(list)
    if method == "cdhit":
        with open(cluster_file, 'r') as file:
            current_cluster = None
            for line in file:
                if line.startswith('>Cluster'):
                    current_cluster = line.strip().split()[-1]
                elif current_cluster:
                    sequence_id = line.split('>')[1].split('|')[0]
                    clusters[current_cluster].append(sequence_id)
    elif method == "mmseq":
        with open(cluster_file, 'r') as file:
            cluster_ids = {}
            current_cluster_id = 1
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                cluster_key = parts[0]
                if cluster_key not in cluster_ids:
                    cluster_ids[cluster_key] = current_cluster_id
                    current_cluster_id += 1
                sequence_id = parts[1].split('|')[0]
                clusters[cluster_ids[cluster_key]].append(sequence_id)
    return clusters

# Step 3: Partitioning Clusters into Train and Test Sets
def split_train_test_by_size(cluster_sizes, train_ratio):
    adjusted_train_ratio = train_ratio - 0.05  # Adjust train ratio to 5% less
    total_sequences = sum(cluster_sizes.values())
    if len(cluster_sizes) == 1:
        return list(cluster_sizes.keys()), []

    sorted_clusters = sorted(cluster_sizes, key=cluster_sizes.get)
    train_clusters, test_clusters = [], []
    current_train_sequences = 0

    if sorted_clusters:
        test_clusters.append(sorted_clusters.pop(0))
        if sorted_clusters:
            train_clusters.append(sorted_clusters.pop(-1))
            current_train_sequences += cluster_sizes[train_clusters[-1]]

    random.shuffle(sorted_clusters)
    for cluster in sorted_clusters:
        if (current_train_sequences / total_sequences) < adjusted_train_ratio:
            train_clusters.append(cluster)
            current_train_sequences += cluster_sizes[cluster]
        else:
            test_clusters.append(cluster)

    return train_clusters, test_clusters

# Writing Sequences to Fasta File Based on Cluster Assignments
def write_sequences_to_fasta(fasta_file, clusters, output_file):
    sequence_ids = {seq_id for cluster_id in clusters for seq_id in clusters[cluster_id]}
    with open(output_file, 'w') as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id.split('|')[0] in sequence_ids:
                SeqIO.write(record, out_fasta, "fasta")

# Step 4: Running DIAMOND Alignment and Iterative Reassignment
def create_diamond_db(train_fasta, db_name):
    os.system(f"diamond makedb --in {train_fasta} -d {db_name}")

def run_diamond_alignment(test_fasta, db_name, output_file, identity_threshold, query_cover):
    diamond_identity = identity_threshold * 100  # Convert identity threshold to percentage
    os.system(f"diamond blastp -d {db_name} -q {test_fasta} --id {diamond_identity} --query-cover {query_cover} -o {output_file} --quiet")
    return os.path.getsize(output_file) != 0

def filter_diamond_output(diamond_output):
    qualifying_sequences, reference_sequences = set(), set()
    with open(diamond_output, 'r') as file:
        for line in file:
            parts = line.strip().split()
            qualifying_sequences.add(parts[0].split('|')[0])
            reference_sequences.add(parts[1].split('|')[0])
    return qualifying_sequences, reference_sequences

def iterative_diamond_alignment_and_reassignment(train_clusters, test_clusters, clusters, fasta_file, label_dir, train_ratio, identity_threshold, query_cover):
    train_output_file = os.path.join(label_dir, "train.fasta")
    test_output_file = os.path.join(label_dir, "test.fasta")
    write_sequences_to_fasta(fasta_file, {cl: clusters[cl] for cl in train_clusters}, train_output_file)
    write_sequences_to_fasta(fasta_file, {cl: clusters[cl] for cl in test_clusters}, test_output_file)

    while True:
        if os.path.getsize(train_output_file) == 0 or os.path.getsize(test_output_file) == 0:
            break

        db_name = os.path.join(label_dir, "train_db")
        diamond_output = os.path.join(label_dir, "diamond_matches.m8")
        create_diamond_db(train_output_file, db_name)

        success = run_diamond_alignment(test_output_file, db_name, diamond_output, identity_threshold, query_cover)
        if not success:
            break

        qualifying_sequences, reference_sequences = filter_diamond_output(diamond_output)
        if not qualifying_sequences:
            break

        clusters_to_reassign = {cluster_id for cluster_id in test_clusters if set(clusters[cluster_id]).intersection(qualifying_sequences)}
        new_test_clusters = set(test_clusters) - clusters_to_reassign
        new_train_clusters = set(train_clusters).union(clusters_to_reassign)

        total_sequences = sum(len(clusters[cl]) for cl in clusters)
        new_train_set_size = sum(len(clusters[cl]) for cl in new_train_clusters)

        if new_train_set_size / total_sequences > train_ratio or len(new_test_clusters) < 2:
            break

        train_clusters, test_clusters = list(new_train_clusters), list(new_test_clusters)
        write_sequences_to_fasta(fasta_file, {cl: clusters[cl] for cl in train_clusters}, train_output_file)
        write_sequences_to_fasta(fasta_file, {cl: clusters[cl] for cl in test_clusters}, test_output_file)

    return train_clusters, test_clusters

# Function to remove sequences from a FASTA file
def remove_sequences_from_fasta(fasta_file, sequences_to_remove, output_file):
    sequences = (record for record in SeqIO.parse(fasta_file, "fasta") if record.id.split('|')[0] not in sequences_to_remove)
    SeqIO.write(sequences, output_file, "fasta")

def main():
    args = get_args()

    # Automatically switch to mmseq mode if identity threshold is below 0.4
    if args.identity_threshold < 0.4:
        args.method = "mmseq"
    
    start_time = time.time()
    os.makedirs(args.output_dir, exist_ok=True)
    temp_dir = os.path.join(args.output_dir, "Temp")
    os.makedirs(temp_dir, exist_ok=True)

    # Define summary_file early so it can be used in the cleanup step later
    summary_file = os.path.join(args.output_dir, "summary.txt")
    
    # Write configuration settings to summary file
    with open(summary_file, 'w') as summary:
        summary.write("Configuration Summary:\n")
        summary.write(f"Clustering Method: {args.method}\n")
        summary.write(f"Identity Threshold: {args.identity_threshold}\n")
        summary.write(f"Cores: {args.cdhit_cores}\n")
        summary.write(f"Train Ratio (adjusted): {args.train_ratio - 0.05}\n")
        summary.write(f"Query Cover: {args.query_cover}\n\n")

    sequences_by_label = parse_and_separate_fasta(args.fasta_file)
    initial_sequence_count = len(list(SeqIO.parse(args.fasta_file, "fasta")))

    combined_train_sequences = []
    combined_test_sequences = []
    label_sequence_counts = {}
    train_label_counts = defaultdict(int)
    test_label_counts = defaultdict(int)

    for label, fasta_file in write_sequences_by_label(sequences_by_label, args.output_dir):
        label_dir = os.path.join(args.output_dir, label)
        output_file = os.path.join(temp_dir, f"{label}_{args.method}")
        
        # Select clustering method
        if args.method == "cdhit":
            run_cdhit(fasta_file, output_file, args.identity_threshold, args.cdhit_cores)
            cluster_file = f"{output_file}.clstr"
        elif args.method == "mmseq":
            tmp_dir = os.path.join(label_dir, "tmp")
            run_mmseqs2_clustering(fasta_file, output_file, tmp_dir, args.identity_threshold)
            cluster_file = f"{output_file}.tsv"
        
        clusters = parse_clusters(args.method, cluster_file)
        cluster_sizes = {cluster: len(ids) for cluster, ids in clusters.items()}

        label_sequence_counts[label] = len(sequences_by_label[label])

        train_clusters, test_clusters = split_train_test_by_size(cluster_sizes, args.train_ratio)
        final_train_clusters, final_test_clusters = iterative_diamond_alignment_and_reassignment(
            train_clusters, test_clusters, clusters, fasta_file, temp_dir, args.train_ratio, args.identity_threshold, args.query_cover)

        train_output_file = os.path.join(temp_dir, f"{label}_train.fasta")
        test_output_file = os.path.join(temp_dir, f"{label}_test.fasta")
        write_sequences_to_fasta(fasta_file, {cl: clusters[cl] for cl in final_train_clusters}, train_output_file)
        write_sequences_to_fasta(fasta_file, {cl: clusters[cl] for cl in final_test_clusters}, test_output_file)

        train_label_counts[label] += len(list(SeqIO.parse(train_output_file, "fasta")))
        test_label_counts[label] += len(list(SeqIO.parse(test_output_file, "fasta")))

        combined_train_sequences.extend(SeqIO.parse(train_output_file, "fasta"))
        combined_test_sequences.extend(SeqIO.parse(test_output_file, "fasta"))

    # Write combined train and test sets
    train_file = os.path.join(args.output_dir, "Train.fasta")
    test_file = os.path.join(args.output_dir, "Test.fasta")
    SeqIO.write(combined_train_sequences, train_file, "fasta")
    SeqIO.write(combined_test_sequences, test_file, "fasta")

    # Calculate removed sequences
    result_file = os.path.join(temp_dir, "result_ARG.txt")
    create_diamond_db(train_file, os.path.join(temp_dir, "train_db"))
    run_diamond_alignment(test_file, os.path.join(temp_dir, "train_db"), result_file, args.identity_threshold, args.query_cover)
    removed_sequences = list(filter_diamond_output(result_file)[1])

    updated_train_file = os.path.join(args.output_dir, "Train.fasta")
    remove_sequences_from_fasta(train_file, removed_sequences, updated_train_file)

    # Find missing labels in train and test sets
    all_labels = set(label_sequence_counts.keys())
    train_labels = set(train_label_counts.keys())
    test_labels = set(test_label_counts.keys())
    missing_train_labels = all_labels - train_labels
    missing_test_labels = all_labels - test_labels

    # Remove temporary files and directories
    for root, dirs, files in os.walk(temp_dir, topdown=False):
        for file in files:
            os.remove(os.path.join(root, file))
        for dir in dirs:
            os.rmdir(os.path.join(root, dir))
    os.rmdir(temp_dir)

    # Delete all remaining files in output directory except Train, Test, and summary.txt
    for item in os.listdir(args.output_dir):
        item_path = os.path.join(args.output_dir, item)
        if item_path not in {train_file, test_file, summary_file}:
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                for root, dirs, files in os.walk(item_path, topdown=False):
                    for file in files:
                        os.remove(os.path.join(root, file))
                    for dir in dirs:
                        os.rmdir(os.path.join(root, dir))
                os.rmdir(item_path)

    end_time = time.time()
    elapsed_time = end_time - start_time

    # Write summary file with added information
    with open(summary_file, 'a') as summary:
        summary.write(f"Initial total sequence count: {initial_sequence_count}\n")
        summary.write("Sequence counts per label:\n")
        for label, count in label_sequence_counts.items():
            summary.write(f"  {label}: {count}\n")
        summary.write(f"\nFinal sequence counts per label in Train set:\n")
        for label, count in train_label_counts.items():
            percentage = (count / label_sequence_counts[label]) * 100 if label_sequence_counts[label] > 0 else 0
            summary.write(f"  {label}: {count} ({percentage:.2f}%)\n")
        summary.write(f"\nFinal sequence counts per label in Test set:\n")
        for label, count in test_label_counts.items():
            percentage = (count / label_sequence_counts[label]) * 100 if label_sequence_counts[label] > 0 else 0
            summary.write(f"  {label}: {count} ({percentage:.2f}%)\n")
        summary.write(f"\nTotal sequences removed: {len(removed_sequences)}\n")
        summary.write("Removed sequence IDs:\n")
        for seq_id in removed_sequences:
            summary.write(f"  {seq_id}\n")
        summary.write(f"\nNumber of labels in Train set: {len(train_labels)}\n")
        summary.write(f"Number of labels in Test set: {len(test_labels)}\n")
        summary.write("Missing labels in Train set:\n")
        for label in missing_train_labels:
            summary.write(f"  {label}\n")
        summary.write("Missing labels in Test set:\n")
        for label in missing_test_labels:
            summary.write(f"  {label}\n")
        summary.write(f"\nExecution time: {elapsed_time:.2f} seconds\n")

    print(f"Updated combined train file is saved to {updated_train_file}")
    print(f"Updated combined test file is saved to {test_file}")
    print(f"Summary file is saved to {summary_file}")
    print(f"Execution time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
