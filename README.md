# Fast_Part

**Fast_Part** is an efficient tool for partitioning and clustering sequences. It uses CD-HIT or MMseqs2 for clustering, along with DIAMOND for alignments, to enable streamlined handling of biological sequences, automated partitioning into training and test sets, and iterative reassignment based on DIAMOND alignments. 

## Features

- **Clustering Options**: Supports clustering with **CD-HIT** or **MMseqs2**.
- **Automatic Word Length Adjustment**: For CD-HIT, automatically adjusts word length based on the identity threshold.
- **Default to MMseqs2**: Automatically defaults to MMseqs2 if the identity threshold is set below 0.4.
- **Train-Test Partitioning**: Automatically applies a 5% adjustment to the specified train ratio.
- **Detailed Summary**: Outputs a summary file with configuration details, partition sizes, removed sequences, and execution time.

## Requirements

- **Python 3.6+**
- **CD-HIT**
- **MMseqs2**
- **DIAMOND**
- **Biopython** library (`pip install biopython`)

Ensure that **CD-HIT**, **MMseqs2**, and **DIAMOND** are installed and accessible in your system’s PATH.

## Installation Instructions

Before using **Fast_Part**, ensure that **CD-HIT**, **MMseqs2**, and **DIAMOND** are installed on your system. Below are installation instructions for each tool.

### CD-HIT

CD-HIT is a clustering tool for protein and nucleotide sequences.

1. **Download CD-HIT**:
   - Go to the [CD-HIT website](https://github.com/weizhongli/cdhit) and download the latest release.
   
2. **Install CD-HIT**:
   - Extract the downloaded file and navigate to the extracted folder.
   - Run the following commands:
     ```bash
     make
     sudo make install
     ```
   - Ensure the `cd-hit` binary is accessible in your PATH. You can test it with:
     ```bash
     cd-hit --help
     ```

### MMseqs2

MMseqs2 (Many-against-Many sequence searching) is used for fast and sensitive sequence similarity searches and clustering.

1. **Download MMseqs2**:
   - Visit the [MMseqs2 GitHub repository](https://github.com/soedinglab/MMseqs2) and download the appropriate release for your operating system.
   
2. **Install MMseqs2**:
   - Extract the downloaded file.
   - Add MMseqs2 to your PATH or move it to a directory already in your PATH.
   - You can verify the installation with:
     ```bash
     mmseqs --help
     ```

### DIAMOND

DIAMOND is a fast sequence aligner for protein and translated DNA searches.

1. **Download DIAMOND**:
   - Go to the [DIAMOND GitHub page](https://github.com/bbuchfink/diamond) and download the latest release for your platform.
   
2. **Install DIAMOND**:
   - Extract the downloaded file.
   - Move the `diamond` executable to a directory in your PATH.
   - Test the installation with:
     ```bash
     diamond --help
     ```


## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/Shafayat115/Fast-Part.git
   cd Fast_Part
2. Install dependencies:
   ```bash
   pip install biopython
## Usage
Run Fast_Part with a specified fasta_file and output_dir, along with other optional parameters. By default, the tool uses cdhit with an identity threshold of 0.8, 128 CD-HIT cores, a train ratio of 0.8 (adjusted by 5% internally), and a query coverage of 60 for DIAMOND alignment.

## Command-Line Arguments

| Argument             | Description                                                 | Default    |
|----------------------|-------------------------------------------------------------|------------|
| `--fasta_file`       | Input FASTA file containing sequences                       | *Required* |
| `--output_dir`       | Directory to store output files                             | *Required* |
| `--method`           | Clustering method (`cdhit` or `mmseq`)                      | *Required* |
| `--identity_threshold` | Identity threshold for clustering                          | `0.8`      |
| `--cdhit_cores`      | Number of cores for CD-HIT (only if `cdhit` method is used) | `128`      |
| `--train_ratio`      | Train set ratio for splitting clusters (adjusted by 5%)     | `0.8`      |
| `--query_cover`      | Query coverage for DIAMOND alignment                        | `60`       |



## Example 1: Using CD-HIT with default settings
```bash
python fast_part.py --fasta_file path/to/input.fasta --output_dir path/to/output --method cdhit
```
## Example 2: Using MMseqs2 with custom identity threshold

```bash
python fast_part.py --fasta_file path/to/input.fasta --output_dir path/to/output --method mmseq --identity_threshold 0.7
```
## Example 3: Specifying custom cores for CD-HIT
```bash
python fast_part.py --fasta_file path/to/input.fasta --output_dir path/to/output --method cdhit --cdhit_cores 64
```
## Example 4: Running with a custom train ratio and query cover
```bash
python fast_part.py --fasta_file path/to/input.fasta --output_dir path/to/output --method cdhit --train_ratio 0.85 --query_cover 70
```
## Output Files

**Train.fasta:** Contains sequences assigned to the training set.

**Test.fasta:** Contains sequences assigned to the test set.

**summary.txt:** A summary file with configuration details, initial and final sequence counts per label, removed sequences, missing labels, and execution time.
