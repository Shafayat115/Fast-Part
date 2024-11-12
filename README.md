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

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/Fast_Part.git
   cd Fast_Part
