# FASTA Sequence Analysis Pipeline

Automated Bash pipeline for fetching NCBI FASTA gene sequences, calculating sequence statistics, and generating reverse complement sequences for multiple gene accession IDs.

## Features
- Automatically downloads FASTA sequences from NCBI using Accession IDs.
- Calculates GC/AT percentage, N50 value, and sequence length.
- Computes window-based GC content statistics.
- Generates reverse complement sequences.

## Usage
Run the script with any NCBI gene accession ID (e.g.,'NM_000059', `NM_0000545`, `NM_000786`, `NM_000799`):
> cat << 'EOF' > README.md
# FASTA Sequence Analysis Pipeline

Automated Bash pipeline for fetching NCBI FASTA gene sequences, calculating sequence statistics, and generating reverse complement sequences for multiple gene accession IDs.

## Features
- Automatically downloads FASTA sequences from NCBI using Accession IDs.
- Calculates GC/AT percentage, N50 value, and sequence length.
- Computes window-based GC content statistics.
- Generates reverse complement sequences.

## Usage
Run the script with any NCBI gene accession ID (e.g.,'NM_000059', `NM_0000545`, `NM_000786`, `NM_000799`):
