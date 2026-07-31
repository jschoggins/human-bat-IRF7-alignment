# Human–bat IRF7 protein sequence alignment

This repository contains the R script and input protein sequences used to
generate a global pairwise alignment of human IRF7 (`NP_001563.2`) and bat
IRF7 (`NP_001307207.1`). The script saves the aligned sequences and renders the
alignment as SVG, PDF, and PNG figures.

## Repository contents

- `run_IRF7_alignment.R`: alignment, validation, and figure-generation script
- `IRF7_human_bat_input.fasta`: two unaligned IRF7 protein sequences

## Requirements

- A recent version of R
- RStudio (recommended for novice users)
- An internet connection the first time the required packages are installed

The script checks for the required CRAN and Bioconductor packages and installs
any that are missing.

## Run the analysis

1. Download the R script and FASTA file from this repository.
2. Open `run_IRF7_alignment.R` in RStudio.
3. At the beginning of the script, replace the example `input_fasta` and
   `output_folder` paths with their locations on your computer.
4. Run the script from top to bottom in RStudio.

The input FASTA must contain exactly two uniquely named, unaligned protein
sequences. The script creates the output folder automatically if it does not
already exist.

## Output files

The output folder will contain:

- `IRF7_human_bat_alignment.fasta`: aligned protein sequences
- `IRF7_human_bat_alignment_ggmsa.svg`: editable vector figure
- `IRF7_human_bat_alignment_ggmsa.pdf`: PDF figure
- `IRF7_human_bat_alignment_ggmsa.png`: 400-dpi PNG figure
- `sessionInfo.txt`: R and package versions used for the analysis

## Alignment and visualization

The two proteins are globally aligned with `pwalign::pairwiseAlignment()` using
the BLOSUM62 substitution matrix and affine-gap parameters of 12 for gap
opening and 4 for gap extension. The script validates the alignment by removing
all introduced gaps and confirming that both original sequences are recovered
exactly.

The alignment is visualized with `ggmsa` using its `Chemistry_AA` color scheme,
which groups amino acids according to side-chain chemistry. The figure is
wrapped into blocks of 60 alignment columns.

The script includes a small compatibility correction that restores the
`ggmsa` color legend with modern versions of `ggplot2`.
