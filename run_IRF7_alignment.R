## Human-bat IRF7 global protein alignment


### Step 1: Enter the file locations

input_fasta <- "/Users/jschoggins/5-Reference/SchogginsLab/ExperimentsData/Bat ISGs/bat IRF7/orthologs/Human vs Palecto align for paper/IRF7_human_bat_input.fasta"
#input_fasta <- "/Users/yourname/path/to/IRF7_human_bat_input.fasta"

output_folder <- "/Users/jschoggins/5-Reference/SchogginsLab/ExperimentsData/Bat ISGs/bat IRF7/orthologs/Human vs Palecto align for paper/output"
#output_folder <- "/Users/yourname/path/to/alignment_results"

### Confirm that the input FASTA exists

if (!file.exists(input_fasta)) {
  stop("Input FASTA file not found. Check the input_fasta path.")
}

### Create the output folder if it does not already exist

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

print(paste("Input FASTA:", input_fasta))
print(paste("Output folder:", output_folder))


### Step 2: Install the required packages if they are missing

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

### Install missing CRAN packages

cran_packages <- c("ggplot2", "svglite")

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

### Install missing Bioconductor packages

bioconductor_packages <- c("Biostrings", "pwalign", "ggmsa")

missing_bioconductor <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioconductor) > 0) {
  BiocManager::install(
    missing_bioconductor,
    ask = FALSE,
    update = FALSE
  )
}

### Confirm that every required package can be loaded

required_packages <- c(cran_packages, bioconductor_packages)

still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(still_missing) > 0) {
  stop(
    "The following packages are still unavailable: ",
    paste(still_missing, collapse = ", ")
  )
}

print("Required packages are available.")


### Step 3: Read the two protein sequences

print("Step 3: Reading the protein sequences.")

sequences <- Biostrings::readAAStringSet(input_fasta)

### Confirm that the FASTA contains exactly two named sequences

if (length(sequences) != 2) {
  stop("The FASTA file must contain exactly two protein sequences.")
}

if (any(is.na(names(sequences))) || any(names(sequences) == "")) {
  stop("Both FASTA sequences must have names.")
}

if (anyDuplicated(names(sequences))) {
  stop("The two FASTA sequence names must be different.")
}

### Confirm that the input contains unaligned sequences

if (any(grepl("-", as.character(sequences), fixed = TRUE))) {
  stop(
    "The input FASTA contains '-' gap characters. ",
    "Use the original unaligned protein sequences."
  )
}

print(paste("Sequence 1:", names(sequences)[1]))
print(paste("Sequence 1 length:", length(sequences[[1]])))
print(paste("Sequence 2:", names(sequences)[2]))
print(paste("Sequence 2 length:", length(sequences[[2]])))


### Step 4: Perform the global protein alignment

print("Step 4: Performing the global protein alignment.")

### Load the BLOSUM62 amino-acid substitution matrix

utils::data("BLOSUM62", package = "pwalign")

### Perform the alignment
### Method: global pairwise alignment
### Substitution matrix: BLOSUM62
### Gap opening penalty: 12
### Gap extension penalty: 4

alignment <- pwalign::pairwiseAlignment(
  pattern = sequences[[1]],
  subject = sequences[[2]],
  type = "global",
  substitutionMatrix = BLOSUM62,
  gapOpening = 12,
  gapExtension = 4
)

### Extract the aligned amino-acid sequences

aligned_first <- as.character(
  pwalign::alignedPattern(alignment)
)

aligned_second <- as.character(
  pwalign::alignedSubject(alignment)
)

### Confirm that both outputs contain the same number of alignment columns

if (nchar(aligned_first) != nchar(aligned_second)) {
  stop("The aligned outputs contain different numbers of alignment columns.")
}

### Remove alignment gaps and confirm that the original sequences are intact

restored_first <- gsub("-", "", aligned_first, fixed = TRUE)
restored_second <- gsub("-", "", aligned_second, fixed = TRUE)

if (!identical(restored_first, as.character(sequences[[1]]))) {
  stop("The first aligned sequence failed validation.")
}

if (!identical(restored_second, as.character(sequences[[2]]))) {
  stop("The second aligned sequence failed validation.")
}

print(paste("Number of alignment columns:", nchar(aligned_first)))
print("Alignment completed and validated successfully.")


### Step 5: Save the aligned FASTA file

print("Step 5: Saving the aligned FASTA file.")

aligned_fasta_file <- file.path(
  output_folder,
  "IRF7_human_bat_alignment.fasta"
)

### Wrap one sequence at 80 characters per FASTA line

wrap_sequence <- function(sequence, line_width = 80) {
  starts <- seq.int(1, nchar(sequence), by = line_width)

  vapply(
    starts,
    function(start) {
      substr(
        sequence,
        start,
        min(nchar(sequence), start + line_width - 1)
      )
    },
    character(1)
  )
}

### Assemble and write the aligned FASTA

aligned_fasta_lines <- c(
  paste0(">", names(sequences)[1]),
  wrap_sequence(aligned_first),
  paste0(">", names(sequences)[2]),
  wrap_sequence(aligned_second)
)

writeLines(aligned_fasta_lines, aligned_fasta_file)

if (!file.exists(aligned_fasta_file)) {
  stop("The aligned FASTA file was not written successfully.")
}

print(paste("Aligned FASTA saved:", aligned_fasta_file))


### Step 6: Create a standard protein-alignment figure

print("Step 6: Creating the protein-alignment figure with ggmsa.")

### Display amino acids according to side-chain chemistry using
### the Chemistry_AA color scheme and wrap the alignment into
### blocks of 60 columns

alignment_plot <- ggmsa::ggmsa(
  msa = aligned_fasta_file,
  font = "DroidSansMono",
  color = "Chemistry_AA",
  char_width = 0.7,
  by_conservation = FALSE,
  seq_name = TRUE,
  show.legend = TRUE
) +
  ggmsa::facet_msa(field = 60) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      title = "Amino-acid chemistry",
      nrow = 1,
      byrow = TRUE
    )
  ) +
  ggplot2::labs(
    title = "Human and bat IRF7 protein sequence alignment",
    subtitle = "Global alignment | BLOSUM62 | gap opening 12 | gap extension 4"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 6.5),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7)
  )

### Compatibility fix for ggmsa with modern versions of ggplot2:
### restore the fill mapping required to display the color legend


alignment_plot$layers[[1]]$mapping$fill <-
  ggplot2::aes(fill = color)$fill

print("Protein-alignment figure created.")

### Step 7: Save the figure as SVG, PDF, and PNG files

print("Step 7: Saving the alignment figures.")

for (figure_format in c("svg", "pdf", "png")) {
  figure_file <- file.path(
    output_folder,
    paste0("IRF7_human_bat_alignment_ggmsa.", figure_format)
  )

  ggplot2::ggsave(
    filename = figure_file,
    plot = alignment_plot,
    width = 7,
    height = 8,
    units = "in",
    dpi = 400,
    bg = "white"
  )

  if (!file.exists(figure_file)) {
    stop("The figure was not written successfully: ", figure_file)
  }

  print(paste("Figure saved:", figure_file))
}


### Step 8: Save information about the R session

session_info_file <- file.path(
  output_folder,
  "sessionInfo.txt"
)

utils::capture.output(
  utils::sessionInfo(),
  file = session_info_file
)

if (!file.exists(session_info_file)) {
  stop("The session information file was not written successfully.")
}

print(paste("Session information saved:", session_info_file))


### Analysis complete

print("SUCCESS: The alignment analysis is complete.")
print(paste("All results were saved in:", output_folder))
