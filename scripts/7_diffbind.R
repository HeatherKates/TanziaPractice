library(DiffBind)
library(edgeR)
library(BiocParallel)

# Set up parallelization to use 4 cores
bp_param <- MulticoreParam(workers = 6)

samples <- data.frame(
  SampleID = c("MSH2KO-1", "MSH2KO-2", "MSH2KO-3", "MSH2R4-1", "MSH2R4-2", "MSH2R4-3"),
  Tissue = c(rep("KO", 3), rep("MSH2", 3)),  # You can change 'Tissue' to represent conditions or leave this column out.
  Factor = c(rep("MSH2KO", 3), rep("MSH2R4", 3)),
  Condition = c(rep("MSH2KO", 3), rep("MSH2", 3)),  # The experimental conditions
  Treatment = rep("None", 6),  # Since there's no treatment in this setup
  Replicate = c(1, 2, 3, 1, 2, 3),  # Indicating replicates
  bamReads = c("../6_bigwig/MSH2KO-1_sort_n.bam", "../6_bigwig/MSH2KO-2_sort_n.bam", "../6_bigwig/MSH2KO-3_sort_n.bam",
               "../6_bigwig/MSH2R4-1_sort_n.bam", "../6_bigwig/MSH2R4-2_sort_n.bam", "../6_bigwig/MSH2R4-3_sort_n.bam"),
  Peaks = c("../5_genrich/MSH2KO-1.narrowPeak", "../5_genrich/MSH2KO-2.narrowPeak", "../5_genrich/MSH2KO-3.narrowPeak",
            "../5_genrich/MSH2R4-1.narrowPeak", "../5_genrich/MSH2R4-2.narrowPeak", "../5_genrich/MSH2R4-3.narrowPeak"),
  PeakCaller = rep("narrow", 6),  # Since you're using Genrich, 'narrow' for narrowPeak format
  stringsAsFactors = FALSE
)


# Initialize DiffBind object
dbaObj <- NULL
dbaObj <- dba(sampleSheet=samples)
# Count reads in peaks

# Count reads in consensus peaks using 6 cores
dbaObj <- dba.count(dbaObj, summits=250, bParallel = TRUE,minOverlap = 2)

# Define contrasts (KO vs HLA)
dbaObj <- dba.contrast(dbaObj, categories = DBA_CONDITION)

# Perform differential binding analysis
dbaObj <- dba.analyze(dbaObj)

# Generate and inspect the report of differentially bound regions
diffPeaks <- dba.report(dbaObj)
head(diffPeaks)
write.csv(diffPeaks, file = "../7_diffbind/differential_binding_results.csv")
diffPeaksIn <- read.csv("../7_diffbind/differential_binding_results.csv")


# Create a BED file by selecting the appropriate columns (assuming columns are named 'chr', 'start', 'end')
bed <- data.frame(
  chr = diffPeaksIn$seqnames,       # Replace with the actual column name for chromosome
  start = diffPeaksIn$start,   # Replace with the actual column name for start position
  end = diffPeaksIn$end,       # Replace with the actual column name for end position
  name = diffPeaksIn$X,   # Optionally include peak ID or use `.` if no name is required
  score = diffPeaksIn$Fold   # You can use log2 fold-change or a custom score column
)

# Write the BED file
write.table(bed, file = "../7_diffbind/differential_binding_results.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
