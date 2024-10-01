library(DiffBind)
library(edgeR)
library(BiocParallel)

# Set up parallelization to use 4 cores
bp_param <- MulticoreParam(workers = 6)

samples <- data.frame(
  SampleID = c("KO1", "KO2", "KO3", "MSH2R4_1", "MSH2R4_2", "MSH2R4_3"),
  Tissue = c(rep("KO", 3), rep("MSH2", 3)),  # You can change 'Tissue' to represent conditions or leave this column out.
  Factor = c(rep("MSH2KO", 3), rep("MSH2R4", 3)),
  Condition = c(rep("KO", 3), rep("MSH2", 3)),  # The experimental conditions
  Treatment = rep("None", 6),  # Since there's no treatment in this setup
  Replicate = c(1, 2, 3, 1, 2, 3),  # Indicating replicates
  bamReads = c("../4a_samtools/MSH2KO-1.sorted.bam", "../4a_samtools/MSH2KO-2.sorted.bam", "../4a_samtools/MSH2KO-3.sorted.bam",
               "../4a_samtools/MSH2R4-1.sorted.bam", "../4a_samtools/MSH2R4-2.sorted.bam", "../4a_samtools/MSH2R4-3.sorted.bam"),
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
dbaObj <- dba.count(dbaObj, summits=250, bParallel = TRUE)

# Define contrasts (KO vs HLA)
dbaObj <- dba.contrast(dbaObj, categories = DBA_CONDITION)

# Perform differential binding analysis
dbaObj <- dba.analyze(dbaObj)

# Generate and inspect the report of differentially bound regions
diffPeaks <- dba.report(dbaObj)
head(diffPeaks)
write.csv(diffPeaks, file = "differential_binding_results.csv")



