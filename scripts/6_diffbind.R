library(DiffBind)
library(edgeR)
library(BiocParallel)

# Set up parallelization to use 4 cores
bp_param <- MulticoreParam(workers = 6)

# Create metadata df with BAM file paths added
samples <- data.frame(
  SampleID = c("KO1", "KO2", "KO3", "HLA1", "HLA2", "HLA3"),
  Condition = c(rep("KO", 3), rep("HLA", 3)),
  PeakFile = c("../5_genrich/MSH2KO-1.narrowPeak", "../5_genrich/MSH2KO-2.narrowPeak", "../5_genrich/MSH2KO-3.narrowPeak",
               "../5_genrich/MSH2R4-1.narrowPeak", "../5_genrich/MSH2R4-2.narrowPeak", "../5_genrich/MSH2R4-3.narrowPeak"),
  BamFile = c("../4a_samtools/MSH2KO-1.sorted.bam", "../4a_samtools/MSH2KO-2.sorted.bam", "../4a_samtools/MSH2KO-3.sorted.bam",
              "../4a_samtools/MSH2R4-1.sorted.bam", "../4a_samtools/MSH2R4-2.sorted.bam", "../4a_samtools/MSH2R4-3.sorted.bam"),
  stringsAsFactors = FALSE
)

# Initialize DiffBind object
dbaObj <- NULL

# Load peaks and BAM files for each sample
for (i in 1:nrow(samples)) {
  dbaObj <- dba.peakset(dbaObj, peaks = samples$PeakFile[i], bamReads = samples$BamFile[i], 
                        sampID = samples$SampleID[i], tissue = samples$Condition[i], peak.format = "narrow")
}

# Count reads in consensus peaks using 6 cores
dbaObj <- dba.count(dbaObj, bParallel = TRUE)

# Define contrasts (KO vs HLA)
dbaObj <- dba.contrast(dbaObj, categories = DBA_CONDITION)

# Perform differential binding analysis
dbaObj <- dba.analyze(dbaObj)

# Generate and inspect the report of differentially bound regions
diffPeaks <- dba.report(dbaObj)
head(diffPeaks)
write.csv(diffPeaks, file = "differential_binding_results.csv")



