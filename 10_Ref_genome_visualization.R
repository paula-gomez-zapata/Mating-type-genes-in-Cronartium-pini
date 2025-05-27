################################################################################
# 09 Reference genome. Visualization.
# Plot Circular representation of gene and repeat distribution across scaffolds that 
# cumulatively contain 85% of the reference genome of C. pini.
# Plot zoom in scaffolds where MAT regions are found.
################################################################################
# Clear global environment
rm(list = ls())

# Load required libraries
library(circlize)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(GenomicFeatures)
library(ggplot2)
library(rtracklayer)

# Set working directory
setwd("C:/Users/pago0001/Documents/Personal Documents/Paula's notes Meetings Project Progress Follow-up/Report 28")

# Load scaffold length data
scaffold_data <- read.table("scaffold_length_nucleargenome.txt", header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(scaffold_data) <- c("Chrom", "stop")

scaffold_data$Chrom <- as.character(scaffold_data$Chrom)
scaffold_data$stop <- as.numeric(scaffold_data$stop)

# Sort scaffolds by size in descending order
scaffold_data <- scaffold_data %>% arrange(desc(stop))

# Calculate cumulative genome coverage
scaffold_data_cum <- scaffold_data %>%
  mutate(CumulativeLength = cumsum(stop),
         TotalGenomeSize = sum(stop),
         CumulativePercentage = (CumulativeLength / TotalGenomeSize) * 100)

# Find the number of scaffolds needed to reach 85% of the genome
num_scaffolds_85 <- scaffold_data_cum %>% filter(CumulativePercentage <= 85) %>% nrow()

# Plot cumulative genome coverage
cumulative_plot <- ggplot(scaffold_data_cum, aes(x = 1:nrow(scaffold_data_cum), y = CumulativePercentage)) +
  geom_line(color = "blue", size = 1.0) +
  geom_hline(yintercept = 85, linetype = "dashed", color = "black") +
  geom_vline(xintercept = num_scaffolds_85, linetype = "dashed", color = "black") +
  labs(x = "Number of Scaffolds",
       y = "Cumulative Length (%)") +
  theme_minimal()

ggsave("cumulative_genome_coverage.svg", plot = cumulative_plot, width = 10, height = 15, dpi = 300)


head(scaffold_data)
scaffold_data$start <- 1
scaffold_data <- scaffold_data %>% dplyr::select(Chrom, start, stop)


# Select only the top 34 scaffolds (change this number as needed)
top_scaffolds <- scaffold_data %>% head(34)

# Load Repeat Annotations
repeats <- read.table('Cpini.filteredRepeats_nucleargenome.sorted.gff', sep = '\t', header = FALSE)
repeats <- repeats[c('V1', 'V4', 'V5')]
colnames(repeats) <- c('Chrom', 'start', 'stop')
repeats <- repeats %>% filter(Chrom %in% top_scaffolds$Chrom)
repeats$Chrom <- as.character(repeats$Chrom)
repeats$start <- as.numeric(repeats$start)
repeats$stop <- as.numeric(repeats$stop)
head(repeats)

# Load Gene Annotations
genes <- read.table('C_flacc.ipr.sorted_nucleargenome_edited.gff3', sep = '\t', header = FALSE)
head(genes)
genes <- genes[c('V1', 'V3', 'V4', 'V5')]
# Extract genes only
genes_only <- genes %>%
  filter(V3 %in% c("gene"))

head(genes_only)
genes_only <- genes_only[c('V1', 'V4', 'V5')]
colnames(genes_only) <- c('Chrom', 'start', 'stop')

genes_only <- genes_only %>% filter(Chrom %in% top_scaffolds$Chrom)
genes_only$Chrom <- as.character(genes_only$Chrom)
genes_only$start <- as.numeric(genes_only$start)
genes_only$stop <- as.numeric(genes_only$stop)
str(genes_only)

# Create a function to calculate percentage coverage of genes/repeats per window
calculate_coverage <- function(features, window_size = 10000) {
  coverage_data <- data.frame()
  
  # Ensure scaffold names are character strings
  scaffold_data$Chrom <- as.character(scaffold_data$Chrom)
  features$Chrom <- as.character(features$Chrom)
  
  # Ensure scaffolds in both datasets match
  valid_scaffolds <- intersect(scaffold_data$Chrom, features$Chrom)
  features <- features %>% filter(Chrom %in% valid_scaffolds)
  scaffold_data <- scaffold_data %>% filter(Chrom %in% valid_scaffolds)
  
  for (chrom in unique(features$Chrom)) {
    chrom_length <- scaffold_data$stop[scaffold_data$Chrom == chrom]
    
    if (length(chrom_length) == 0) next  # Skip if scaffold is not found
    
    # Generate windowed regions for the scaffold
    windows <- data.frame(
      Chrom = chrom,
      start = seq(1, chrom_length, by = window_size),
      stop = pmin(seq(1, chrom_length, by = window_size) + window_size - 1, chrom_length)
    )
    
    # Ensure features exist for this scaffold
    features_filtered <- features %>% filter(Chrom == chrom)
    if (nrow(features_filtered) == 0) next  
    
    # Convert features and windows into GRanges objects
    feature_gr <- GRanges(seqnames = features_filtered$Chrom, 
                          ranges = IRanges(start = features_filtered$start, end = features_filtered$stop))
    
    window_gr <- GRanges(seqnames = windows$Chrom, 
                         ranges = IRanges(start = windows$start, end = windows$stop))
    
    # Ensure seqlevels match exactly
    seqlevels(feature_gr) <- seqlevels(window_gr) <- unique(c(seqlevels(feature_gr), seqlevels(window_gr)))
    
    # Find overlaps between features and windows
    overlaps <- findOverlaps(window_gr, feature_gr)
    
    # Compute base coverage for each window
    coverage_values <- rep(0, length(window_gr))
    
    for (i in seq_along(window_gr)) {
      overlapping_features <- subjectHits(overlaps)[queryHits(overlaps) == i]
      
      if (length(overlapping_features) > 0) {
        # Fix: Process each feature separately using a loop
        covered_bases <- 0
        for (feat in overlapping_features) {
          covered_bases <- covered_bases + sum(width(pintersect(window_gr[i], feature_gr[feat])))
        }
        
        # Fix: Normalize coverage correctly per window
        window_size_actual <- width(window_gr[i])  # Actual window size (handles last window cases)
        coverage_values[i] <- min((covered_bases / window_size_actual) * 100, 100)  # Now capped at 100%
      }
    }
    
    # Store coverage data
    windows$coverage <- coverage_values
    coverage_data <- rbind(coverage_data, windows)
  }
  
  return(coverage_data)
}


# Calculate coverage for genes and repeats
gene_coverage <- calculate_coverage(genes_only, 10000)
head(gene_coverage)
max(gene_coverage$coverage)

repeat_coverage <- calculate_coverage(repeats, 10000)
head(repeat_coverage)
max(repeat_coverage$coverage)

# Normalize coverage values
gene_coverage$coverage <- gene_coverage$coverage / max(gene_coverage$coverage, na.rm = TRUE)
repeat_coverage$coverage <- repeat_coverage$coverage / max(repeat_coverage$coverage, na.rm = TRUE)

# Load telomere repeats
telomeres <- read.table('telomere_repeats_summary.csv', sep = ',', header = TRUE)
colnames(telomeres) <- c("Chrom", "start_repeats", "end_repeats")
telomeres <- telomeres %>% filter(Chrom %in% top_scaffolds$Chrom)
str(telomeres)
telomeres$Chrom <- as.character(telomeres$Chrom)
telomeres$start_repeats <- as.numeric(telomeres$start_repeats)
telomeres$end_repeats <- as.numeric(telomeres$end_repeats)

# Merge telomere data with scaffold_data to get scaffold lengths (stop positions)
telomeres_merge <- telomeres %>%
  left_join(top_scaffolds, by = "Chrom")

# Reshape data to have each chromosome twice (start and end telomere)
telomeres_long <- telomeres_merge %>%
  pivot_longer(cols = c(start_repeats, end_repeats), 
               names_to = "telomere_end", values_to = "value") %>%
  mutate(
    start = if_else(telomere_end == "start_repeats", 1, stop-1),
    stop  = if_else(telomere_end == "start_repeats", 2, stop)
  ) %>%
  dplyr::select(Chrom, start, stop, value)

head(telomeres_long)

# Remove NAs and ensure numeric data
telomeres_long <- na.omit(telomeres_long)
telomeres_long$Chrom <- as.character(telomeres_long$Chrom)
telomeres_long$start <- as.numeric(telomeres_long$start)
telomeres_long$stop <- as.numeric(telomeres_long$stop)
telomeres_long$value <- as.numeric(telomeres_long$value)

# Filter out rows where 'value' is 0
telomere_repeats_filter <- telomeres_long %>%
  filter(value != 0)

# Print the filtered dataframe
print(telomere_repeats_filter)

# Define output file
svg("circos_BASEcoverage_cuml-genome.svg", width = 10, height = 10)

# Generate Circular Circos Plot
circos.clear()
circos.par("track.height" = 0.1)
circos.genomicInitialize(top_scaffolds, plotType = NULL)

# Add chromosome labels
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ycenter, 
              CELL_META$sector.index, cex = 1, niceFacing = TRUE)
}, track.height = 0.1, cell.padding = c(0, 0, 0, 0), bg.border = NA)

#add telomere repeats
circos.genomicLabels(telomere_repeats_filter, labels.column = 4, 
                     niceFacing = TRUE, cex=0.8, connection_height = mm_h(1), 
                     side = "outside")

# Add genes (blue) using coverage percentage
max_gene_coverage <- max(gene_coverage$coverage, na.rm = TRUE)
circos.genomicTrackPlotRegion(gene_coverage, ylim = c(0, max_gene_coverage * 1.1), 
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, col = "#0000FF80", lwd = 2, type = "h")
                                # Add axis labels
                                if(CELL_META$sector.index == "1") {  # Only add labels in one sector
                                  circos.yaxis(side = "left", at = seq(0, max_gene_coverage, length.out = 4), 
                                               labels = round(seq(0, max_gene_coverage, length.out = 4), 2), 
                                               tick.length = 0.02, labels.cex = 0.6)
                                }
                              })

# Add repeats (red) using coverage percentage
max_repeat_coverage <- max(repeat_coverage$coverage, na.rm = TRUE)
circos.genomicTrackPlotRegion(repeat_coverage, ylim = c(0, max_repeat_coverage * 1.1), 
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, col = "#FF000080", lwd = 2, type = "h")
                                # Add axis labels
                                if(CELL_META$sector.index == "1") {  
                                  circos.yaxis(side = "left", at = seq(0, max_repeat_coverage, length.out = 4), 
                                               labels = round(seq(0, max_repeat_coverage, length.out = 4), 2), 
                                               tick.length = 0.02, labels.cex = 0.6)
                                }
                              })
# Save and close the plot
dev.off()

#Zoom in Scaffolds where MAT regions were found.
library(karyoploteR)

# Filter scaffolds 13 and 20 from `scaffold_data`
selected_scaffolds <- scaffold_data %>% filter(Chrom %in% c("13", "20"))

# Convert scaffold data into a GRanges object
custom.genome <- GRanges(seqnames = selected_scaffolds$Chrom, 
                         ranges = IRanges(start = selected_scaffolds$start, 
                                          end = selected_scaffolds$stop))

#Generate Linear Karyotype Plot
#Sets up karyotype plot layout (pp).
#Plots scaffold ideograms (plotKaryotype).
#Adds repeat and gene densities (kpPlotRibbon).
pp <- getDefaultPlotParams(plot.type = 2) #change to 4
pp$data1height <- 150
pp$data2outmargin <- 10
pp$data1outmargin <- 10
pp$data1inmargin <- 2
pp$leftmargin <- 0.2
pp$topmargin <- 100
pp$bottommargin <- 20
pp$ideogramheight <- 2

# Initialize karyoploteR plot for scaffold 13
kp <- plotKaryotype(genome = custom.genome, chromosomes = "13", 
                    zoom=toGRanges("13:1-2000000"), 
                    plot.params = pp, plot.type = 2)

kpAddBaseNumbers(kp, tick.dist=200000, tick.len =8, tick.col="black", 
                 add.units = TRUE, cex = 0.8) 

# Add Y-axis for genes (blue)
kpAxis(kp, ymin = 0, ymax = max_gene_coverage, side = 2, r0 = 0.5, r1 = 1, cex = 0.8)
kpAddLabels(kp, labels = "Genes", r0 = 0.5, r1 = 1, label.margin = 0.05, cex = 0.8)

# Add Y-axis for repeats (red)
kpAxis(kp, ymin = 0, ymax = max_repeat_coverage, side = 2, r0 = 0.5, r1 = 0, cex = 0.8)
kpAddLabels(kp, labels = "Repeats", r0 = 0, r1 = 0.5, label.margin = 0.05, cex = 0.8)

# Plot gene coverage as a shaded area (blue)
kpPlotRibbon(kp, chr = gene_coverage$Chrom, 
             x0 = gene_coverage$start, x1 = gene_coverage$stop, 
             y1 = gene_coverage$coverage, y0 = 0, 
             col = "#0000FF80", border = "blue", r0 = 0.5, r1 = 1)

# Plot repeat coverage as a shaded area (red)
kpPlotRibbon(kp, chr = repeat_coverage$Chrom, 
             x0 = repeat_coverage$start, x1 = repeat_coverage$stop, 
             y1 = repeat_coverage$coverage, y0 = 0, 
             col = "#FF000080", border = "red", r0 = 0.5, r1 = 0)

#mark MAT regions
kpRect(kp, chr="13", x0=259771, x1=407776, y0=0, y1=1, col=NA, data.panel="1", border="black")


# Initialize karyoploteR plot for scaffold 20
kp <- plotKaryotype(genome = custom.genome, chromosomes = "20", 
                    zoom=toGRanges("20:1-2000000"), 
                    plot.params = pp, plot.type = 2)

kpAddBaseNumbers(kp, tick.dist=200000, tick.len =8, tick.col="black", 
                 add.units = TRUE, cex = 0.8) 

# Add Y-axis for genes (blue)
kpAxis(kp, ymin = 0, ymax = max_gene_coverage, side = 2, r0 = 0.5, r1 = 1, cex = 0.8)
kpAddLabels(kp, labels = "Genes", r0 = 0.5, r1 = 1, label.margin = 0.05, cex = 0.8)

# Add Y-axis for repeats (red)
kpAxis(kp, ymin = 0, ymax = max_repeat_coverage, side = 2, r0 = 0.5, r1 = 0, cex = 0.8)
kpAddLabels(kp, labels = "Repeats", r0 = 0, r1 = 0.5, label.margin = 0.05, cex = 0.8)

# Plot gene coverage as a shaded area (blue)
kpPlotRibbon(kp, chr = gene_coverage$Chrom, 
             x0 = gene_coverage$start, x1 = gene_coverage$stop, 
             y1 = gene_coverage$coverage, y0 = 0, 
             col = "#0000FF80", border = "blue", r0 = 0.5, r1 = 1)

# Plot repeat coverage as a shaded area (red)
kpPlotRibbon(kp, chr = repeat_coverage$Chrom, 
             x0 = repeat_coverage$start, x1 = repeat_coverage$stop, 
             y1 = repeat_coverage$coverage, y0 = 0, 
             col = "#FF000080", border = "red", r0 = 0.5, r1 = 0)

#mark MAT regions
kpRect(kp, chr="20", x0=1689356, x1=1743567, y0=0, y1=1, col=NA, data.panel="1", border="black")


#kpPlotDensity(kp, data=genes_kp, col = "#0000FF80", r0 = 0.5, r1 = 1, window.size = 10000)
#kpAddLabels(kp, "Genes", r0=0.5, r1=0.7, label.margin = 0.05)
#kpPlotDensity(kp, data=repeats_kp, col = "#FF000080", r0 = 0.5, r1 =0, window.size = 10000)
#kpAddLabels(kp, "Repeats", r0=0.4, r1=0.2, label.margin = 0.05)

#scaffold 20
#zoom plot scaffolds as bars
kp <- plotKaryotype(genome = custom.genome, chromosomes = "20", zoom=toGRanges("20:1400000-2200000"), plot.params = pp, plot.type = 2)
kpAddBaseNumbers(kp, tick.dist=50000, tick.len =10, tick.col="black", add.units = TRUE, cex = 0.7) 

#mark MAT regions
kpRect(kp, chr="13", x0=259771, x1=407776, y0=0, y1=1, col=NA, data.panel="1", border="black")
kpRect(kp, chr="20", x0=1689356, x1=1743567, y0=0, y1=1, col=NA, data.panel="1", border="black")


kpPlotDensity(kp, data=genes_kp, col = "#0000FF80", r0 = 0.5, r1 = 1, window.size = 1000)
kpAddLabels(kp, "Genes", r0=0.5, r1=0.7, label.margin = 0.05)
kpPlotDensity(kp, data=repeats_kp, col = "#FF000080", r0 = 0.5, r1 =0, window.size = 1000)
kpAddLabels(kp, "Repeats", r0=0.4, r1=0.2, label.margin = 0.05)

kp <- plotKaryotype(genome=custom.genome, chromosomes="all")

genes_cv <- toGRanges(gene_coverage)
kpPlotCoverage(kp, data=gene_coverage)


