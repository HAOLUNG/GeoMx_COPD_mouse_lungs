library(SpatialDecon)
library(GeomxTools)
library(ggplot2)
library(reshape2)
library(dplyr)

# load geomx obj
load("/Users/99171821/Desktop/Nanostring_timecourse/nanostring_object.Rdata")

# load MLCA sig matrix
MLCA <- read.delim('/Users/99171821/Desktop/Nanostring_timecourse/MLCA_sig_matrix/Mouse_Atlas_Signature.txt')
rownames(MLCA) <- MLCA[,1]
MLCA <- MLCA[,-1]
MLCA <- as.matrix(MLCA)
heatmap(sweep(MLCA, 1, apply(MLCA, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

# define bg of nanostring data
per.observation.mean.neg <- target_nanostring@assayData[["q_norm"]]["NegProbe-WTX", ]
bg <- sweep(target_nanostring@assayData[["q_norm"]] * 0, 2, per.observation.mean.neg, "+")

# deconv
deconv <- spatialdecon(norm = target_nanostring@assayData[["q_norm"]],          
                       bg = bg,                       
                       X = MLCA, 
                       align_genes = TRUE)     

# save & load
save(deconv, file = "/Users/99171821/Desktop/Nanostring_timecourse/SpatialDecon_deconv.Rdata")
load("/Users/99171821/Desktop/Nanostring_timecourse/SpatialDecon_deconv.Rdata")

heatmap(deconv$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))
heatmap(sweep(deconv$X, 1, apply(deconv$X, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

# plotting for QC based on cell component per roi
levels(x=factor(target_nanostring@protocolData@data$roi))
png(file="/Users/99171821/Desktop/Nanostring_timecourse/plots/SpatialDecon_airway.png",
    width = 12, height = 6, units = "in", res = 1000)
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 4))
TIL_barplot(deconv$prop_of_nontumor[,target_nanostring@protocolData@data$roi == 'Airway'], 
                 draw_legend = TRUE, cex.names = 0.75)
dev.off()
png(file="/Users/99171821/Desktop/Nanostring_timecourse/plots/SpatialDecon_artery.png",
    width = 12, height = 6, units = "in", res = 1000)
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 4))
TIL_barplot(deconv$prop_of_nontumor[,target_nanostring@protocolData@data$roi == 'Artery'], 
            draw_legend = TRUE, cex.names = 0.75)
dev.off()
png(file="/Users/99171821/Desktop/Nanostring_timecourse/plots/SpatialDecon_LymphFoll.png",
    width = 12, height = 6, units = "in", res = 1000)
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 4))
TIL_barplot(deconv$prop_of_nontumor[,target_nanostring@protocolData@data$roi == 'LymphFoll'], 
            draw_legend = TRUE, cex.names = 0.75)
dev.off()
png(file="/Users/99171821/Desktop/Nanostring_timecourse/plots/SpatialDecon_ParencFar.png",
    width = 12, height = 6, units = "in", res = 1000)
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 4))
TIL_barplot(deconv$prop_of_nontumor[,target_nanostring@protocolData@data$roi == 'ParencFar'], 
            draw_legend = TRUE, cex.names = 0.75)
dev.off()
png(file="/Users/99171821/Desktop/Nanostring_timecourse/plots/SpatialDecon_ParencNear.png",
    width = 12, height = 6, units = "in", res = 1000)
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 4))
TIL_barplot(deconv$prop_of_nontumor[,target_nanostring@protocolData@data$roi == 'ParencNear'], 
            draw_legend = TRUE, cex.names = 0.75)
dev.off()

# plotting barplot of cell deconv per roi
results.matrix <- deconv$beta
results.matrix <- results.matrix/colSums(results.matrix)
results.matrix.df  <- as.data.frame(results.matrix)
results.matrix.df <- cbind(Sample = colnames(results.matrix), t(results.matrix.df))
results.matrix.df <- cbind(ROI = target_nanostring@protocolData@data$roi, results.matrix.df)
results.matrix.df <- as.data.frame(results.matrix.df)

df.long <- melt(results.matrix.df, 
                id.vars = c('Sample','ROI'),
                variable.name = "CellType", value.name = "Count")
df.long$Count <- as.numeric(df.long$Count)
df.summary <- df.long %>%
  group_by(ROI, CellType) %>%
  summarise(TotalCount = sum(Count))
normalize_specific_column <- function(data, column_index, row_group_size) {
  data[, column_index] <- unlist(lapply(seq(1, nrow(data), by = row_group_size), function(i) {
    group <- data[i:(i + row_group_size - 1), column_index]
    group / sum(group, na.rm = TRUE)
  }))
  return(data)
}
df.summary <- normalize_specific_column(df.summary, 3, 16)

cell_type_colors <- c(
  "AT.I.epithelial" = "#1f77b4",         # Blue
  "AT.II.epithelial" = "#ff7f0e",        # Orange
  "Alveolar.macrophages" = "#2ca02c",    # Green
  "B.cells" = "#d62728",                 # Red
  "Capillary.endothelial" = "#9467bd",   # Purple
  "Cd4.T.cells" = "#8c564b",             # Brown
  "Cd8.T.cells" = "#e377c2",             # Pink
  "Dendritic.cells" = "#7f7f7f",         # Gray
  "Fibroblast" = "#bcbd22",              # Olive
  "Granulocyte" = "#17becf",             # Cyan
  "Interstitial.macrophages" = "#aec7e8",# Light Blue
  "Monocytes" = "#ffbb78",               # Light Orange
  "Multiciliated_Deuterostome" = "#98df8a",# Light Green
  "NK.cells" = "#ff9896",                # Light Red
  "Secretory" = "#c5b0d5",               # Light Purple
  "Venous.endothelial" = "#8c564b"       # Replaced with another color to avoid similarity with Cd4.T.cells
)
p <- ggplot(df.summary, aes(x = ROI, y = TotalCount, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell composition by ROI", x = "ROI", y = "Predicted cell composition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = cell_type_colors)
p
ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/SpatialDecon_barplot_avg_int.png', p, 
       width = 6, height = 5, units = "in", dpi = 1000, bg = 'white')

# plotting barplot of cell deconv for each roi separated by timepoint
results.matrix.df <- cbind(Timecourse = target_nanostring@protocolData@data$timecourse, results.matrix.df)
df.long <- melt(results.matrix.df, id.vars = c('Sample', 'ROI', 'Timecourse'), 
                variable.name = "CellType", value.name = "Count")
df.long$Count <- as.numeric(df.long$Count)
df.summary <- df.long %>%
  group_by(ROI, Timecourse, CellType) %>%
  summarise(TotalCount = sum(Count))
normalize_specific_column <- function(data, column_index, row_group_size) {
  data[, column_index] <- unlist(lapply(seq(1, nrow(data), by = row_group_size), function(i) {
    group <- data[i:(i + row_group_size - 1), column_index]
    group / sum(group, na.rm = TRUE)
  }))
  return(data)
}
df.summary <- normalize_specific_column(df.summary, 4, 16)

unique_rois <- unique(df.summary$ROI)
output_directory <- '/Users/99171821/Desktop/Nanostring_timecourse/plots/'
for (roi in unique_rois) {
  df.roi <- df.summary %>% filter(ROI == roi)
  p <- ggplot(df.roi, aes(x = Timecourse, y = TotalCount, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Cell composition over time for", roi), x = "Timecourse", y = "Predicted cell composition") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = cell_type_colors)
  
  print(p)
  ggsave(filename = paste0(output_directory, "SpatialDecon_barplot_time_", roi, "_NEW.png"), plot = p, 
         width = 6, height = 5, units = "in", dpi = 1000, bg = 'white')
}

