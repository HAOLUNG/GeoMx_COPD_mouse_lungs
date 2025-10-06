library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# load smoke timecourse lung scRNA-seq data
lung_timecourse_scRNAseq <- readRDS('/Users/99171821/Desktop/Nanostring_timecourse/lung_timecourse_smoking/mouse_COPD.rds')
Idents(lung_timecourse_scRNAseq) <- lung_timecourse_scRNAseq$cluster_labels_res.0.3
DimPlot(lung_timecourse_scRNAseq, label = T) + NoLegend()

# load lymphoid follicle marker geneset
LF_marker <- read.csv('/Users/99171821/Desktop/Nanostring_timecourse/LymphFoll_marker_GOOD.csv')
LF_marker <- LF_marker[,2]
LF_marker <- list(LF_marker)

# Gene set testing in scRNAseq
lung_timecourse_scRNAseq <- AddModuleScore(object = lung_timecourse_scRNAseq, 
                                           features = LF_marker, 
                                           name = "LymphFoll")
colnames(lung_timecourse_scRNAseq@meta.data)[47] <- 'LymphFoll'
FeaturePlot(lung_timecourse_scRNAseq, 
            features = 'LymphFoll',
            split.by = 'Timepoint')
DotPlot(lung_timecourse_scRNAseq, 
        features = 'LymphFoll',
        split.by = 'Timepoint',
        cols=c(rep("blue",6), "white"))

# subset to B cells in scRNA-seq data
B_cell_lung_timecourse <- subset(x = lung_timecourse_scRNAseq, idents = c("B cells", "B cells_2"))

# Test LymphFoll gene set on B cells
B_cell_lung_timecourse <- AddModuleScore(object = B_cell_lung_timecourse, 
                                         features = LF_marker, 
                                         name = "LymphFoll")
colnames(B_cell_lung_timecourse@meta.data)[47] <- 'LymphFoll'
FeaturePlot(B_cell_lung_timecourse, 
            features = 'LymphFoll',
            split.by = 'Timepoint')
DotPlot(B_cell_lung_timecourse, 
        features = 'LymphFoll',
        split.by = 'Timepoint',
        cols=c(rep("blue",6), "white"))

# plot cell numbers for B cells across timepoint
metadata <- lung_timecourse_scRNAseq@meta.data
metadata <- metadata %>%
  mutate(Timepoint = str_replace_all(Timepoint, c("CS" = "smoke", "Week" = "week", "Air" = 'air')))
metadata <- metadata %>%
  mutate(Timepoint = factor(Timepoint, levels = c("12 week air", 
                                                  "2 week smoke", 
                                                  "4 week smoke", 
                                                  "6 week smoke", 
                                                  "8 week smoke", 
                                                  "12 week smoke")))
filtered_metadata <- metadata %>%
  filter(cluster_labels_res.0.3 == "B cells")
cell_counts <- filtered_metadata %>%
  group_by(Timepoint, cluster_labels_res.0.3) %>%
  summarise(Count = n()) %>%
  ungroup()
cell_counts$Timepoint <- as.factor(cell_counts$Timepoint)
ggplot(cell_counts, aes(x = Timepoint, y = Count, fill = cluster_labels_res.0.3)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of B cells per timepoint",
       y = "Number of B cells") +
  theme_minimal() +
  theme(legend.position = "none") 

