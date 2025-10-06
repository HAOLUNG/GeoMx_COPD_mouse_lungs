library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(dplyr)
library(ggforce)
library(ggplot2)
library(knitr)
library(scales) 
library(cowplot) 
library(reshape2) 
library(umap)
library(Rtsne)
library(pheatmap) 
library(ggrepel)
library(tidyr)

# load obj
load("/Users/99171821/Desktop/Nanostring_timecourse/nanostring_object.Rdata")

# load between-slide results for volcano plots
Smoke2wks_vs_Air2wks <- read.csv('/Users/99171821/Desktop/Nanostring_timecourse/between_slide_DEG/results_Smoke2wks_vs_Air2wks.csv')
Smoke12wks_vs_Air12wks <- read.csv('/Users/99171821/Desktop/Nanostring_timecourse/between_slide_DEG/results_Smoke12wks_vs_Air12wks.csv')
Smoke12wks_vs_Smoke4wks <- read.csv('/Users/99171821/Desktop/Nanostring_timecourse/between_slide_DEG/results_Smoke12wks_vs_Smoke4wks.csv')

# make labels for volcano plots
Smoke2wks_vs_Air2wks$label <- NA
Smoke2wks_vs_Air2wks$label[abs(Smoke2wks_vs_Air2wks$Estimate) > 1 & Smoke2wks_vs_Air2wks$Pr...t.. < 0.05] <- Smoke2wks_vs_Air2wks$Gene[abs(Smoke2wks_vs_Air2wks$Estimate) > 1 & Smoke2wks_vs_Air2wks$Pr...t.. < 0.05]
levels(factor(Smoke2wks_vs_Air2wks$Subset))

Smoke12wks_vs_Air12wks$label <- NA
Smoke12wks_vs_Air12wks$label[abs(Smoke12wks_vs_Air12wks$Estimate) > 1 & Smoke12wks_vs_Air12wks$Pr...t.. < 0.05] <- Smoke12wks_vs_Air12wks$Gene[abs(Smoke12wks_vs_Air12wks$Estimate) > 1 & Smoke12wks_vs_Air12wks$Pr...t.. < 0.05]
levels(factor(Smoke12wks_vs_Air12wks$Subset))

Smoke12wks_vs_Smoke4wks$label <- NA
Smoke12wks_vs_Smoke4wks$label[abs(Smoke12wks_vs_Smoke4wks$Estimate) > 1 & Smoke12wks_vs_Smoke4wks$Pr...t.. < 0.05] <- Smoke12wks_vs_Smoke4wks$Gene[abs(Smoke12wks_vs_Smoke4wks$Estimate) > 1 & Smoke12wks_vs_Smoke4wks$Pr...t.. < 0.05]
levels(factor(Smoke12wks_vs_Smoke4wks$Subset))

# subset to ROI
Airway_Smoke2wks_vs_Air2wks <- Smoke2wks_vs_Air2wks %>% 
  filter(Subset == 'Airway')
Artery_Smoke2wks_vs_Air2wks <- Smoke2wks_vs_Air2wks %>% 
  filter(Subset == 'Artery')
ParencFar_Smoke2wks_vs_Air2wks <- Smoke2wks_vs_Air2wks %>% 
  filter(Subset == 'ParencFar')
ParencNear_Smoke2wks_vs_Air2wks <- Smoke2wks_vs_Air2wks %>% 
  filter(Subset == 'ParencNear')

Airway_Smoke12wks_vs_Air12wks <- Smoke12wks_vs_Air12wks %>% 
  filter(Subset == 'Airway')
Artery_Smoke12wks_vs_Air12wks <- Smoke12wks_vs_Air12wks %>% 
  filter(Subset == 'Artery')
ParencFar_Smoke12wks_vs_Air12wks <- Smoke12wks_vs_Air12wks %>% 
  filter(Subset == 'ParencFar')
ParencNear_Smoke12wks_vs_Air12wks <- Smoke12wks_vs_Air12wks %>% 
  filter(Subset == 'ParencNear')

Airway_Smoke12wks_vs_Smoke4wks <- Smoke12wks_vs_Smoke4wks %>% 
  filter(Subset == 'Airway')
Artery_Smoke12wks_vs_Smoke4wks <- Smoke12wks_vs_Smoke4wks %>% 
  filter(Subset == 'Artery')
ParencFar_Smoke12wks_vs_Smoke4wks <- Smoke12wks_vs_Smoke4wks %>% 
  filter(Subset == 'ParencFar')
ParencNear_Smoke12wks_vs_Smoke4wks <- Smoke12wks_vs_Smoke4wks %>% 
  filter(Subset == 'ParencNear')
LymphFoll_Smoke12wks_vs_Smoke4wks <- Smoke12wks_vs_Smoke4wks %>% 
  filter(Subset == 'LymphFoll')

# volcano plots
p <- ggplot(ParencNear_Smoke12wks_vs_Air12wks, aes(x = Estimate, y = -log10(Pr...t..), label = label)) +
  geom_point(aes(color = ifelse(abs(Estimate) > 1 & Pr...t.. < 0.05, "significant", "non-significant")), size = 3) +
  geom_label_repel() + 
  scale_color_manual(values = c("non-significant" = "grey", "significant" = "red")) +
  theme_minimal() +
  theme_classic() +
  labs(
    title = "ParencNear -- Smoke12wks vs Air12wks",
    x = "log2 fold change",
    y = "-log10(p-value)",
    color = "Significant"
  )
p
ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/DEG_volcano_ParencNear_Smoke12wks_vs_Air12wks.png', 
       p, width = 5, height = 7, units = "in", dpi = 1000,
       bg='white')

# load within-slide results for marker gene plots
marker_for_plots <- read.csv('/Users/99171821/Desktop/Nanostring_timecourse/within_slide_DEG/results_LymphFoll_vs_ParencFar.csv')
test <- marker_for_plots %>% arrange(Estimate) %>% filter(FDR <0.05)

# scatter-violin plots for marker genes
p <- ggplot(pData(target_nanostring),
       aes(x = roi, fill = roi,
           y = as.numeric(assayDataElement(target_nanostring["Hdc", ],
                                           elt = "q_norm")))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "Hdc expression") +
  scale_y_continuous(trans = "log2") +
  theme_bw() +
  theme_classic()
p
ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/violin_marker_parencFar_Hdc.png', 
       p, width = 5, height = 4, units = "in", dpi = 1000,
       bg='white')

# violin plot for multiple marker genes
genes_of_interest <- c("Scgb3a2", "Acta2", "Cxcl13", "Ager")
genes_of_interest <- c("Tmtc1", "Cxcl17")
expression_data <- assayDataElement(target_nanostring[genes_of_interest, ], elt = "q_norm")
data <- as.data.frame(pData(target_nanostring))
data <- as.data.frame(data$roi)
colnames(data)[1] <- 'roi'
long_data <- data %>% cbind(t(expression_data)) %>% gather(key = "gene", value = "expression", -roi)
long_data$gene <- factor(long_data$gene, levels = genes_of_interest)
p <- ggplot(long_data, aes(x = roi, fill = roi, y = as.numeric(expression))) +
  geom_violin() +
  geom_jitter(width = .2, size = 0.1) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(y = "Expression", x = "ROI", title = "Marker genes per ROI") +
  scale_y_continuous(trans = "log2") +
  theme_bw() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)  # Adjust angle and alignment
  )
p
ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/violin_markers.png', 
       p, width = 5, height = 5, units = "in", dpi = 1000,
       bg='white')

