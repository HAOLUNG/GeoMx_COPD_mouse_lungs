library(SpatialDecon)
library(GeomxTools)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)
library(multcomp)

# load geomx obj
load("/Users/99171821/Desktop/Nanostring_timecourse/nanostring_object.Rdata")

# load deconv matrix
load("/Users/99171821/Desktop/Nanostring_timecourse/SpatialDecon_deconv.Rdata")

# make deconv results readable
results.matrix <- deconv$beta
results.matrix <- results.matrix/colSums(results.matrix)
results.matrix.df  <- as.data.frame(results.matrix)
results.matrix.df <- cbind(Sample = colnames(results.matrix), t(results.matrix.df))
results.matrix.df <- cbind(ROI = target_nanostring@protocolData@data$roi, results.matrix.df)
results.matrix.df <- as.data.frame(results.matrix.df)
results.matrix.df <- cbind(Timecourse = target_nanostring@protocolData@data$timecourse, results.matrix.df)
results.matrix.df <- as.data.frame(results.matrix.df)
results.matrix.df <- results.matrix.df %>%
  mutate_at(vars(4:ncol(results.matrix.df)), as.numeric)
results.matrix.df <- results.matrix.df %>%
  mutate(across(4:ncol(results.matrix.df), as.numeric)) %>% 
  rowwise() %>%
  mutate(across(4:ncol(results.matrix.df), ~ . / sum(across(4:ncol(results.matrix.df)), na.rm = TRUE))) %>%
  ungroup() # make cell composition proportion in per ROI

df.long <- melt(results.matrix.df, id.vars = c('Sample','ROI','Timecourse'), 
                variable.name = "CellType", value.name = "Count")
df.long$Count <- as.numeric(df.long$Count)

# filter to the wanted ROI (optional: with cell type)
df_Airway <- df.long %>%
  filter(ROI == "Airway" & CellType == "B.cells")
area_Airway <- target_nanostring@phenoData@data %>%
  filter(roi == "Airway") %>% dplyr::select(area, timecourse)
df_Airway$Count <- df_Airway$Count / area_Airway$area # normalise to the area of each ROI (optional)
write.csv(df_Airway, file = '/Users/99171821/Downloads/df_airway.csv')

levels(x=as.factor(df_Airway$Timecourse))
df_Airway %>% filter(Timecourse == "Smoke_12wks") %>% dplyr::select(Count)

# stats
anova_Airway <- aov(Count ~ Timecourse, data = df_Airway)
summary(anova_Airway)
tukey_Airway <- TukeyHSD(anova_Airway)
print(tukey_Airway)

