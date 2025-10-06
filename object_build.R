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

# assemble object
dcc <- dir(file.path('/Users/99171821/Desktop/Nanostring_timecourse/DCC_files'), 
           pattern = ".dcc$", 
           full.names = TRUE, recursive = TRUE)
pkc <- c('/Users/99171821/Desktop/Nanostring_timecourse/Mm_R_NGS_WTA_v1.0.pkc')
sample_annot <- c('/Users/99171821/Desktop/Nanostring_timecourse/DSP_annot.xlsx')
nanostring <- readNanoStringGeoMxSet(dccFiles = dcc,
                                     pkcFiles = pkc,
                                     phenoDataFile = sample_annot,
                                     phenoDataSheet = "annot",
                                     phenoDataDccColName = "Sample_ID",
                                     protocolDataColNames = c("aoi", "roi"),
                                     experimentDataColNames = c("panel"))
pkcs <- annotation(nanostring)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

# correct metadata for sample and ROI ID
levels(x=as.factor(nanostring@protocolData@data$roi))

nanostring@protocolData@data$ID <- gsub('21-3025', 'Air_2wks', nanostring@protocolData@data$roi)
nanostring@protocolData@data$ID <- gsub('21-3026', 'Air_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3027', 'Air_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3028', 'Air_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3029', 'Air_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3030', 'Air_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3031', 'Air_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3032', 'Air_2wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('21-3033', 'Smoke_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3034', 'Smoke_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3035', 'Smoke_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3036', 'Smoke_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3037', 'Smoke_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3038', 'Smoke_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3039', 'Smoke_2wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('21-3040', 'Air_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3041', 'Air_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3042', 'Air_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3043', 'Air_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3044', 'Air_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3045', 'Air_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3046', 'Air_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3047', 'Air_4wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('21-3048', 'Smoke_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3049', 'Smoke_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3050', 'Smoke_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3051', 'Smoke_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3052', 'Smoke_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3053', 'Smoke_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3054', 'Smoke_4wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21-3055', 'Smoke_4wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('22-0167', 'Air_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0168', 'Air_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0169', 'Air_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0170', 'Air_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0171', 'Air_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0172', 'Air_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0173', 'Air_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0174', 'Air_6wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('22-0175', 'Smoke_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0176', 'Smoke_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0177', 'Smoke_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0178', 'Smoke_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0179', 'Smoke_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0180', 'Smoke_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0181', 'Smoke_6wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0182', 'Smoke_6wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('22-0183', 'Air_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0184', 'Air_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0185', 'Air_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0186', 'Air_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0187', 'Air_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0188', 'Air_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0189', 'Air_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0190', 'Air_8wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('22-0191', 'Smoke_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0192', 'Smoke_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0193', 'Smoke_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0194', 'Smoke_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0195', 'Smoke_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0196', 'Smoke_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0197', 'Smoke_8wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0198', 'Smoke_8wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('22-0199', 'Air_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0200', 'Air_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0201', 'Air_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0202', 'Air_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0203', 'Air_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0204', 'Air_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0205', 'Air_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0206', 'Air_12wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('22-0207', 'Smoke_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0208', 'Smoke_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0209', 'Smoke_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0210', 'Smoke_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0211', 'Smoke_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0212', 'Smoke_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0213', 'Smoke_12wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('22-0214', 'Smoke_12wks', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub('\\b001\\b', 'Random', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('\\b002\\b', 'Random', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('\\b22-3\\b', 'Random', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('001Smoke', 'Smoke', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('0Air', 'Air', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('21_3028', 'Air_2wks', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('004.', '004', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('001]', '001', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('001 Random', '001', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('002.', '002', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub('002Random', '002', nanostring@protocolData@data$ID)
nanostring@protocolData@data$ID <- gsub("001'", '001', nanostring@protocolData@data$ID)

nanostring@protocolData@data$ID <- gsub("Airway", 'Airway', nanostring@protocolData@data$ID, ignore.case = TRUE)

nanostring@protocolData@data$ID <- gsub("Artery", 'Artery', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("Airtery", 'Artery', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("Arteryt", 'Artery', nanostring@protocolData@data$ID, ignore.case = TRUE)

nanostring@protocolData@data$ID <- gsub("ParencFar", 'ParencFar', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("ParenFar", 'ParencFar', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("ParecFar", 'ParencFar', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("PaencFar", 'ParencFar', nanostring@protocolData@data$ID, ignore.case = TRUE)

nanostring@protocolData@data$ID <- gsub("ParencNear", 'ParencNear', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("ParenNear", 'ParencNear', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("PrencNear", 'ParencNear', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("ParencfNear", 'ParencNear', nanostring@protocolData@data$ID, ignore.case = TRUE)
nanostring@protocolData@data$ID <- gsub("ParecNear", 'ParencNear', nanostring@protocolData@data$ID, ignore.case = TRUE)

nanostring@protocolData@data$ID <- gsub("_LF_", '_LymphFoll_', nanostring@protocolData@data$ID, ignore.case = TRUE)

ID_check <- levels(x=as.factor(nanostring@protocolData@data$ID))

# correct for slide name
nanostring@protocolData@data$slide <- gsub("([0-9]+([-_])[0-9]+).*", "\\1", nanostring@protocolData@data$roi)
nanostring@protocolData@data$slide <- gsub("-", "_", nanostring@protocolData@data$slide)
nanostring@protocolData@data$slide <- gsub("00122_0213", "22_0213", nanostring@protocolData@data$slide)
nanostring@protocolData@data$slide <- gsub("022_0201", "22_0201", nanostring@protocolData@data$slide)
slide_check <- levels(x=as.factor(nanostring@protocolData@data$slide))

# extract ID info from ID
nanostring@protocolData@data$roi <- sub(".*_(\\d+wks_)([^_]+)_\\d+", "\\2", nanostring@protocolData@data$ID)
levels(x=as.factor(nanostring@protocolData@data$roi))

# extract timecourse info from ID
nanostring@protocolData@data$timecourse <- sub("(_\\d+wks).*", "\\1", nanostring@protocolData@data$ID)
levels(x=as.factor(nanostring@protocolData@data$timecourse))

# QC
nanostring <- shiftCountsOne(nanostring, useDALogic = TRUE)

# segment QC
QC_params <-
  list(minSegmentReads = 1000, 
       percentTrimmed = 80,    
       percentStitched = 80,   
       percentAligned = 80,    
       percentSaturation = 50,
       minNegativeCount = 1,   
       maxNTCCount = 1000,    
       minNuclei = 20,         
       minArea = 5000)      
nanostring <- setSegmentQCFlags(nanostring, qcCutoffs = QC_params)        

# collate QC results
QCResults <- protocolData(nanostring)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# graphical summaries of QC statistics plot function
col_by <- 'roi'
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}
p<-QC_histogram(sData(nanostring), "Trimmed (%)", col_by, 80)
p<-QC_histogram(sData(nanostring), "Stitched (%)", col_by, 80)
p<-QC_histogram(sData(nanostring), "Aligned (%)", col_by, 80)
p<-QC_histogram(sData(nanostring), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
p<-QC_histogram(sData(nanostring), "area", col_by, 5000, scale_trans = "log10")

ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/QC_saturated.png', 
       p, width = 3, height = 5, units = "in", dpi = 1000,
       bg='white')

# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(nanostring), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(nanostring)[["NegGeoMean"]] <- negativeGeoMeans
negCols <- paste0("NegGeoMean_", modules)
pData(nanostring)[, negCols] <- sData(nanostring)[["NegGeoMean"]]

# detatch neg_geomean columns ahead of aggregateCounts call
pData(nanostring) <- pData(nanostring)[, !colnames(pData(nanostring)) %in% negCols]
kable(table(NTC_Count = sData(nanostring)$NTC),
      col.names = c("NTC Count", "# of Segments"))

# plot all QC metrics
kable(QC_Summary, caption = "QC Summary Table for each Segment")
dim(nanostring)

# remove flagged ID
nanostring <- nanostring[, QCResults$QCStatus == "PASS"]
dim(nanostring)

# probe QC
nanostring <- setBioProbeQCFlags(nanostring, 
                                 qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                                 removeLocalOutliers = TRUE)
ProbeQCResults <- fData(nanostring)[["QCFlags"]]
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df # the number of global and local outlier probes

# exclude outlier probes
ProbeQCPassed <- subset(nanostring, 
                        fData(nanostring)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
                          fData(nanostring)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
dim(nanostring)
nanostring <- ProbeQCPassed 

# create gene-level count data
length(unique(featureData(nanostring)[["TargetName"]]))
target_nanostring <- aggregateCounts(nanostring) # collapse to targets
dim(target_nanostring)
exprs(target_nanostring)[1:5, 1:2]

# calculate limit of quantification (LOQ) per segment
cutoff <- 2
minLOQ <- 2
LOQ <- data.frame(row.names = colnames(target_nanostring))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_nanostring)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_nanostring)[, vars[1]] * 
             pData(target_nanostring)[, vars[2]] ^ cutoff)
  }
}
pData(target_nanostring)$LOQ <- LOQ

# determine the number of genes detected per segment
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_nanostring)$Module == module
  Mat_i <- t(esApply(target_nanostring[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
LOQ_Mat <- LOQ_Mat[fData(target_nanostring)$TargetName, ]

# filter out segments with very low signal
pData(target_nanostring)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_nanostring)$GeneDetectionRate <-
  pData(target_nanostring)$GenesDetected / nrow(target_nanostring)
pData(target_nanostring)$DetectionThreshold <- 
  cut(pData(target_nanostring)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
p<-ggplot(pData(target_nanostring),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = target_nanostring@protocolData@data$roi)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/QC_segment_gene_detect.png', 
       p, width = 4, height = 3, units = "in", dpi = 1000,
       bg='white')

dim(target_nanostring)
target_nanostring <- target_nanostring[, pData(target_nanostring)$GeneDetectionRate >= 0.1] # filter out segments with gene detetion <10%
dim(target_nanostring)

# gene detection rate
LOQ_Mat <- LOQ_Mat[, colnames(target_nanostring)]
fData(target_nanostring)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_nanostring)$DetectionRate <-
  fData(target_nanostring)$DetectedSegments / nrow(pData(target_nanostring))
goi <- c('Epcam','Scgb3a2','Pecam1','Ager','Sftpc','Sftpa1','Msln','Col1a1','Siglecf') # gene of interest detection table
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_nanostring)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_nanostring)[goi, "DetectionRate"]))
goi_df

# gene filtering
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <- unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                                    function(x) {sum(fData(target_nanostring)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_nanostring))
rownames(plot_detect) <- plot_detect$Freq
ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")
negativeProbefData <- subset(fData(target_nanostring), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
dim(target_nanostring)
target_nanostring <- target_nanostring[fData(target_nanostring)$DetectionRate >= 0.1 |
                                         fData(target_nanostring)$TargetName %in% neg_probes, ] # focus on the genes detected in >=10% of segments
dim(target_nanostring)
goi <- goi[goi %in% rownames(target_nanostring)] # retain only detected genes of interest

# normalisation
ann_of_interest <- "roi"
Stat_data <- data.frame(row.names = colnames(exprs(target_nanostring)),
                        Segment = colnames(exprs(target_nanostring)),
                        Annotation = target_nanostring@protocolData@data[, ann_of_interest],
                        Q3 = unlist(apply(exprs(target_nanostring), 2, quantile, 0.75, na.rm = TRUE)),
                        NegProbe = exprs(target_nanostring)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")
plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")
plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
p<-plot_grid(plt1, btm_row, ncol = 1, labels = c("A", "")) # good separation between Q3 normalised data and neg probe counts
ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/QC_Q3_NegProbe_separation.png', 
       p, width = 8, height = 5, units = "in", dpi = 1000,
       bg='white')

# Q3 norm (75th percentile)
target_nanostring <- normalize(target_nanostring, 
                               norm_method = "quant", 
                               desiredQuantile = 0.75,
                               toElt = "q_norm")

# background normalization 
target_nanostring <- normalize(target_nanostring,
                               norm_method = "neg", 
                               fromElt = "exprs",
                               toElt = "neg_norm")

# visualise the first 10 segments with each normalisation method
par(mfrow = c(1, 3))

p1<-boxplot(exprs(target_nanostring)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")
p2<-boxplot(assayDataElement(target_nanostring[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")
p3<-boxplot(assayDataElement(target_nanostring[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")

# UMAP
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
umap_out <- umap(t(log2(assayDataElement(target_nanostring , elt = "q_norm"))),  
                 config = custom_umap)
pData(target_nanostring)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
target_nanostring@protocolData@data$roi <- gsub("ParencFar|ParencNear", "Parenchyma", 
                                                target_nanostring@protocolData@data$roi) # integrate both parenchyma naming
p<-ggplot(pData(target_nanostring),
       aes(x = UMAP1, y = UMAP2, color = target_nanostring@protocolData@data$roi)) +
  theme_bw() + geom_point(size = 2) + theme(legend.title = element_blank()) # good separation between ROI
p<-ggplot(pData(target_nanostring),
       aes(x = UMAP1, y = UMAP2, color = target_nanostring@protocolData@data$timecourse)) +
  theme_bw() + geom_point(size = 2) + theme(legend.title = element_blank())
p<-ggplot(pData(target_nanostring),
          aes(x = UMAP1, y = UMAP2, color = target_nanostring@protocolData@data$slide)) +
  theme_bw() + geom_point(size = 2) + theme(legend.title = element_blank())

ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/umap_int.png', 
       p, width = 5, height = 3, units = "in", dpi = 1000,
       bg='white')

# tSNE
set.seed(42)
tsne_out <- Rtsne(t(log2(assayDataElement(target_nanostring , elt = "q_norm"))),
                  perplexity = ncol(target_nanostring)*.15)
pData(target_nanostring)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
ggplot(pData(target_nanostring),
       aes(x = tSNE1, y = tSNE2, color = target_nanostring@protocolData@data$roi)) +
  geom_point(size = 3) +
  theme_bw() # also good separation between IDs

# PCA
pca_result <- prcomp(t(log2(assayDataElement(target_nanostring, elt = "q_norm"))))
pca_scores <- as.data.frame(pca_result$x)
pData(target_nanostring)[, c("PC1", "PC2")] <- pca_scores[, 1:2]
custom_colors <- c("Airway" = "red", "Artery" = "green", "LymphFoll" = "orange", 
                   "ParencFar" = "blue", "ParencNear" = "purple")
pData(target_nanostring)$roi <- target_nanostring@protocolData@data$roi
p <- ggplot(pData(target_nanostring),
            aes(x = PC1, y = PC2, color = roi)) +
  geom_point(size = 2) +
  scale_color_manual(values = custom_colors) +  # Specify custom colors
  theme_bw() +
  theme(legend.title = element_blank())
p
ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/pca_roi.png', 
       p, width = 5, height = 3, units = "in", dpi = 1000,
       bg='white')

# clustering high coefficient of variation (CV) genes
assayDataElement(object = target_nanostring, elt = "log_q") <- 
  assayDataApply(target_nanostring, 2, FUN = log, base = 2, elt = "q_norm") # create log2 transform of the data for analysis
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_nanostring,
                         elt = "log_q", MARGIN = 1, calc_CV)
sort(CV_dat, decreasing = TRUE)[1:5] # genes with highest CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
pData(target_nanostring)[, c("ROI")] <- target_nanostring@protocolData@data$roi
pData(target_nanostring)[, c("Timecourse")] <- target_nanostring@protocolData@data$timecourse
p<-pheatmap(assayDataElement(target_nanostring[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_nanostring)[, c("ROI", "Timecourse")])
  
p
ggsave('/Users/99171821/Desktop/Nanostring_timecourse/plots/heatmap_ROI_timecourse.png', 
       p, width = 5, height = 8, units = "in", dpi = 1000,
       bg='white')

# save & load
save(target_nanostring, file = "/Users/99171821/Desktop/Nanostring_timecourse/nanostring_object.Rdata")
load("/Users/99171821/Desktop/Nanostring_timecourse/nanostring_object.Rdata")

# within-slide DE
pData(target_nanostring)$roi <- factor(pData(target_nanostring)$roi)
pData(target_nanostring)$testRegion <- factor(pData(target_nanostring)$roi)
pData(target_nanostring)[["slide"]] <- factor(pData(target_nanostring)[["slide name"]])
assayDataElement(object = target_nanostring, elt = "log_q") <-
  assayDataApply(target_nanostring, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM
timecourse <- unique(target_nanostring@protocolData@data$timecourse)
results <- c()
for(time in timecourse) {
  ind <- pData(target_nanostring)$timecourse == time
  mixedOutmc <-
    mixedModelDE(target_nanostring[, ind],
                 elt = "log_q",
                 modelFormula = ~ testRegion + (1 + testRegion | slide),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- time
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}















