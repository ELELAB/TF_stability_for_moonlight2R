### Libraries
library(maftools)

### Load data
maf <- read.csv("../raw_data/mutations.csv")
maf <- read.maf(maf)
TF <- unique(TRUSST_data$TF)

maf_tf <- subsetMaf(maf = maf, isTCGA = TRUE, genes = TF)

## Summary statistics
png("../figures/maf_summary.png")
plotmafSummary(maf = maf_tf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, 
               log_scale = FALSE, top = 10)
dev.off()

png("../figures/maf_oncoplot.png")
oncoplot(maf = maf_tf, top = 15)
dev.off()

oncoplot(maf = maf, top = 15)
## TP53
png("../figures/loll_tp53.png")
lollipopPlot(
  maf,
  data = NULL,
  gene = 'TP53',
  AACol = NULL,
  labelPos = c(251, 135, 245, 274, 111,157, 158, 126, 234, 257),
  labPosSize = 0.3,
  showMutationRate = TRUE,
  showDomainLabel = TRUE,
  cBioPortal = FALSE,
  refSeqID = NULL,
  proteinID = NULL,
  roundedRect = TRUE,
  repel = TRUE,
  collapsePosLabel = TRUE,
  showLegend = TRUE,
  legendTxtSize = 0.8,
  labPosAngle = 45,
  domainLabelSize = 0.8,
  axisTextSize = c(1, 1),
  printCount = TRUE,
  colors = NULL,
  domainAlpha = 1,
  domainBorderCol = "black",
  bgBorderCol = "black",
  labelOnlyUniqueDoamins = FALSE,
  defaultYaxis = FALSE,
  titleSize = c(1.2, 1),
  pointSize = 1.5
)
dev.off()

#########
png("../figures/loll_stat3.png")
lollipopPlot(
  maf,
  data = NULL,
  gene = 'STAT3',
  AACol = NULL,
  labelPos = 564,
  labPosSize = 0.5,
  showMutationRate = TRUE,
  showDomainLabel = TRUE,
  cBioPortal = FALSE,
  refSeqID = NULL,
  proteinID = NULL,
  roundedRect = TRUE,
  repel = TRUE,
  collapsePosLabel = TRUE,
  showLegend = TRUE,
  legendTxtSize = 0.8,
  labPosAngle = 45,
  domainLabelSize = 0.8,
  axisTextSize = c(1, 1),
  printCount = TRUE,
  colors = NULL,
  domainAlpha = 1,
  domainBorderCol = "black",
  bgBorderCol = "black",
  labelOnlyUniqueDoamins = FALSE,
  defaultYaxis = FALSE,
  titleSize = c(1.2, 1),
  pointSize = 1.5
)
dev.off()
######
png("../figures/loll_runx1.png")
lollipopPlot(
  maf,
  data = NULL,
  gene = 'RUNX1',
  AACol = NULL,
  labelPos = 141,
  labPosSize = 0.5,
  showMutationRate = TRUE,
  showDomainLabel = TRUE,
  cBioPortal = FALSE,
  refSeqID = NULL,
  proteinID = NULL,
  roundedRect = TRUE,
  repel = TRUE,
  collapsePosLabel = TRUE,
  showLegend = TRUE,
  legendTxtSize = 0.8,
  labPosAngle = 45,
  domainLabelSize = 0.8,
  axisTextSize = c(1, 1),
  printCount = TRUE,
  colors = NULL,
  domainAlpha = 1,
  domainBorderCol = "black",
  bgBorderCol = "black",
  labelOnlyUniqueDoamins = FALSE,
  defaultYaxis = FALSE,
  titleSize = c(1.2, 1),
  pointSize = 1.5
)
dev.off()
######
png("../figures/loll_foxa1.png")
lollipopPlot(
  maf,
  data = NULL,
  gene = 'FOXA1',
  AACol = NULL,
  labelPos = 214,
  labPosSize = 0.5,
  showMutationRate = TRUE,
  showDomainLabel = TRUE,
  cBioPortal = FALSE,
  refSeqID = NULL,
  proteinID = NULL,
  roundedRect = TRUE,
  repel = TRUE,
  collapsePosLabel = TRUE,
  showLegend = TRUE,
  legendTxtSize = 0.8,
  labPosAngle = 45,
  domainLabelSize = 0.8,
  axisTextSize = c(1, 1),
  printCount = TRUE,
  colors = NULL,
  domainAlpha = 1,
  domainBorderCol = "black",
  bgBorderCol = "black",
  labelOnlyUniqueDoamins = FALSE,
  defaultYaxis = FALSE,
  titleSize = c(1.2, 1),
  pointSize = 1.5
)
dev.off()
