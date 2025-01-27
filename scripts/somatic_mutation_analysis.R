# This script runs the Somatic Mutational Analysis of the thesis
### Libraries
library(maftools)

### Load data
maf <- read.csv("../raw_data/mutations.csv")
maf <- read.maf(maf)
TF <- unique(TRUSST_data$TF)

maf_tf <- subsetMaf(maf = maf, isTCGA = TRUE, genes = TF)

## MAF Summary & Oncoplot
png("../figures/maf_summary.png")
plotmafSummary(maf = maf_tf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, 
               log_scale = FALSE, top = 10)
dev.off()

png("../figures/maf_oncoplot.png")
oncoplot(maf = maf_tf, top = 15)
dev.off()


### Lollipop Plots ----
##-- TP53
png("../figures/loll_tp53.png")
lollipopPlot(
  maf,
  gene = 'TP53',
  labelPos = c(157, 126, 234, 257),
  labPosSize = 0.6,
  showDomainLabel = FALSE
  )
dev.off()

##-- STAT3
png("../figures/loll_stat3.png")
lollipopPlot(
  maf,
  gene = 'STAT3',
  labelPos = 564,
  labPosSize = 0.8,
  showDomainLabel = FALSE
)
dev.off()

##-- RUNX1
png("../figures/loll_runx1.png")
lollipopPlot(
  maf,
  gene = 'RUNX1',
  labelPos = c(76,141),
  labPosSize = 0.8,
  showDomainLabel = FALSE
)
dev.off()

##-- FOXA1
png("../figures/loll_foxa1.png")
lollipopPlot(
  maf,
  gene = 'FOXA1',
  labelPos = 214,
  labPosSize = 0.8,
  showDomainLabel = FALSE
)
dev.off()

##-- ATM
png("../figures/loll_atm.png")
lollipopPlot(
  maf,
  gene = 'ATM',
  labelPos = 2780,
  labPosSize = 0.8,
  showDomainLabel = FALSE
)
dev.off()

##-- NF1
png("../figures/loll_nf1.png")
lollipopPlot(
  maf,
  gene = 'NF1',
  AACol = NULL,
  labelPos = 529,
  labPosSize = 0.8,
  showDomainLabel = FALSE
  )
  dev.off()

##-- SMARCA4
png("../figures/loll_smarca4.png")
lollipopPlot(
  maf,
  gene = 'SMARCA4',
  AACol = NULL,
  labelPos = 1500,
  labPosSize = 0.8,
  showDomainLabel = FALSE
  )
  dev.off()
