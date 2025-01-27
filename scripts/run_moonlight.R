# This script runs the Moonlight2R pipeline of the thesis

# SET UP -------------------------------
# Library ----------------
library(TCGAbiolinks)
library(Moonlight2R) 
library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Data -------------------
DEA_lumA <- get(load("../rawdata/BRCA_BRCA.LumA_dataDEGs.rda"))
Filt_lumA_raw <- get(load("../rawdata/BRCA_BRCA.LumA_dataFilt_HUGO.rda"))
Maf_lumA <- read.csv("../rawdata/mutations.csv")
TRUSST_data <- read_tsv("../rawdata/trrust_rawdata.human.tsv", col_names = c("TF", "Target", "InteractionType", "DOI")) %>%
  filter(InteractionType %in% c("Activation", "Repression"))

# Wrangle data -------------------------
# Removing non-tumor data points: 
tum.samples <- TCGAquery_SampleTypes(barcode = colnames(Filt_lumA_raw), 
                                     typesample = "TP")
Filt_lumA <- Filt_lumA_raw[,which(colnames(Filt_lumA_raw) %in% tum.samples)]
print('everything is loaded')

# # Call Moonlight functions -----------------------
# # FEA -------
FEA_lumA <- FEA(DEGsmatrix = DEA_lumA)
save(FEA_lumA, file = "../results/FEA_lumA.rda")
# Plot result ---------------------------
png(filename = "../results/plotFEA_lumA.png",
	width = 8, height = 4, units = "in", pointsize = 4, res = 1200)
plotFEA(FEA_lumA, topBP = 20)
dev.off()
print('FEA is finished')


# GRN -------
GRN_lumA <- GRN(TFs = rownames(DEA_lumA),
                 DEGsmatrix = DEA_lumA,
                 normCounts = Filt_lumA,
                 kNearest = 3,
                 DiffGenes = TRUE,
                 nGenesPerm = 1000,
                 nBoot = 100)

#filter to include only targets
common_targets <- intersect(TRUSST_data$Target, rownames(GRN_lumA$miTFGenes))
GRN_lumA$miTFGenes <- GRN_lumA$miTFGenes[common_targets, , drop = FALSE]
GRN_lumA$maxmi <- GRN_lumA$maxmi[common_targets] 
print('GRN is finished')
save(GRN_lumA, file = "../results/GRN_lumA.rda")



# URA ------
URA_lumA <- URA(dataGRN = GRN_lumA,
                  DEGsmatrix = DEA_lumA,
                  nCores = 4,
                  BPname = c("proliferation of cells","apoptosis"))
save(URA_lumA, file = "../results/URA_lumA.rda")
print('URA is done')

# PRA ------
PRA_lumA <- PRA(dataURA = URA_lumA,
                  BPname = c("proliferation of cells", "apoptosis"),
                  thres.role = 0)
save(PRA_lumA, file = "../results/PRA_lumA.rda")
drivers <- PRAtoTibble(PRA_lumA)
write.csv(drivers, file = "../results/PRA_drivers.csv", row.names = FALSE)
print('PRA is done')


# Load MAVISp -------
dataMAVISp <- loadMAVISp(mavispDB = "../rawdata/mavisp_csv/10012025_ALL",
                       proteins_of_interest = NULL,
                       mode = 'simple',
                       ensemble = 'md')

# TF -------
TF_lumA <- TFinfluence(dataPRA = PRA_lumA,
            dataTRRUST = TRUSST_data,
            dataMAF = Maf_lumA,
            dataDEGs = DEA_lumA,
            dataMAVISp = dataMAVISp,
            stabClassMAVISp = 'rasp'
            )
save(TF_lumA, file =  "../results/TF_lumA.rda")

# Final Table ----
final_table <- TF_lumA %>%
  filter(
    stab_class == "Destabilizing", 
    mutation_available == TRUE, 
    Moonlight_Oncogenic_Mediator %in% c("TSG", "OG"),
    InteractionType %in% c("Activation", "Repression")
  ) %>%
  group_by(TF, Moonlight_Oncogenic_Mediator, InteractionType) %>%
  summarize(
    targets = paste(unique(Target), collapse = "; "),
    destabilizing_mutations = paste(unique(tf_mutation), collapse = "; "),
    .groups = 'drop'
  )
write.csv(final_table, file = "../results/TFinfluence_summary_table.csv", row.names = FALSE)
print('TFinfluence is done')

# Plots -----
# MAVISp coverage ----
TF_lumA_classified <- TF_lumA %>%
  mutate(
    stab_class_dynamic = case_when(
      is.na(stab_class) & !in_MAVISp & !mutation_available ~ "TF not in MAVISP yet",
      is.na(stab_class) & in_MAVISp & !mutation_available ~ "mutation in IDR",
      TRUE ~ stab_class 
    )
  ) %>%
  group_by(TF, stab_class_dynamic) %>%
  summarise(unique_mutations = n_distinct(tf_mutation), .groups = 'drop') %>% 
  arrange(desc(unique_mutations))

# Compute total mutations to find top  TFs
top_TFs <- TF_lumA_classified %>%
  group_by(TF) %>%
  summarise(total_mutations = sum(unique_mutations)) %>%
  arrange(desc(total_mutations))
write.csv(top_TFs, "../results/TFinfluence_TFmutations.csv", row.names = FALSE)

top_TFs <- top_TFs %>% 
  slice_head(n = 20) %>% 
  pull(TF)
  
# Filter the top TFs
mutation_count_top <- TF_lumA_classified %>%
  filter(TF %in% top_TFs) %>% 
  mutate(TF = factor(TF, levels = top_TFs))

png(filename = "../figures/plot_TFinfluence_MAVISp_coverage.png",
    width = 15, height = 6, units = "in", pointsize = 4, res = 1200)
ggplot(mutation_count_top, aes(x = TF, y = unique_mutations, fill = stab_class_dynamic)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "MAVISp Coverage of TFs & Mutations",
    x = "TF",
    y = "Number of Mutations found in LumA cancer patients",
    fill = "MAVISp stability prediction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Pie charts ----
total_TFs <- nrow(TF_lumA_classified)
total_mutations <- sum(TF_lumA_classified$unique_mutations)

not_in_mavisp <- TF_lumA_classified %>%
  filter(stab_class_dynamic == "TF not in MAVISP yet")
not_in_mavisp_count <- nrow(not_in_mavisp)
not_in_mavisp_mutations <- sum(not_in_mavisp$unique_mutations)

percentage_not_in_mavisp <- (not_in_mavisp_count / total_TFs) * 100
percentage_mutations_not_in_mavisp <- (not_in_mavisp_mutations / total_mutations) * 100

cat("Percentage of TFs not in MAVISP yet:", percentage_not_in_mavisp, "%\n")
cat("Percentage of total mutations by these TFs:", percentage_mutations_not_in_mavisp, "%\n")

pie_data <- tibble(
  category = c("TFs not in MAVISP", "Other TFs"),
  count = c(not_in_mavisp_count, total_TFs - not_in_mavisp_count),
  percentage = c(percentage_not_in_mavisp, 100 - percentage_not_in_mavisp)
) %>%
  mutate(label = paste0(category, "\n", count, " (", round(percentage, 1), "%)"))

mutation_data <- tibble(
  category = c("Mutations by TFs not in MAVISP", "Mutations by other TFs"),
  count = c(not_in_mavisp_mutations, total_mutations - not_in_mavisp_mutations),
  percentage = c(percentage_mutations_not_in_mavisp, 100 - percentage_mutations_not_in_mavisp)
) %>%
  mutate(label = paste0(category, "\n", count, " (", round(percentage, 1), "%)"))

tf_pie <- ggplot(pie_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(label = pie_data$label), position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y") +
  labs(
    title = "MAVISp coverage of TFs",
    fill = "Category",
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

mutation_pie <- ggplot(mutation_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(
    title = "Distribution of Mutations",
    fill = "Category",
    x = NULL, y = NULL
  ) +
  scale_fill_manual(values = c("#604cf5", "#f2b25e"),
                    labels = paste0(mutation_data$category, ": ", mutation_data$count, " (", round(mutation_data$percentage, 1), "%)")) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))
png(filename = "../figures/plot_TFinfluence_total_coverage.png",
	width = 8, height = 4, units = "in", pointsize = 4, res = 1200)
grid.arrange(tf_pie, mutation_pie, ncol = 2)

