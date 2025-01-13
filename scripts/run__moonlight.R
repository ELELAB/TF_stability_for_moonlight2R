# This script runs the moonlight pipeline

# SET UP -------------------------------
# Library ----------------
library(TCGAbiolinks)
library(Moonlight2R) 
library(tidyverse)
library(data.table)

options(bitmapType='cairo')
# Sys.setenv('OMP_NUM_THREADS'=3)
# print(Sys.getenv('OMP_NUM_THREADS'))

# Data -------------------
DEA_lumA <- get(load("../rawdata/BRCA_BRCA.LumA_dataDEGs.rda"))
Filt_lumA_raw <- get(load("../rawdata/BRCA_BRCA.LumA_dataFilt_HUGO.rda"))
Maf_lumA <- read.csv("../rawdata/mutations.csv")
TRUSST_data <- read_tsv("../rawdata/trrust_rawdata.human.tsv", col_names = c("TF", "Target", "InteractionType", "DOI")) %>%
  filter(InteractionType %in% c("Activation", "Repression")) #######)

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
print('PRA is done')


# DMA -------
Maf_lumA_filtered <- Maf_lumA %>%
  mutate(VAF = t_alt_count/t_depth) %>% 
  filter(VAF > 0.05, 
         t_alt_count >= 3,
         t_depth > 30,
         n_depth >10,
         Variant_Classification == "Missense_Mutation")

DMA_lumA <- DMA(dataMAF = Maf_lumA_filtered,
                 dataDEGs = DEA_lumA,
                 dataPRA = PRA_lumA,
                 coding_file = "../rawdata/css_coding.vcf.gz",
                 noncoding_file ="../rawdata/css_noncoding.vcf.gz")
save(DMA_lumA, file =  "../results/DMA_lumA.rda")


# Load MAVISp -------
loadMAVISp <- function(mavispDB = NULL,
                       proteins_of_interest = NULL,
                       mode = 'simple',
                       ensemble = 'md'){
    # Look in simple mode index.csv if the protein is in the database
    if (file.exists(str_c(mavispDB,'/dataset_info.csv')) == FALSE){
        stop('MAVISp database file not found at the provided path, or the database_info.csv file is missing.')
    } else if (mode == 'simple'){
        table_location <- str_c(mavispDB,'/simple_mode/dataset_tables/')
    } else if (mode == 'ensemble'){
        if (length(ensemble) > 1){
            stop('Only one ensemble can be specified for ensemble mode.')
        }
        table_location <- str_c(mavispDB,'/ensemble_mode/dataset_tables/')
    } else {
        stop('Mode not specified correctly. Accepts strings "simple" or "ensemble"')
    }

    # Load data for proteins of interest or all proteins
    if (is.null(proteins_of_interest)){
        rawFiles <- list.files(table_location,
                               full.names = TRUE)

    } else {
        proteins_of_interest <- str_c(proteins_of_interest,'-')
        rawFiles <- list.files(table_location,
                               full.names = TRUE) |>
                               as_tibble_col(column_name = 'filepath') |>
                    filter(grepl(paste(proteins_of_interest, 
                                        collapse = '|'),
                                        filepath)) |>
                    pull(filepath)
    }

    # Double check there are files
    if (length(rawFiles) == 0){
        stop('No MAVISp files matching the criteria were found.')
    }

    mavispData <- rawFiles |>
        set_names(str_split_i(basename(rawFiles), '-', 1)) |>
        # Supress non-fatal warnings
        map(function(x) withr::with_options(
                        list(rlib_name_repair_verbosity = 'quiet'),
                        suppressWarnings(
                        classes = 'vroom_parse_issue',
                        read_csv(file = x,
                                progress = FALSE,
                                show_col_types = FALSE))))
    
    # Combine stability results based on user specification
    if (mode == 'ensemble'){
        mavispData <- mavispData |>
            map(function(x) rename(x, 'Stability classification, (Rosetta, FoldX)' = matches(paste0('Stability classification, [A-Za-z0-9, ]*\\(Rosetta, FoldX\\)( \\[',ensemble,'\\])?')),
                                    'Stability classification, (RaSP, FoldX)' = matches(paste0('Stability classification, [A-Za-z0-9, ]*\\(RaSP, FoldX\\)( \\[',ensemble,'\\])?'))) |>
                            mutate(stab_class_ros_source = paste0('ensemble_mode_',ensemble),
                                   stab_class_rasp_source = paste0('ensemble_mode_',ensemble)))
    } else if (mode == 'simple'){
        mavispData <- mavispData |>
            map(function(x) mutate(x, stab_class_data_type = 'simple_mode'))
    }


    return(mavispData)
}



dataMAVISp <- loadMAVISp(mavispDB = "../rawdata/mavisp_csv/11102024_ALL",
                       proteins_of_interest = NULL,
                       mode = 'simple',
                       ensemble = 'md')

# TF -------
TFinfluence <- function(dataPRA,
                        dataDEGs,
                        dataTRRUST,
                        dataMAF,
                        dataMAVISp,
                        stabClassMAVISp = 'rosetta'){ 
    # Control user input -------------
    # dataPRA
    if (all(names(dataPRA) %in% c("TSG", "OCG")) == FALSE) {
        stop("The two list elements in PRA data must be named TSG and OCG")
    }

    # dataDEGs
    if (is.null(dim(dataDEGs)) | nrow(dataDEGs) == 0) {
        stop("The DEG data must be a non-empty table")
    }

    # dataTRRUST
    if (is.null(dim(dataTRRUST)) | nrow(dataTRRUST) == 0) {
        stop("The transcription factor data must be a non-empty table")
    }

    trrust_columns <- c('TF', 'Target', 'InteractionType', 'DOI')

    if (all(trrust_columns %in% names(dataTRRUST)) == FALSE) {
        stop("TRRUST dataframe does not contain the correct columns")
    }

    # dataMAF
    if (is.null(dim(dataMAF)) | nrow(dataMAF) == 0) {
        stop("The mutation data must be a non-empty table")
    }

    maf_columns <- c("Hugo_Symbol",
                "HGVSp_Short")

    if (all(maf_columns %in% names(dataMAF)) == FALSE) {
        stop("MAF file does not contain the correct columns")
    }

    # stabClassMAVISp
    if (!(stabClassMAVISp %in% c('rasp', 'rosetta'))){
        stop("The protocol for stability classification is not specified correctly. Accepts strings 'rasp' or 'rosetta'")
    } else {
        stabClassMAVISp <- paste0('stab_class_', stabClassMAVISp)
    }

    # Load data --------------------------------
    drivers <- PRAtoTibble(dataPRA)

    # Read maf and add ID number to each mutation
    dataMAFFiltered <- dataMAF |>
        #filter(Variant_Classification == "Missense_Mutation") |>  # Only keep rows with Missense mutations
        select(c("Hugo_Symbol",
                "HGVSp_Short")) |>
        mutate(HGVSp_Short = str_extract(HGVSp_Short, pattern = "[A-Z]\\d+[A-Z]")) |>
        rename(mutation = HGVSp_Short) |>
        drop_na()

    # Filter MAVISp data 
    # Keep only the stability classification
    dataMAVISpFiltered <- dataMAVISp |>
                       map(function(x) rename(x, 'stab_class_rosetta' = matches('Stability classification, [A-Za-z0-9, ]*\\(Rosetta, FoldX\\)'),
                                                 'stab_class_rasp' = matches('Stability classification, [A-Za-z0-9, ]*\\(RaSP, FoldX\\)'))) |>
                       keep(function(x) stabClassMAVISp %in% colnames(x)) |>
                       map(function(x) rename(x, 
                                              'mutation' = 1,
                                              'stab_class' = stabClassMAVISp) |>
                                       select(mutation, stab_class)) |>
                       rbindlist(idcol = 'protein') |>
                       as_tibble()
      
    # Add a flag to dataMAVISpFiltered to indicate presence
    dataMAVISpFiltered <- dataMAVISpFiltered |> 
      mutate(in_MAVISp = TRUE)  # Add a flag indicating the TF is present in MAVISp
    
    

    # Analysis -------------------
    # Convert rownames to column for DEGs
    dataDEGs <- dataDEGs |>
        rownames_to_column(var = "GENE") |> # Target
        select(GENE, logFC, adj.P.Val) 

    # Join TFs drivers with expression
    drivers_expr <- drivers |>
        left_join(dataDEGs, 
                  by = join_by(Hugo_Symbol == GENE))

    # Join drivers with TF
    drivers_expr_tf <- drivers_expr |>
        left_join(dataTRRUST,
                  by = join_by(Hugo_Symbol == Target)) |> # 
        rename('logFC_target' = logFC,
               'Target' = Hugo_Symbol) #

    # Map mutation file to TF
    drivers_TF_mut <- drivers_expr_tf |>
        left_join(dataMAFFiltered,
                  by = join_by(TF == Hugo_Symbol),
                  relationship = "many-to-many") |>
        rename('tf_mutation' = mutation)

    # Match TF-mut with mavisp to see the stability effect
    drivers_mut_mavisp <- drivers_TF_mut |>
      left_join(dataMAVISpFiltered, 
                by = join_by(TF == protein, tf_mutation == mutation)) |> 
      # Add helper column to check if TF is in MAVISp
      mutate(in_MAVISp = ifelse(TF %in% unique(dataMAVISpFiltered$protein), TRUE, FALSE)) |> 
      # Determine stab_class based on helper column and mutation presence
      mutate(stab_class = case_when(
        !in_MAVISp ~ 'not in MAVISp yet',                  # TF not in MAVISp
        in_MAVISp & is.na(stab_class) ~ 'mutation in IDR',  # TF in MAVISp but mutation missing
        TRUE ~ stab_class                                  # Keep existing stab_class if no issue
      )) |> 
      # Drop the helper column
      select(-in_MAVISp)
    
    # Check biological implications
    # Activation -> destabilising mutation -> decrease
    # Repression -> destabilising mutation -> increase
    filtered_effect <- drivers_mut_mavisp |>
      filter((InteractionType == 'Activation' & stab_class == 'Destabilizing' & logFC_target < 0) |
               (InteractionType == 'Repression' & stab_class == 'Destabilizing' & logFC_target > 0) |
               (stab_class == 'Uncertain') |
               (stab_class == 'Neutral'))
    
    return(drivers_mut_mavisp)
}



TF_lumA <- TFinfluence(dataPRA = PRA_lumA,
            dataTRRUST = TRUSST_data,
            dataMAF = Maf_lumA,
            dataDEGs = DEA_lumA,
            dataMAVISp = dataMAVISp
            )
save(TF_lumA, file =  "../results/TF_lumA.rda")

# Plots -----
mutation_count <- TF_lumA %>%
  group_by(TF, stab_class) %>%         
  summarise(unique_mutations = n_distinct(tf_mutation), 
            .groups = 'drop') 

top_TFs <- mutation_count %>%
  group_by(TF) %>%
  summarise(total_mutations = sum(unique_mutations)) %>%
  arrange(desc(total_mutations))  
write.csv(top_TFs, file = "../results/TFinfluence_TFmutations.csv", row.names = FALSE)

top_15 <- top_TFs %>%
  slice_head(n = 15) %>%
  pull(TF) 

mutation_count_top <- mutation_count %>%
  filter(TF %in% top_15)

png(filename = "../figures/plot_TFinfluence_top15.png",
    width = 15, height = 6, units = "in", pointsize = 4, res = 1200)
ggplot(mutation_count_top, aes(x = TF, y = unique_mutations, fill = stab_class)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Top 15 TFs: Mutations affecting the stability of TFs",
    x = "TF",
    y = "Number of Mutations",
    fill = "MAVISp stability class"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


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






