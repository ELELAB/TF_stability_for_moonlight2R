###             tf_maf.py
 Usage: program.py <trrust_rawdata.human.tsv> <mutations.csv>

 mutations.csv in this folder is a copy of /data/raw_data/computational_data/TCGA_data_081220/maf_tools/all_muts/cancer_subtypes/BRCA/BRCA.LumA/mutations.csv

 Output: filtered_data which contains:
 # trrust_processed.human.tsv (TF-target-Acti/Rep-DOI)
 # TF_mutations_processed.csv (MAF of TFs with known regulation containing missnense mutations)and their mutational load
 # top_TFS.csv tf of highest mutational load
 # histogram of mutational load
 # barplot of top tfs - mutational load

# Run: 
 tsp bash run.sh




###             MAF_mutations_check.py 
usage: program.py <TF_mutations.csv> <clinvar_metatable.csv> <GENE>

# e.g for TBX3 Run: 
python MAF_mutations_check.py ../filtered_data/TF_mutations_processed.csv ../../../../mavisp/TBX3/cancermuts_25102023/metatable_pancancer_TBX3.csv TBX3 

 #Output: mutations_not_covered_in_cancermuts TBX3


