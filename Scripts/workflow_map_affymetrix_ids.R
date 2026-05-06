library(biomaRt)

testmode<-0

if (testmode==1){
  ensembl <- useEnsembl(biomart = "genes")
  searchDatasets(mart = ensembl, pattern = "mmusculus")
  
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
  
  # example from the biomaRt tutorial
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  affyids=c("202763_at","209310_s_at","207500_at")
  getBM(attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol', 'chromosome_name',
                       'start_position', 'end_position', 'band'),
        filters = 'affy_hg_u133_plus_2', 
        values = affyids, 
        mart = ensembl)
  
  
  ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  affyids=c("17478195", "17491205", "17414747")
  affytypes=c("affy_hc_g110", "affy_hg_focus", "affy_hg_u133a_2", "affy_hg_u133b",
              "affy_hg_u133_plus_2", "affy_hg_u95a", "affy_hg_u95av2", "affy_hg_u95b",
              "affy_hg_u95c", "affy_hg_u95d", "affy_hg_u95e", "affy_hta_2_0", 
              "affy_ht_hg_u133_plus_pm", "affy_huex_1_0_st_v2", "affy_hugenefl",
              "affy_hugene_1_0_st_v1", "affy_hugene_2_0_st_v1", "affy_hugene_2_1_st_v1",
              "affy_primeview", "affy_u133_x3p")
  for (affytype in affytypes){
    print(affytype)
    print(getBM(attributes = c(affytype, 'ensembl_gene_id', 'external_gene_name',
                       'hgnc_id'),
        filters = affytype, 
        values = affyids, 
        mart = ensembl))
  }
}

ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#testattr <- listAttributes(ensembl)
#write.table(testattr, 'mouse_attributes_ensembl.csv', row.names=FALSE)

getBM(attributes = c('affy_mogene_2_1_st_v1', 'ensembl_gene_id', 'external_gene_name'),
      filters = 'affy_mogene_2_1_st_v1', 
      values = affyids, 
      mart = ensembl)

# read table with affy data
inputFile <- './Data/Duparc_et_al/GSE73489_series_matrix_dataonly.txt'
inputData <- read.table(inputFile, header=1)
affyids = inputData$ID_REF
inputData_IDs <- getBM(attributes = c('affy_mogene_2_1_st_v1', 'ensembl_gene_id', 'external_gene_name'),
      filters = 'affy_mogene_2_1_st_v1', 
      values = affyids, 
      mart = ensembl)
write.table(inputData_IDs, 
            "./Data/Duparc_et_al/GSE73489_probe_ids",
            row.names=FALSE)




