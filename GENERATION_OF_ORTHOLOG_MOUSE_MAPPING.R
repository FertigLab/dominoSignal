#Modify the ortholog table from Ensembl. The steps to obtain the appropriate inputs from Ensembl are described in "Creation_of_the_human_mouse_ortholog_map_using_Ensembl.txt"

library(dplyr)
library(Seurat)

#Read in the ortholog table from Ensembl
mart_export <- read.delim("mart_export_ENSEMBL_115.txt", stringsAsFactors = FALSE, header=TRUE)

#Reordering columns
mart_export <- mart_export[c(3,4,2,1)]
colnames(mart_export) <- c("hs.ens", "hgnc", "mgi", "mm.ens")

#Saving as RDS
saveRDS(mart_export, "mart_export_ENSEMBL_115.rds")

#Read in RDS
mart_export_ENSEMBL_115 <- readRDS("mart_export_ENSEMBL_115.rds")

##Before you continue, add uniprot to your conversion table?

mouse_uniprot <- read.delim("mart_export_uniprot.txt", stringsAsFactors = FALSE) #This file has colnames Gene stable ID	UniProtKB/Swiss-Prot ID	UniProtKB/TrEMBL ID

colnames(mouse_uniprot) <- c("Gene.stable.ID", "UniProt_SwissProt", "UniProt_TrEMBL")

#Keep the SwissProt ID if it exists, otherwise use the TrEMBL ID. If neither exists, leave as NA. Keep only one UniProt ID per gene, prioritizing SwissProt IDs over TrEMBL IDs when both are available.

mouse_uniprot_1to1 <- mouse_uniprot %>%
  mutate(
    UniProt_SwissProt = na_if(UniProt_SwissProt, ""),
    UniProt_TrEMBL    = na_if(UniProt_TrEMBL, "")
  ) %>%
  mutate(
    uniprot_best = coalesce(UniProt_SwissProt, UniProt_TrEMBL)
  ) %>%
  filter(!is.na(uniprot_best)) %>%
  group_by(Gene.stable.ID) %>%
  summarise(
    uniprot_best = first(uniprot_best),
    .groups = "drop"
  )


ortho_mouse_uniprot <- mart_export_ENSEMBL_115 %>%
  left_join(
    mouse_uniprot_1to1,
    by = c("mm.ens" = "Gene.stable.ID")
  )

human_uniprot_choice <- ortho_mouse_uniprot %>%
  mutate(uniprot_best = na_if(uniprot_best, "")) %>%  # <<< critical
  group_by(hs.ens) %>%
  summarise(
    uniprot_best_chosen = first(na.omit(uniprot_best)),
    .groups = "drop"
  )


final_ortholog_table <- ortho_mouse_uniprot %>%
  left_join(
    human_uniprot_choice,
    by = "hs.ens"
  )

colnames(final_ortholog_table)[colnames(final_ortholog_table)=="uniprot_best_chosen"] <- "uniprot"

saveRDS(final_ortholog_table, "final_mart_ortholog_table_ENSEMBL_115.rds")