read_if_char <- function(x, ...) { #This is a necessary function that isn't in my install of dominoSignal
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    ext <- tools::file_ext(x)
    
    if (ext %in% c("tsv", "txt")) {
      return(readr::read_tsv(x, ...))
    } else if (ext == "csv") {
      return(readr::read_csv(x, ...))
    } else if (ext == "rds") {
      return(readRDS(x))
    } else {
      stop("Unsupported file extension: ", ext)
    }
  }
  
  return(x)
}


create_rl_map_matricom <- function(
    genes,
    proteins,
    interactions,
    kegg_lr, #This is new, and it was added in because OmniPath seems to lack many COL-ITG interactions present in KEGG
    ortholog_table,
    complexes = NULL,
    database_name = "MatriComDB",
    save_path = NULL
) {
  
  ## ------------------------------------------------------------------
  ## Stage 0 — Input checks & loading. Changed to BaseR to avoid package issues with dreamerr
  ## ------------------------------------------------------------------
  
  stopifnot(
    is.character(genes)       || is.data.frame(genes),
    is.character(proteins)    || is.data.frame(proteins),
    is.character(interactions)|| is.data.frame(interactions),
    is.character(kegg_lr)     || is.data.frame(kegg_lr),
    is.character(ortholog_table) || is.data.frame(ortholog_table),
    is.character(database_name) && length(database_name) == 1
  )
  
  if (!is.null(save_path)) {
    stopifnot(is.character(save_path), length(save_path) == 1)
  }
  
  genes          <- read_if_char(genes) #This section is unchanged
  proteins       <- read_if_char(proteins)
  interactions   <- read_if_char(interactions)
  kegg_lr        <- read_if_char(kegg_lr)
  ortholog_table <- read_if_char(ortholog_table)
  
  ## Explicitly disable complexes. These will have to be added back in later
  complexes <- NULL
  
  ## ------------------------------------------------------------------
  ## Stage 1 — Build base RL map (HUMAN)
  ## ------------------------------------------------------------------
  gene_features <- c("gene", "category", "family", "label", "uniprot")
  
  annotation_tbl <- genes %>%
    dplyr::select(any_of(gene_features)) %>% #left joining to avoid issues with factor levels and names being converted into numbers
    dplyr::left_join(
      proteins %>% dplyr::select(any_of(gene_features)),
      by = "gene"
    ) %>%
    tidyr::unite(
      "annotation_strategy",
      label, category, family,
      sep = ";",
      na.rm = TRUE,
      remove = TRUE
    ) %>%
    dplyr::distinct(gene, .keep_all = TRUE) #keep only distinct entries
  
  rl_map <- NULL #constructing the rl map skeleton
  
  safe_first <- function(x) { #I need to fill blank fields with NA prior to filling these in from the conversion table
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(NA_character_)
    }
    x[1]
  }
  
  for (i in seq_len(nrow(interactions))) { #Follows logic from the create_rl_map_cellphonedb function
    
    inter <- interactions[i, ]
    
    partner_a <- inter[["V1"]]
    partner_b <- inter[["V2"]]
    
    if (is.na(partner_a) || is.na(partner_b)) next
    
    anno_a <- annotation_tbl[annotation_tbl$gene == partner_a, ]
    anno_b <- annotation_tbl[annotation_tbl$gene == partner_b, ]
    
    if (nrow(anno_a) == 0 || nrow(anno_b) == 0) next #Skip if annotation missing
    
    row_i <- data.frame(
      name_A = partner_a,
      gene_A = toupper(partner_a),
      uniprot_A = NA_character_, #These have to come from the ortholog table, since the original annotation table has no uniprot information. We will fill these in later.
      type_A = NA_character_,
      name_B = partner_b,
      gene_B = toupper(partner_b),
      uniprot_B = NA_character_,
      type_B = NA_character_,
      annotation_strategy = safe_first(anno_a$annotation_strategy),
      source = as.character(inter[["source"]])[1],
      database_name = database_name,
      stringsAsFactors = FALSE
    )
    
    rl_map <- rbind(row_i, rl_map)
  }
  
  rl_map <- rl_map %>%
    dplyr::distinct(gene_A, gene_B, .keep_all = TRUE) # <<< Critical, and we need this upfront to prevent interaction space from blowing up due to gene duplication
  
  ## ------------------------------------------------------------------
  ## Stage 2 — OmniPath annotation (HUMAN)
  ## ------------------------------------------------------------------
  lr_network <- OmnipathR::intercell_network( #Takes high-confidence human lr interactions from OmniPath
    intercell_consensus_percentile = 33,
    ligand_receptor = TRUE
  )
  
  op_lr <- lr_network %>%
    dplyr::select(source_genesymbol, target_genesymbol) %>%
    tidyr::separate_longer_delim(source_genesymbol, "_") %>% #Split complexes by "_"
    tidyr::separate_longer_delim(target_genesymbol, "_") %>%
    dplyr::transmute(
      gene_A = toupper(source_genesymbol), #Harmonize gene symbols across datasets (safe) for easy merge
      gene_B = toupper(target_genesymbol),
      type_A = "L",
      type_B = "R"
    ) %>%
    dplyr::bind_rows( #The order of the genes in my original map is arbitrary, so I need to add in the flipped interactions to make sure I can merge with the OmniPath data regardless of the order of genes in the OmniPath interactions. This also allows me to annotate interactions that are present in my original map but flipped in OmniPath.
      dplyr::transmute(
        .,
        gene_A = gene_B,
        gene_B = gene_A,
        type_A = "R",
        type_B = "L"
      )
    ) %>%
    dplyr::distinct()
  
  rl_map <- rl_map %>%
    dplyr::left_join(op_lr, by = c("gene_A", "gene_B"), suffix = c(".old", "")) %>%
    dplyr::mutate(
      type_A = dplyr::coalesce(type_A, type_A.old), #Coalesce to avoid having multiple columns with type_A and type_B. We want to overwrite the original type_A and type_B columns with the new ones from OmniPath, but if there is no annotation from OmniPath, we want to keep the original NA values rather than overwriting them with NA from the join.
      type_B = dplyr::coalesce(type_B, type_B.old),
      source = ifelse(
        !is.na(type_A) & !grepl("OmniPath", source),
        paste(source, "OmniPath_annotated", sep = ";"),
        source
      )
    ) %>%
    dplyr::select(-ends_with(".old")) #Drop unnecessary columns
  
  ## ------------------------------------------------------------------
  ## Stage 3 — KEGG annotation (HUMAN). Needed bc OmniPath seems to lack many COL-ITG interactions present in KEGG, and these are important for my downstream analyses. 
  ## ------------------------------------------------------------------
  kegg_lr <- kegg_lr %>%
    dplyr::mutate(
      gene_A = toupper(Gene1), #Gene symbol harmonization for safe merging
      gene_B = toupper(Gene2),
      type_A = ifelse(dir == "LR", "L", "R"),
      type_B = ifelse(dir == "LR", "R", "L") #The other orientation should have type_A as R and type_B as L
    ) %>%
    dplyr::select(gene_A, gene_B, type_A, type_B) %>%
    dplyr::distinct()
  
  rl_map <- rl_map %>%
    dplyr::left_join(kegg_lr, by = c("gene_A", "gene_B"), suffix = c("", ".kegg")) %>%
    dplyr::mutate(
      type_A = dplyr::coalesce(type_A, type_A.kegg),
      type_B = dplyr::coalesce(type_B, type_B.kegg),
      source = ifelse(
        !is.na(type_A.kegg) & !grepl("KEGG_hsa04512", source),
        paste(source, "KEGG_hsa04512", sep = ";"), #Add KEGG annotation to source column if the interaction is annotated by KEGG and the source column doesn't already indicate that it's been annotated by KEGG. We want to keep existing annotations in the source column, so we concatenate the new annotation with the existing source rather than overwriting it.
        source
      )
    ) %>%
    dplyr::select(-ends_with(".kegg")) #Remove unneeded columns from the join
  
  ## ------------------------------------------------------------------
  ## Stage 4 — Resolve unannotated interactions (HUMAN)
  ## ------------------------------------------------------------------
  annotated <- rl_map %>% dplyr::filter(!is.na(type_A) & !is.na(type_B))
  unannotated <- rl_map %>% dplyr::filter(is.na(type_A) & is.na(type_B))
  
  if (nrow(unannotated) > 0) {
    unannotated <- unannotated %>%
      dplyr::slice(rep(1:n(), each = 2)) %>% #Duplicate rows that don't already have LR annotations, and make such that either partner can be L or R
      dplyr::group_by(gene_A, gene_B) %>%
      dplyr::mutate(
        type_A = ifelse(dplyr::row_number() %% 2 == 1, "R", "L"),
        type_B = ifelse(dplyr::row_number() %% 2 == 1, "L", "R")
      ) %>%
      dplyr::ungroup()
  }
  
  rl_map_human <- dplyr::bind_rows(annotated, unannotated) #Stitch everything together
  
  ## ------------------------------------------------------------------
  ## Stage 5 — Human → Mouse ortholog conversion (PARALOG-AWARE)
  ## ------------------------------------------------------------------
  
  
  orth_map <- ortholog_table %>% #Gather columns we need from our ortholog map
    dplyr::select(
      hgnc,
      mgi,
      uniprot
    ) %>%
    dplyr::distinct()
  
  
  ## ------------------------------------------------------------------
  ## Stage 5a — Populate UniProt IDs from ortholog_table
  ## ------------------------------------------------------------------
  
  # Join UniProt for gene_A. This follows from GENERATION_OF_ORTHOLOG_MOUSE_MAPPING.R, where we prioritized SwissProt IDs over TrEMBL IDs, and thus have only one UniProt ID per gene in the ortholog table. We can therefore directly join and coalesce without worrying about multiple UniProt IDs per gene.
  rl_map_human <- rl_map_human %>%
    dplyr::left_join(
      orth_map %>%
        dplyr::select(hgnc, uniprot) %>%
        dplyr::rename(uniprot_join = uniprot),
      by = c("gene_A" = "hgnc")
    ) %>%
    dplyr::mutate(
      uniprot_A = dplyr::coalesce(uniprot_join, uniprot_A) #coalesce to avoid having multiple columns with uniprot. We want to overwrite
    ) %>%
    dplyr::select(-uniprot_join)#Remove extraneous column
  
  # Join UniProt for gene_B
  rl_map_human <- rl_map_human %>%
    dplyr::left_join(
      orth_map %>%
        dplyr::select(hgnc, uniprot) %>%
        dplyr::rename(uniprot_join = uniprot),
      by = c("gene_B" = "hgnc")
    ) %>%
    dplyr::mutate(
      uniprot_B = dplyr::coalesce(uniprot_join, uniprot_B)
    ) %>%
    dplyr::select(-uniprot_join)
  
  
  ## Expand paralog combinations. In mice, there are multiple expansion events for homologous genes (interaction-ready mouse map)
  rl_map_mouse <- rl_map_human %>%
    # map gene_A
    dplyr::left_join(
      orth_map %>% dplyr::select(hgnc, mgi),
      by = c("gene_A" = "hgnc")
    ) %>%
    dplyr::rename(gene_A_mouse = mgi) %>%
    # map gene_B
    dplyr::left_join(
      orth_map %>% dplyr::select(hgnc, mgi),
      by = c("gene_B" = "hgnc")
    ) %>%
    dplyr::rename(gene_B_mouse = mgi) %>%
    # drop rows that cannot be mapped to mouse
    dplyr::filter(
      !is.na(gene_A_mouse),
      !is.na(gene_B_mouse)
    )
  
  rl_map_mouse <- rl_map_mouse %>% #We need to rename our columns to match names in conventional domino rl map
    dplyr::transmute(
      name_A = gene_A_mouse,
      gene_A = gene_A_mouse,
      uniprot_A = uniprot_A,
      name_B = gene_B_mouse,
      gene_B = gene_B_mouse,
      uniprot_B = uniprot_B,
      annotation_strategy = annotation_strategy,
      source = source,
      database_name = database_name,
      type_A = type_A,
      type_B = type_B
    )
  
  stopifnot(ncol(rl_map_mouse) == 11) #Just making sure the mouse map has the same number of columns as the human map
  
  rl_map_human <- rl_map_human %>% #int_pair needs to be created separately for both human and mouse datasets from name_A and name_B
    dplyr::mutate(int_pair = paste(name_A, name_B, sep = " & ")) %>%
    dplyr::relocate(int_pair)
  
  rl_map_mouse <- rl_map_mouse %>%
    dplyr::mutate(int_pair = paste(name_A, name_B, sep = " & ")) %>%
    dplyr::relocate(int_pair)
  
  rl_map_human <- rl_map_human %>% #Distinct is used to prevent artificial inflation of int due to Cartesian multiplication of duplications
    dplyr::distinct(gene_A, gene_B, type_A, type_B, .keep_all = TRUE)
  
  rl_map_human <- rl_map_human[, c(
       "int_pair", "name_A", "uniprot_A", "gene_A", "type_A", "name_B", "uniprot_B",
       "gene_B", "type_B", "annotation_strategy", "source", "database_name"
     )]
  
  rl_map_mouse <- rl_map_mouse %>%
    dplyr::distinct(gene_A, gene_B, type_A, type_B, .keep_all = TRUE)
  
  rl_map_mouse <- rl_map_mouse[, c(
    "int_pair", "name_A", "uniprot_A","gene_A", "type_A", "name_B", "uniprot_B",
    "gene_B", "type_B", "annotation_strategy", "source", "database_name"
  )]
  
  ## ------------------------------------------------------------------
  ## Step 6 — Save & return human and mouse maps (RDS), since we already generated the human as intermediate file for mouse rl generation.
  ## ------------------------------------------------------------------
  
  if (!is.null(save_path)) {
    
    saveRDS(
      rl_map_human,
      file = file.path(save_path, "rl_map_human.rds")
    )
    
    saveRDS(
      rl_map_mouse,
      file = file.path(save_path, "rl_map_mouse.rds")
    )
  }
  
  return(
    list(
      rl_map_human = rl_map_human,
      rl_map_mouse = rl_map_mouse
    )
  )
}

#The below was just a sample of running the code

save_path=getwd()

library(dominoSignal)
library(dplyr)
library(tidyr)
library(Seurat)
library(OmnipathR)

#Read in my KEGG ligand-receptor file, which is a csv downloaded from the MatriCom Github (hsa04512_LR-RL.csv)

kegg_lr <- read.csv("KEGG-hsa04512_LR-RL.csv")

#Read in the mouse ortholog mapping, which is a tsv I generated from Ensembl (final_mart_ortholog_table_ENSEMBL_115.txt)

ortholog_table <- readr::read_tsv("final_mart_ortholog_table_ENSEMBL_115.txt")

#Read in the reference files, which are rda objects. The function does not yet have complexes, and these need to be added back in.

genes<- load(CCgenes2.rda)
proteins <- load(mlist.rda)
interactions <- load(intlist.rda)


create_rl_map_matricom(
  genes = CCgenes2,
  proteins = mlist,
  interactions = intlist,
  kegg_lr = kegg_lr,
  ortholog_table = ortholog_table,
  database_name = "MatriComDB",
  save_path = save_path
)

#I also wrote the output to tsv for easier inspection, but the RDS files are the main output for downstream use

write.table(rl_map_human, file = "rl_map_human.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(rl_map_mouse, file = "rl_map_mouse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(intlist, file = "intlist.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(op_lr, file = "op_lr.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#The results of some external checking is below.

###Unique genes in mouse RL map: 5,189
###Genes shared between Seurat (Batch3_Full) & RL map 4,853
###% of Seurat genes present in RL map 19.06%
###% of RL map genes present in Seurat object 93.52%
###The above fits into the biologically plausible range (doi: 10.1093/bioinformatics/btab036; https://doi.org/10.1093/nar/gkaf1108)

#Now, I just want to make sure some COL-ITG interactions have correct annotations in mouse map and some of the less well-known interactions are annotated ambiguously (i.e. with both orientations). 

#From KEGG hsa04512, we have the following interactions:
#COL9A3	ITGA11	COL9A3_ITGA11	LR 1.)
#LAMB2	ITGA6	LAMB2_ITGA6	LR 2.)
#FN1	ITGA2B	FN1_ITGA2B	LR 3.)

#From OmniPath, we have the following interactions:

#DCN TLR4 LR 4.)
#VTN PLAUR LR 5.)
#THBS1 LRP1 LR 6.)

#Present in neither, we have the following interactions:

#NID2 – LAMA2 7.)
#EFEMP2 – ELN 8.)
#FN1 – DCN 9.)


#Checking my human map to see if items 1-9 are annotated correctly. 

check_vector <- c("COL9A3 & ITGA11", "LAMB2 & ITGA6", "FN1 & ITGA2B", "DCN & TLR4", "VTN & PLAUR", "THBS1 & LRP1", "NID2 & LAMA2", "EFEMP2 & ELN", "FN1 & DCN")
rl_map_human %>%
  filter(int_pair %in% check_vector) %>%
  select(int_pair, type_A, type_B, source)

#Ok, so "DCN & TLR4", "VTN & PLAUR" not present in human map. Checking to see if these are in intlist

intlist %>%
  filter(
    (V1 == "DCN" & V2 == "TLR4") |
      (V1 == "TLR4" & V2 == "DCN")
  )

#Returns 0 rows. So, these are interactions that are present in OmniPath but not in my original MatriCom interaction list, which is why they weren't annotated in the human map. This is expected.
#I think I should be okay here. One other thing to check is that the str of the mouse map resembles previous domino RL maps. Here, I am using the cpdb map with expanded complexes

str(rl_map_mouse) #All chr variables with 12 columns
str(cpdb_rl_map_expanded) #All chr variables with 12 columns. cpdb_rl_map_expanded is the expanded version of the CellPhoneDB RL map that I generated in CREATE_RL_MAP_CELLPHONEDB.R. Note that this is a tibble and not a data frame, but the structure is the same otherwise.

#Check with Kavita. Do we need a tibble to run domino, or can the rl map be a data frame? If we need a tibble, then I can convert the mouse map to a tibble and save it as an RDS. If we can run domino with a data frame, then I can just save the mouse map as an RDS without converting to a tibble.

#Next step: send the source files and code to Kavita, opening a new sub-branch for dominoSignal named 'Varun'. Push to GitHub.