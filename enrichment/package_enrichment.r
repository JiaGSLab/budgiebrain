# Load the required packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(AnnotationForge)
library(pkgbuild)
# Create a directory for the annotation package
if (!dir.exists("Annotation_Package")) {
  dir.create("Annotation_Package")
}

#### ---- Annotation ---- ####

# Read the annotation file
eggno <- readr::read_delim(
  file = "/data/input/Files/data/MM_s5dlclua.emapper.annotations.tsv",
  col_names = TRUE,
  delim = "\t",
  comment = "##",
  show_col_types = FALSE
)

# Verify the structure of the data
print(head(eggno))

# Select relevant columns for annotation
emapper <- eggno %>%
  dplyr::select(
    GID = Preferred_name,  # Use Preferred_name as Gene ID
    GO = GOs,              # GO terms
    KO = KEGG_ko,          # KEGG Orthology
    Pathway = KEGG_Pathway # KEGG Pathway
  )

# Prepare gene information with GID as the first column
gene_info <- emapper %>%
  dplyr::select(GID) %>% # Ensure GID is the first column
  dplyr::mutate(Gene_Name = GID) %>%
  dplyr::distinct() # Remove duplicate rows

# Print the first few rows to verify structure
print(head(gene_info))

# Prepare gene to GO mapping with GID as the first column
gene2go <- emapper %>%
  dplyr::select(GID, GO) %>% # Ensure GID is the first column
  tidyr::separate_rows(GO, sep = ",", convert = FALSE) %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::filter(stringr::str_detect(GO, pattern = "GO")) %>%
  dplyr::mutate(EVIDENCE = "IEA") %>%
  dplyr::distinct()

# Print the first few rows to verify structure
print(head(gene2go))

# Ensure that no duplicate column names exist
print(names(gene_info))
print(names(gene2go))

# Create an annotation package using AnnotationForge
AnnotationForge::makeOrgPackage(
  gene_info = gene_info,
  go = gene2go,
  maintainer = 'Hongcheng Shan <shanhongcheng@cibr.ac.cn>',
  author = 'Hongcheng Shan',
  outputDir = ".",
  tax_id = 200540, # Update this with the correct taxonomy ID for Melopsittacus undulatus
  genus = 'Melopsittacus',
  species = 'undulatus',
  goTable = "go",
  version = "1.0"
)

# Build the package
pkgbuild::build('./org.Mundulatus.eg.db/', dest_path = "Annotation_Package")

# Specify the path to the package tar.gz file
package_path <- "/data/work/Annotation_Package/org.Mundulatus.eg.db_1.0.tar.gz"
# Install the package using install.packages
install.packages(
  package_path,
  repos = NULL,
  type = "source"
)

# Verify the package installation
library(org.Mundulatus.eg.db)

# Read the annotation file
eggno <- readr::read_delim(
  file = "/data/work/enrichment_test/MM_s5dlclua.emapper.annotations.tsv",
  col_names = TRUE,
  delim = "\t",
  comment = "##",
  show_col_types = FALSE
)

# Verify the structure of the data
print(head(eggno))

# Select relevant columns for annotation
emapper <- eggno %>%
  dplyr::select(
    GID = Preferred_name,  # Use Preferred_name as Gene ID
    GO = GOs,              # GO terms
    KO = KEGG_ko,          # KEGG Orthology
    Pathway = KEGG_Pathway # KEGG Pathway
  )
# Function to get KEGG pathway ID to name mapping from a local file
get_path2name_local <- function(file_path) {
  # Read the KEGG pathway data from the local file
  tryCatch({
    keggpathid2name.df <- read.delim(file_path, stringsAsFactors = FALSE, header = FALSE)
    
    # Ensure that the data frame has the correct structure
    # Assuming the first column is the pathway ID and the second is the pathway name
    colnames(keggpathid2name.df) <- c("path_id", "path_name")
    
    # Modify path_id to replace 'map' with 'ko'
    keggpathid2name.df <- keggpathid2name.df %>%
      mutate(path_id = str_replace(path_id, pattern = "map", replacement = "ko"))
    
    return(keggpathid2name.df)
  }, error = function(e) {
    message("Failed to read KEGG data from the local file. Please check the file path and format.")
    return(NULL)
  })
}

# Path to the local KEGG pathway file
kegg_file_path <- "/data/work/enrichment_test/kegg_pathway.txt"

# Get the pathway to name mapping from the local file
pathway2name <- get_path2name_local(kegg_file_path)

# Check if the data was loaded successfully
if (!is.null(pathway2name)) {
  # Assume emapper is a data frame already defined with columns 'Pathway' and 'GID'
  # Use dplyr::select to ensure the correct function is used
  pathway2gene <- emapper %>%
    dplyr::select(Pathway, GID) %>%
    separate_rows(Pathway, sep = ",", convert = FALSE) %>%
    filter(str_detect(Pathway, "ko"))
  
  # View the results
  print(head(pathway2name))
  print(head(pathway2gene))
} else {
  message("Pathway to name mapping could not be completed due to file reading failure.")
}
#enrichment example
library(clusterProfiler)
# Load the marker gene list(see DATA_S2)
marker_file <- "/data/work/5month_TeO/Glu/TeO_Glu_level_subcelltype_annotation_markers.csv"
markers <- read_csv(marker_file)

# Check the structure of the markers data
str(markers)

# Filter for the cluster of interest
CABP7_MEIS2_Ex_markers <- markers %>%
  filter(cluster == "CABP7_MEIS2_Ex")

# Inspect the filtered data
head(CABP7_MEIS2_Ex_markers)
# Extract gene symbols or IDs (assuming the column is named 'gene' or 'symbol')
CABP7_MEIS2_Ex_genes <- CABP7_MEIS2_Ex_markers$gene
print(head(CABP7_MEIS2_Ex_genes))
GO_up <- enrichGO(gene = CABP7_MEIS2_Ex_genes,
                  OrgDb = org.Mundulatus.eg.db,
                  keyType = 'GID',
                  ont = 'ALL')
print(dotplot(GO_up,
              label_format = 50,
              showCategory =7))