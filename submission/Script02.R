### ************************************
### PREFACE ----
### ************************************

# this script houses the entirety of code for analysis performed in R^
# ^after sequence processing performed with QIIME 2 (see Script01.R)
# ** note for KDP: PREFACE version 0.5 ** #

# configure the script
# NOTE: my machine runs MacOS; paths have not been tested on other platforms
# NOTE: this information and more is stored in section I below

setwd("~/Desktop/") # *gasp*
name_scrpt <- "Script02" # script filename
name_projd <- "TransFaunation" # name of project directory
path_projd <- paste(sep = "", getwd(), "/", name_projd) # path to project dir
path_vault <- paste(sep = "", path_projd, "/vault") # path to storage vault
PREFACE <- c("name_scrpt", "name_projd", "path_projd", "path_vault")

# paths to universal inputs stored in the project directory (path_projd)
ifp_smp_dat <- paste(sep = "", path_projd, "/sampledata_TF.txt")

# paths to section specific inputs stored in the central vault (path_vault)
A.ifv_fts_tab <- paste(sep = "", path_vault, "/q2_fts_tab.tsv")
A.ifv_rep_seq <- paste(sep = "", path_vault, "/q2_rep_seq.fasta")
A.ifv_int_ggs <- paste(sep = "", path_vault, "/q2_int_ggs.tsv")
A.ifv_int_slv <- paste(sep = "", path_vault, "/q2_int_slv.tsv")

# paths for outputs stored in the central vault (path_vault)
ofv_COSMOS_wksp <- paste(sep = "", path_vault, "/WS_COSMOS.RData")

I.ofv_info <- paste(sep = "", path_vault, "/info_", name_scrpt, "_SectionI.txt")
I.ofv_wksp <- paste(sep = "", path_vault, "/WS_", name_scrpt, "_SectionI.RData")

A.ofv_EBTKS_raw <- paste(sep = "", path_vault, "/table_EBTKS_abs_raw.txt")
A.ofv_EBTKS_pro <- paste(sep = "", path_vault, "/table_EBTKS_abs_processed.txt")
A.ofv_info <- paste(sep = "", path_vault, "/info_", name_scrpt,"_SectionA.txt")
A.ofv_prov <- paste(sep = "", path_vault, "/prov_", name_scrpt,"_SectionA.txt")
A.ofv_wksp <- paste(sep = "", path_vault, "/WS_", name_scrpt,"_SectionA.RData")

S.ofv_info <- paste(sep = "", path_vault, "/info_", name_scrpt,"_SectionS.txt")
S.ofv_wksp <- paste(sep = "", path_vault, "/WS_", name_scrpt,"_SectionS.RData")

T.ofv_EBTKS_tree <- paste(sep = "", path_vault, "/tree_EBTKS_processed.nwk")
T.ofv_prov <- paste(sep = "", path_vault, "/prov_", name_scrpt,"_SectionT.txt")
T.ofv_wksp <- paste(sep = "", path_vault, "/WS_", name_scrpt,"_SectionT.RData")

# save PREFACE workspace
PREFACE.lst <- c(ls(pattern = "ifp"), ls(pattern = "ifv"), ls(pattern = "ofv"), 
                 PREFACE, "PREFACE", "PREFACE.lst")
save(list = PREFACE.lst,
     file = paste(sep = "", path_vault, "/WS_", name_scrpt, "_PREFACE.RData"))

### ************************************
### UNIVERSAL OBJECTS ~ COSMOS ----
### ************************************

# this section contains items that are universal across datasets

# read in the sample data (commonly called metadata)
smp_dat <- read.table(ifp_smp_dat, sep = "\t", header = T, as.is = T, 
                      stringsAsFactors = F)

# hex codes for a gradient of black/grey/white (left to right = dark to light)
greydient <- c("#000000", "#252525", "#525252",
               "#969696", "#bbbbbb", "#d9d9d9", "#e0e0e0", 
               "#ffffff")

# Life Aquatic with Steve Zissou custom color palettes
ziss_blue <- c("#218ec4", "#2ca1db", "#67bbe5", "#92ceec", "#e9f5fb", "#ffffff")
ziss_navy <- c("#0f1b3e", "#192e67", "#2d52b9", "#466cd2", "#eaeffa", "#ffffff")
ziss_reds <- c("#b31300", "#f21a00", "#ff604d", "#ffbbb3", "#ffe8e5", "#ffffff")
ziss_ylow <- c("#fec72A", "#fede80", "#ffebb3", "#fff8e6", "#ffffff")
ziss_teal <- c("#285854", "#489e97", "#84c7c2", "#cae7e5", "#edf7f6", "#ffffff") 
ziss_grey <- c("#252525", "#3d3d3d", "#c4c4c4", "#d9d9d9", "#e9e9e9", "#ffffff")

# font family
fnt <- "Courier"

# multi-panel label ~ font parameters
pan_fnt <- list(size = 13, family = fnt, face = "bold")

# plot widths for single and double column
wid_sgl <- 78
wid_dbl <- 165

# vector with column names commonly removed or retained
com_col <- c("FeatureID", "RepSeq", 
             "int.ggs.tax", "int.ggs.Kingdom", "int.ggs.Phylum", 
             "int.ggs.Class", "int.ggs.Order", "int.ggs.Family", 
             "int.ggs.Genus", "int.ggs.Species", "int.ggs.cnf", 
             "int.ggs.lws.txn", "int.ggs.lws.lvl", 
             "int.slv.tax", "int.slv.Kingdom", "int.slv.Phylum", 
             "int.slv.Class", "int.slv.Order", "int.slv.Family", 
             "int.slv.Genus", "int.slv.Species", "int.slv.cnf", 
             "int.slv.lws.txn", "int.slv.lws.lvl")

# ^^cutpoints and symbols to be used for determining P & BH_P value significance
# NOTE: used as argumants for symnum() base function
cutpts <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
symbls <- c("****", "***", "**", "*", " ")

cutpts_P <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
symbls_P <- c("****", "***", "**", "*", " ")
cutpts_BHP <- c(0, 0.0001, 0.001, 0.01, 0.1, 1)
symbls_BHP <- c("****", "***", "**", "*", "ns")

# create a vector naming all of the above (and include ifp_smp_dat from PREFACE)
# NOTE: this is useful when saving workspaces
cosmos <- c("ziss_blue", "ziss_navy", "ziss_reds", "ziss_ylow", "ziss_teal", 
            "ziss_grey",
            "ifp_smp_dat", "smp_dat", "greydient", "fnt", "pan_fnt", "com_col",
            "cutpts", "symbls",
            "cutpts_P", "symbls_P", "cutpts_BHP", "symbls_BHP",
            "wid_sgl", "wid_dbl")

### ************************************
### UNIVERSAL OBJECTS ~ dataset/project specific ----
### ************************************

# this section contains items that are not universal across datasets ...
# ... but are universal across sections within a specific dataset or project

# abbreviations used below:
# shp = shape; drk = dark; med = medium; lyt = light; 
# cec = cecum; prx = proximal colon; dtl = distal colon; fec = feces; 
# ino = inoculum; plm = plasma;
# nmc = no microbial consortium; cmc = control microbial consortium;
# rmc = rice bran modified microbial consortium

# ggplot2 codes for shapes by sample type
shp_cec <- 24 # open triangle (outline dark; fill medium)
shp_prx <- 21 # open circle (outline medium; fill light)
shp_dtl <- 21 # open circle (outline dark; fill medium)
shp_fec <- 25 # open upside triangle (outline medium; fill mixed)
shp_ino <- 32 # blank (will be plotted as text)
shp_plm <- 23 # open diamond (outline medium; fill mixed)

# ggplot2 codes for shapes by murine sex
shp_mur_fem <- 04 # X
shp_mur_mal <- 22 # open square (filled with white)

# hex codes to color no inoculum (no microbial consortium) group
hex_lyt_nmc <- "#ffffff" # white
hex_drk_nmc <- "#252525" # black

# hex codes to color human female inoculum groups (darkest to lightest)
# hex_F01_cmc <- c("#0f1b3e", "#192e67", "#2d52b9", 
#                  "#466cd2", "#eaeffa") # Zissou navy
# hex_F01_rmc <- c("#285854", "#489e97", "#84c7c2", 
#                  "#cae7e5", "#edf7f6") # Zissou teal

hex_F01_cmc <- c("#0f1b3e", "#192e67", "#2d52b9",
                 "#466cd2", "#eaeffa") # Zissou navy
hex_F01_rmc <- c("#b31300", "#f21a00", "#ff604d", "#ffbbb3", 
                 "#ffe8e5") # Zissou red

# hex codes to color human male inoculum groups (darkest to lightest)
hex_M00_cmc <- c("#252525", "#3d3d3d", "#c4c4c4", 
                 "#d9d9d9", "#e9e9e9") # black
hex_M00_rmc <- c("#b28401", "#fec72A", "#fede80", 
                 "#ffebb3", "#fff8e6") # Zissou yellow
# linetypes by consortium group
lne_typ_cmc <- "solid"
lne_typ_rmc <- "dashed"

# create a vector naming all of the above (useful when saving workspaces)
TransFaun <- c("shp_cec", "shp_prx", "shp_dtl", "shp_fec", "shp_ino", "shp_plm",
               "shp_mur_fem", "shp_mur_mal", 
               "hex_lyt_nmc", "hex_drk_nmc", 
               "hex_F01_cmc", "hex_F01_rmc",
               "hex_M00_cmc", "hex_M00_rmc",
               "lne_typ_cmc", "lne_typ_rmc")

# save universal objects workspace
COSMOS <- c(cosmos, "COSMOS", TransFaun)
save(list = COSMOS, file = ofv_COSMOS_wksp)

### ************************************
### ^^^I - MACHINE/PACKAGE/VERSION Info ----
### ************************************

# ^^^ = check "# capture R package-related information:" to ensure it is correct

# NOTE: section I requires objects from the PREFACE to be in the environment

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)

# R packages accessed via namespace:
# benchmarkme
# dplyr
# DECIPHER
# Biostrings
# phangorn
# ape
# ggbiplot
# GUniFrac
# vegan
# zCompositions

# capture R package-related information:
I.Rpac_ctg <- "R package version"
I.Rpackge_a <- data.frame(info = "benchmarkme", section = "I", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_b <- data.frame(info = "dplyr", section = "I; A; S; B; C; D", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_c <- data.frame(info = "ape", section = "T; B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_d <- data.frame(info = "Biostrings", section = "T", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_e <- data.frame(info = "DECIPHER", section = "T", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_f <- data.frame(info = "phangorn", section = "T", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_g <- data.frame(info = "ggplot2", section = "B; C", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_h <- data.frame(info = "ggpubr", section = "B; C", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_i <- data.frame(info = "ggbiplot", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_j <- data.frame(info = "GUniFrac", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_k <- data.frame(info = "vegan", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_l <- data.frame(info = "zCompositions", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_0 <- rbind(I.Rpackge_a, I.Rpackge_b, I.Rpackge_c, I.Rpackge_d, 
                     I.Rpackge_e, I.Rpackge_f, I.Rpackge_g, I.Rpackge_h, 
                     I.Rpackge_i, I.Rpackge_j, I.Rpackge_k, I.Rpackge_l)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: captureVersions() version 0.1 ** #
# function takes an input data.frame with:
# column 'info' listing the package name and returns a data.frame with a ...
# ... new column 'value' listing the version information for the package
# NOTE: input column names are rigid; column 'info' must list package name
captureVersions <- function(data = data.frame) {
  # internal checks to ensure correct input classes
  if (!inherits(data, "data.frame")) {
    stop("input data must be class 'data.frame'")
  }
  # store original input df and format the new df
  new_dat <- data
  new_dat$value <- ""
  for(i in 1:length(new_dat$info)) {
    # create package name variable to be tested to ensure correct inputs/outputs
    pack <- unlist(
      packageDescription(
        new_dat$info[i], fields = c("Package", "Version"))[1])
    # test if package name in input data matches 'pack' variable
    # if TRUE, add package version value to the new data.frame
    if (pack == new_dat$info[i]) {
      new_dat$value[i] <- unlist(
        packageDescription(
          new_dat$info[i], fields = c("Package", "Version"))[2])
    }
    # if FALSE, print error message
    else {
      if (!pack == new_dat$info[i]) {
        stop("'pack' variable returned by packageDescription();
       does not match row value in column 'package' of input data;
       incorrect package version valuermation likely returned")
      }
    }
  }
  return(new_dat)
}
# 
# # example usage:
# new_dataframe <- captureVersions(data = dataframe)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# captureVersions() function
I.Rpackge_1 <- captureVersions(I.Rpackge_0)
I.Rpackge <- dplyr::select(I.Rpackge_1, category, info, value, script, section)

# capture Project-related information:
I.project <- data.frame(category = "project", 
                        info = "name of project directory",
                        value = name_projd,
                        script = name_scrpt,
                        section = "all", stringsAsFactors = F)

# capture PATH-related information:
I.path_ctg <- "filepath"
I.path_sec <- "all"
I.pathsto_a <- data.frame(info = "working directory",
                          value = paste(getwd(), "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_b <- data.frame(info = "path to project directory",
                          value = paste(path_projd, "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_c <- data.frame(info = "path to project's central '/vault/'",
                          value = paste(path_vault, "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_0 <- rbind(I.pathsto_a, I.pathsto_b, I.pathsto_c)
I.pathsto <- dplyr::select(I.pathsto_0, category, info, value, script, section)


# capture Machine-related information:
I.mach_ctg <- "machine"
I.mach_sec <- "all"
I.mach_cvr <- 1073741824 # number of bytes in 1GB (used for conversion of RAM)
I.machine_a <- data.frame(info = "OS",
                          value = sessionInfo()$running,
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_b <- data.frame(info = "processor",
                          value = benchmarkme::get_cpu()$model_name,
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_c <- data.frame(info = "number of cores",
                          value = parallel::detectCores(logical = F),
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_d <- data.frame(info = "number of threads",
                          value = parallel::detectCores(logical = T),
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_e <- data.frame(info = "RAM",
                          value = 
                            paste(
                              as.numeric(benchmarkme::get_ram()) / I.mach_cvr,
                              "GB", sep = ""),
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_0 <- rbind(I.machine_a, I.machine_b, I.machine_c, I.machine_d, 
                     I.machine_e)
I.machine <- dplyr::select(I.machine_0, category, info, value, script, section)

# capture R-related information:
I.rlan_ctg <- "base R"
I.rlan_sec <- "all"
I.baseRpl_a <- data.frame(info = "version",
                          value = R.Version()$version.string,
                          category = I.rlan_ctg, script = name_scrpt,
                          section = I.rlan_sec, stringsAsFactors = F)
I.baseRpl_b <- data.frame(info = "nickname",
                          value = R.Version()$nickname,
                          category = I.rlan_ctg, script = name_scrpt,
                          section = I.rlan_sec, stringsAsFactors = F)
I.baseRpl_c <- data.frame(info = "platform",
                          value = R.Version()$platform,
                          category = I.rlan_ctg, script = name_scrpt,
                          section = I.rlan_sec, stringsAsFactors = F)
I.baseRpl_0 <- rbind(I.baseRpl_a, I.baseRpl_b, I.baseRpl_c)
I.baseRpl <- dplyr::select(I.baseRpl_0, category, info, value, script, section)

# capture RStudio-related information:
I.RStudio <- data.frame(info = "version",
                        value = as.character(RStudio.Version()$version),
                        category = "R Studio", script = name_scrpt,
                        section = "all", stringsAsFactors = F)

# rbind all of the above data.frames together and write outputs
I.info <- rbind(I.project, I.machine, I.pathsto, I.baseRpl, I.RStudio, 
                I.Rpackge)

# outputs to the vault
write.table(sep = "\t", row.names = F, x = I.info, file = I.ofv_info)

I.obj <- ls(pattern = "I.")
I.lst <- c(I.obj[grep(pattern = "I.", x = I.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, "captureVersions")
save(list = I.lst, file = I.ofv_wksp)

### ************************************
### ^^^DESCRIPTIONS FOR EACH SECTION ----
### ************************************

# ^^^ descriptions not complete for sections outside of this script

## Section A - create a master table (EBTKS)
# purpose: take outputs from QIIME 2 and combine them into a master table...
# ... to serve as the entry point for all downstream processing/analysis
# an EBTKS feature table contains ...
# ... hashed FeatureIDs & representative sequences for all ASVs ...
# ... full & truncated Greengenes & SILVA taxonomic lineages for all ASVs ...
# ... absolute counts for all ASVs in each sample of the dataset
# NOTE: EBTKS = Everything But The Kitchen Sink

## Section S - subset the processed EBTKS table to isolate samples of interest
# purpose: subset tables with specific groups are needed for a number of ...
# .. downstream analyses

## Section T - create a rooted phylogenetic tree from processed EBTKS table
# purpose: a rooted tree is needed for phylogenetic-based analyses
# NOTE: be aware that this is a [typically] computationally intensive section

## Section B - beta diversity and microbiota composition
#
## Ssection C - colonization (i.e. ASV/taxon relative abundances)
#
## Section D - differential abundance testing
# purpose:
# this section tests for differentially abundant features (ASVs) b/w two groups
## Section L - tumor outcomes - neoplastic lesions
#

### ************************************
### A - FUNCTIONS ----
### ************************************

# create a vector naming all Section A functions (used when saving workspace)
A.functions <- c("fasta_to_df", "trunc_tax", "testif_fts_lost", 
                 "gather_lost_fts")

# ** note for KDP: fasta_to_df() version 0.2 ** #
# function takes an input .fasta file in two-line format which looks like:
# line 1: >header
# line 2: DNA sequence
# ... and creates a two column dataframe
# fasta_to_df == fasta file to data.frame
fasta_to_df <- function(fasta_file, hdr_colname = "", seq_colname = "") {
  # read in the .fasta file line by line
  fasta <- readLines(fasta_file)
  # identify header lines by finding line numbers with the '>' character
  hlines <- grep(pattern = ">", x = fasta)
  # create a three column df containing:
  # header line numbers (hdr)
  # sequence line number beginning (beg)
  # sequence line number end (end)
  slines <- data.frame(hdr = hlines,
                       beg = hlines+1,
                       end = c((hlines-1)[-1], length(fasta)))
  # create a vector of identical length to hdr_lines to be used for
  # storing values obtained in the loop below
  storage_vec <- rep(NA, length(hlines))
  # loop to obtain sequences
  for(i in 1:length(hlines)) {
    storage_vec[i] <- paste(fasta[slines$beg[i]:slines$end[i]], collapse = "")
  }
  # create a two column df containing:
  # the header with '>' character replaced (hdr)
  # representative sequence associated with the header (seq)
  new_df_v0 <- data.frame("hdr" = gsub(pattern = ">", replacement = "",
                                       x = fasta[hlines]),
                          "seq" = storage_vec)
  # replace generic column names with names specified in the function input
  new_df_v1 <- new_df_v0
  hdr_colnum <- which(names(new_df_v1) == "hdr")
  seq_colnum <- which(names(new_df_v1) == "seq")
  names(new_df_v1)[hdr_colnum] <- hdr_colname
  names(new_df_v1)[seq_colnum] <- seq_colname
  # cols are currently factors, coerce them to characters
  new_df_v2 <- new_df_v1
  new_df_v2[, hdr_colnum] <- as.character(new_df_v2[, hdr_colnum])
  new_df_v2[, seq_colnum] <- as.character(new_df_v2[, seq_colnum])
  return(new_df_v2)
}
# 
# # example usage:
# df_fasta <- fasta_to_df(fasta_file = ifp_fasta,
#                         hdr_colname = "FeatureID",
#                         seq_colname = "RepSeq")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: trunc_tax() version 0.4 ** #
# trunc_tax() function requires an input df with a column containing either:
# greengenes formatted taxonomic strings, e.g.:
# k__; p__ ; c__; o__; f__; g__; s__"
# or silva taxonomic strings, e.g.:
# "D_0__;D_1__;D_2__;D_3__;D_4__;D_5__;D_6__"
# the output df will contain a new column with:
# the desired taxonomic level for the specified database naming convention
# NOTE: if the above strings are not in consecutive order...
# ... the output may not be correct ...
# i.e. if string input is missing D_2__ = "D_0__;D_1__Example;D_3__Example;D_4__"
# output for ;D_1__ = "Example;D_3__Example;D_4__" rather than just "Example"
# output for ;D_2__ = "Unassigned"
# output for ;D_3__ = "Example"
# trunc_tax = truncate taxonomic string
trunc_tax <- function(data = data.frame, tax_type = c("Greengenes", "SILVA"),
                      tax_lvl = c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "Species"),
                      input_col = "", str1 = "", str2 = "", classifier = "") {
  
  # internal checks to ensure correct input classes
  if (!inherits(data, "data.frame")) {
    stop("input data must be class 'data.frame'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA") {
    stop("tax_type should be either 'Greengenes' or 'SILVA'")
  }
  if (!tax_lvl == "Kingdom" && !tax_lvl == "Phylum" && !tax_lvl == "Class" &&
      !tax_lvl == "Order" && !tax_lvl == "Family" && !tax_lvl == "Genus" &&
      !tax_lvl == "Species") {
    stop("tax_lvl should be one of: 'Kingdom', 'Phylum', 'Class', 'Order',
        'Family', 'Genus', 'Species'")
  }
  if (!inherits(input_col, "character")) {
    stop("input_col must be character")
  }
  if (!inherits(str1, "character")) {
    stop("str1 must be character")
  }
  if (!inherits(str2, "character")) {
    stop("str2 must be character")
  }
  if (!inherits(classifier, "character")) {
    stop("classifier must be character")
  }
  
  # if tax_lvl is Species, pasting the string split is handled differently
  num <- ifelse(tax_lvl == "Species", yes = 2, no = 1)
  
  # store original input df and format the new df
  new_df1 <- data
  new_df1$new_col <- 1
  
  # loop through data and truncate the taxonomic strings
  for(row in 1:nrow(new_df1)) {
    # split the string in the specified input column with input of str1
    splt1 <- base::strsplit(new_df1[row, input_col], str1)[[1]][2]
    # goal of first ifelse test:
    # yes = tax level was not assigned
    # no = tax level was likely assigned, split the string with input of str2
    splt2 <- ifelse(is.na(splt1), yes = "Unassigned",
                    no = paste("", base::strsplit(splt1, str2)[[1]][num], 
                               sep = ""))
    # goal of second ifelse test:
    # yes = tax level was not assigned
    # no = tax level was actually assigned
    new_df1[row, "new_col"] <- ifelse(nchar(splt2) == 0 | "NA" %in% splt2,
                                      yes = "Unassigned", no = splt2)
    # if tax_type is Greengenes, tax_lvl is Species and row is not Unassigned,
    # append Genus in front of Species assignment
    if (tax_type == "Greengenes" & tax_lvl == "Species" &
        !new_df1[row, "new_col"] == "Unassigned") {
      new_df1[row, "new_col"] <- paste(base::strsplit(splt1, str2)[[1]][1],
                                       base::strsplit(splt1, str2)[[1]][2], 
                                       sep = " ")
    }
  } # close loop
  
  # internal warning in the event that all rows in the new column are Unassigned
  # this outcome could be correct if the specified tax lvl had no assignemnts;
  # however, if assignments are expected for the specified tax lvl, this outcome
  # is incorrect and the function inputs need to be double checked
  if (isTRUE(all(grepl(pattern = "Unassigned", x = new_df1[, "new_col"],
                       ignore.case = FALSE, fixed = TRUE)))) {
    warning("in new_col all row values in output column are 'Unassigned'")
  }
  
  # format the new df
  new_df2 <- new_df1
  new_col_num <- which(names(new_df2) == 'new_col')
  if (tax_type == "Greengenes") {
    names(new_df2)[new_col_num] <- paste(classifier, "ggs", tax_lvl, sep = ".")
  }
  if (tax_type == "SILVA") {
    names(new_df2)[new_col_num] <- paste(classifier, "slv", tax_lvl, sep = ".")
  }
  return(new_df2)
}
#
# example usage:
# inputs are incorrect, prints warning
#new_df <- trunc_tax(data = a.df, tax_type = "silva",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")

# corrected inputs, no warning
#new_df <- trunc_tax(data = a.df, tax_type = "greengenes",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.1 ** #
# testif_fts_lost() and gather_lost_fts() functions take identical inputs:
# 1: a raw taxonomy data.frame with column names changed (see: A - STEP 1a)
# 2: a processed taxanomy data frame from one of A - Step 1a; 1b; 1c
# 3: the type specificaion for the input to test_df
# 4: the taxonomic database used for the input to test_df
# 5: input column name to join by (e.g. FeatureID or some equivalent)
# NOTE: input data.frames must be from the same taxonomy database
# NOTE: I understand that there are a million ways to do this, I chose one.
# NOTE: both functions are dependent upon package dplyr()
# ** note for KDP: future versions to halt incorrect _type combos ** #

# the testif_fts_lost() function: 
# checks if any features were lost after a processing step was performed
# if features were lost, an error message is printed and ...
# ... an empty data.frame is returned in order to prevent crucial parts ...
# of the downstream code from continuing on
# the gather_lost_fts() function: 
# identifies any lost features after a processing step was performed
# typically, the SILVA database is the issue with features missing when the ...
# _slv_lws variable is created from the _slv_trunc_fix variable in:
# A - STEP 1c - format taxonomy tables: determine lowest assignment level
# this specific error can be overcome by adding to the _slv_trunc_fix variable
# # more details about that process can be found in:
# A - STEP 1b - format taxonomy tables: "fix" truncated lineages

testif_fts_lost <- function(raw_df = data.frame, test_df = data.frame, 
                            test_type = c("trunc", "trunc_fix", "lws", 
                                          "frm", "mrg"),
                            tax_type = c("Greengenes", "SILVA", "merging"),
                            col_join = "") {
  # internal checks to ensure correct input classes
  if (!inherits(raw_df, "data.frame")) {
    stop("input raw_df must be class 'data.frame'")
  }
  if (!inherits(test_df, "data.frame")) {
    stop("input test_df must be class 'data.frame'")
  }
  if (!test_type == "trunc" && !test_type == "trunc_fix" && 
      !test_type == "lws" && !test_type == "frm" && !test_type == "mrg") {
    stop("test_type must be one of: 'trunc', 'trunc_fix', 'lws', 'frm', 'mrg'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA" && 
      !tax_type == "merging") {
    stop("tax_type must be one of: 'Greengenes' or 'SILVA', 'merging'")
  }
  if (!inherits(col_join, "character")) {
    stop("col_join must be character")
  }
  
  # create the error message for !test_type == "mrg" and !tax_type == "merging"
  err_pr1 <- "\ndetected lost fts in '_"
  err_pr2 <- "' data.frame from "
  err <- paste(err_pr1, test_type, err_pr2, tax_type, sep = "")
  
  # create the error message for test_type == "mrg" and tax_type == "merging"
  err_prg <- "' data.frame after database "
  err_mrg <- paste(err_pr1, test_type, err_prg, tax_type, sep = "")
  
  # add the final part to the error message
  err_fpr <- "run: gather_lost_fts() function to identify missing features"
  fin_err <- paste(err, err_fpr, sep = "\n")
  fin_err_mrg <- paste(err_mrg, err_fpr, sep = "\n")
  
  # anti_join the two input data.frames to determine if features were lost
  anti_df <- dplyr::anti_join(x = raw_df, y = test_df, by = col_join)
  
  # using an ifelse statement, create a logical variable to test the anti_df
  LOGICAL_anti_df <- ifelse(nrow(anti_df) == 0, yes = T, no = F)
  
  # test the LOGICAL variable; 
  # if FALSE, return an empty data.frame which will halt downstream code and ...
  # ... print the appropriate error message based upon the input to test_type
  
  if (!test_type == "merging") {  
    if (!isTRUE(LOGICAL_anti_df)) {
      test_df <- data.frame()
      stop(fin_err)
    }
  }
  
  if (test_type == "merging") {  
    if (!isTRUE(LOGICAL_anti_df)) {
      test_df <- data.frame()
      stop(fin_err_mrg)
    }
  }
  return(test_df)
}

gather_lost_fts <- function(raw_df = data.frame, test_df = data.frame, 
                            test_type = c("trunc", "trunc_fix", "lws", 
                                          "frm", "mrg"),
                            tax_type = c("Greengenes", "SILVA", "merging"),
                            col_join = "") {
  # internal checks to ensure correct input classes
  if (!inherits(raw_df, "data.frame")) {
    stop("input raw_df must be class 'data.frame'")
  }
  if (!inherits(test_df, "data.frame")) {
    stop("input test_df must be class 'data.frame'")
  }
  if (!test_type == "trunc" && !test_type == "trunc_fix" && 
      !test_type == "lws" && !test_type == "frm" && !test_type == "mrg") {
    stop("test_type must be one of: 'trunc', 'trunc_fix', 'lws', 'frm', 'mrg'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA" && 
      !tax_type == "merging") {
    stop("tax_type must be one of: 'Greengenes' or 'SILVA', 'merging'")
  }
  if (!inherits(col_join, "character")) {
    stop("col_join must be character")
  }
  
  # anti_join the two input data.frames to determine if features were lost
  anti_df <- dplyr::anti_join(x = raw_df, y = test_df, by = col_join)
  return(anti_df)
}
#
# example usage:
# A.tst_slv_lws <- testif_fts_lost(raw_df = A.int_slv, test_df = A.int_slv_lws, 
#                                   test_type = "lws", tax_type = "SILVA",
#                                   col_join = "FeatureID")
# A.lost_slv_lws <- gather_lost_fts(raw_df = A.int_slv, test_df = A.int_slv_lws, 
#                                   test_type = "lws", tax_type = "SILVA",
#                                   col_join = "FeatureID")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

### ************************************
### A - STEP 1a - format taxonomy tables: truncate lineages ----
### ************************************

# ** note for KDP: AS truncate lineages version 0.4 ** #

# this step truncates Greengenes & SILVA taxonomic lineages, placing them...
# ...into new cols: .Kingdom; .Phylum; .Class; .Order; .Family; .Genus; .Species
# this step also determines the lowest assigned level, placing that info...
# ... into new cols: .lws.txn and .lws.lvl

# read in taxonomy tables
A.raw_int_ggs <- read.table(A.ifv_int_ggs, header = T, sep = "\t", as.is = T,
                            stringsAsFactors = F)
A.raw_int_slv <- read.table(A.ifv_int_slv, header = T, sep = "\t", as.is = T,
                            stringsAsFactors = F, comment.char = "", quote = "")

# rename columns
A.int_ggs <- dplyr::rename(A.raw_int_ggs, FeatureID = Feature.ID,
                           int.ggs.tax = Taxon, int.ggs.cnf = Confidence)
A.int_slv <- dplyr::rename(A.raw_int_slv, FeatureID = Feature.ID,
                           int.slv.tax = Taxon, int.slv.cnf = Confidence)

# truncate taxonomic lineages for each level, and reorder columns
A.int_ggs_L1 <- trunc_tax(data = A.int_ggs, tax_type = "Greengenes",
                          tax_lvl = "Kingdom", input_col = "int.ggs.tax",
                          str1 = "k__", str2 = "; p__", classifier = "int")
A.int_ggs_L2 <- trunc_tax(data = A.int_ggs_L1, tax_type = "Greengenes",
                          tax_lvl = "Phylum", input_col = "int.ggs.tax",
                          str1 = "; p__", str2 = "; c__", classifier = "int")
A.int_ggs_L3 <- trunc_tax(data = A.int_ggs_L2, tax_type = "Greengenes",
                          tax_lvl = "Class", input_col = "int.ggs.tax",
                          str1 = "; c__", str2 = "; o__", classifier = "int")
A.int_ggs_L4 <- trunc_tax(data = A.int_ggs_L3, tax_type = "Greengenes",
                          tax_lvl = "Order", input_col = "int.ggs.tax",
                          str1 = "; o__", str2 = "; f__", classifier = "int")
A.int_ggs_L5 <- trunc_tax(data = A.int_ggs_L4, tax_type = "Greengenes",
                          tax_lvl = "Family", input_col = "int.ggs.tax",
                          str1 = "; f__", str2 = "; g__", classifier = "int")
A.int_ggs_L6 <- trunc_tax(data = A.int_ggs_L5, tax_type = "Greengenes",
                          tax_lvl = "Genus", input_col = "int.ggs.tax",
                          str1 = "; g__", str2 = "; s__", classifier = "int")
A.int_ggs_L7 <- trunc_tax(data = A.int_ggs_L6, tax_type = "Greengenes",
                          tax_lvl = "Species", input_col = "int.ggs.tax",
                          str1 = "; g__", str2 = "; s__", classifier = "int")
A.int_ggs_trunc <- dplyr::select(A.int_ggs_L7, FeatureID, int.ggs.tax,
                                 int.ggs.Kingdom, int.ggs.Phylum, int.ggs.Class,
                                 int.ggs.Order, int.ggs.Family, int.ggs.Genus,
                                 int.ggs.Species, int.ggs.cnf)

A.int_slv_L1 <- trunc_tax(data = A.int_slv, tax_type = "SILVA",
                          tax_lvl = "Kingdom", input_col = "int.slv.tax",
                          str1 = "D_0__", str2 = ";D_1__", classifier = "int")
A.int_slv_L2 <- trunc_tax(data = A.int_slv_L1, tax_type = "SILVA",
                          tax_lvl = "Phylum", input_col = "int.slv.tax",
                          str1 = ";D_1__", str2 = ";D_2__", classifier = "int")
A.int_slv_L3 <- trunc_tax(data = A.int_slv_L2, tax_type = "SILVA",
                          tax_lvl = "Class", input_col = "int.slv.tax",
                          str1 = ";D_2__", str2 = ";D_3__", classifier = "int")
A.int_slv_L4 <- trunc_tax(data = A.int_slv_L3, tax_type = "SILVA",
                          tax_lvl = "Order", input_col = "int.slv.tax",
                          str1 = ";D_3__", str2 = ";D_4__", classifier = "int")
A.int_slv_L5 <- trunc_tax(data = A.int_slv_L4, tax_type = "SILVA",
                          tax_lvl = "Family", input_col = "int.slv.tax",
                          str1 = ";D_4__", str2 = ";D_5__", classifier = "int")
A.int_slv_L6 <- trunc_tax(data = A.int_slv_L5, tax_type = "SILVA",
                          tax_lvl = "Genus", input_col = "int.slv.tax",
                          str1 = ";D_5__", str2 = ";D_6__", classifier = "int")
A.int_slv_L7 <- trunc_tax(data = A.int_slv_L6, tax_type = "SILVA",
                          tax_lvl = "Species", input_col = "int.slv.tax",
                          str1 = ";D_5__", str2 = ";D_6__", classifier = "int")
A.int_slv_trunc <- dplyr::select(A.int_slv_L7, FeatureID, int.slv.tax,
                                 int.slv.Kingdom, int.slv.Phylum, int.slv.Class,
                                 int.slv.Order, int.slv.Family, int.slv.Genus,
                                 int.slv.Species, int.slv.cnf)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost in either of the truncated lineage dfs

A.test_ggs_trunc <- testif_fts_lost(raw_df = A.int_ggs, 
                                    test_df = A.int_ggs_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "Greengenes",
                                    col_join = "FeatureID")
A.test_slv_trunc <- testif_fts_lost(raw_df = A.int_slv, 
                                    test_df = A.int_slv_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "SILVA",
                                    col_join = "FeatureID")
A.lost_ggs_trunc <- gather_lost_fts(raw_df = A.int_ggs, 
                                    test_df = A.int_ggs_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "Greengenes",
                                    col_join = "FeatureID")
A.lost_slv_trunc <- gather_lost_fts(raw_df = A.int_slv, 
                                    test_df = A.int_slv_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "SILVA",
                                    col_join = "FeatureID")

### ************************************
### A - STEP 1b - format taxonomy tables: "fix" truncated lineages ----
### ************************************

# ** note for KDP: AS "fix" truncated lineages version 0.9 ** #

# An aside on the subject of truncated taxonomic lineages:

# truncated taxonomic lineages are great and mostly useful; however...
# ... a variety of SILVA assignments lose useful information when truncated
# e.g. Order = Lineage IV
# e.g. Family = uncultured Rubrobacteria bacterium
# e.g. Genus = AAP99
# this is especially problematic for the SILVA assignment Family = Family XI ...
# ... this family may be in the Order Bacillales or the Order Clostridiales
# the goal here is to add info to make truncated assignments useful/meaningful
# e.g. Family XI has info prepended for the Order
# e.g. Order: Lineage IV is changed to Elusimicrobia Lineage IV
# NOTE: the primary focus here is on SILVA as this database is actually updated
# NOTE: these are "running-lists" and are likely dataset specific...
# ... i.e. not all lineages will always be present in every dataset
# NOTE: dtabase version info is located in COMPUTER/VERSION/PACKAGE INFO

# end of aside, begin comments for the code below

# create vectors to be used to index column names
# NOTE: these vectors will be indexed using numbers; DO NOT ALTER THEIR ORDER
A.ggs_col <- c("int.ggs.Kingdom", "int.ggs.Phylum", "int.ggs.Class",
               "int.ggs.Order", "int.ggs.Family","int.ggs.Genus",
               "int.ggs.Species")

A.slv_col <- c("int.slv.Kingdom", "int.slv.Phylum", "int.slv.Class",
               "int.slv.Order", "int.slv.Family","int.slv.Genus",
               "int.slv.Species")

# the below is performed on one assignment at a time; the process is as follows:
# (1) create a vector to filter the lineage that needs "fixed"
# (2) isolate that lineage from the A.test_slv_trunc df created in A - STEP 1A
# (3) fix the lineage
# (4) rinse and repeat for all lineages that need to be "fixed"
# NOTE: the capital letter combinations indicate taxonomic level...
# ... i.e. PC = Phylum-Class; OFG = Order-Family-Genus, etc.

# single level fixes first:
# Phylum-Class
A.slv_vec_PC_01 <- c("Acidobacteria", "Subgroup 6")
A.slv_flt_PC_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[2]]] == A.slv_vec_PC_01[1],
                                 .data[[A.slv_col[3]]] == A.slv_vec_PC_01[2])
A.slv_fix_PC_01 <- dplyr::mutate(
  A.slv_flt_PC_01,
  int.slv.Class = replace(x = .data[[A.slv_col[2]]],
                          list = .data[[A.slv_col[3]]] == A.slv_vec_PC_01[1] &
                            .data[[A.slv_col[2]]] == A.slv_vec_PC_01[2],
                          values = paste(A.slv_vec_PC_01[1],
                                         A.slv_vec_PC_01[2],
                                         sep = " ")))

# Class-Order
A.slv_vec_CO_01 <- c("Anaerolineae", "SBR1031")
A.slv_flt_CO_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[3]]] == A.slv_vec_CO_01[1],
                                 .data[[A.slv_col[4]]] == A.slv_vec_CO_01[2])
A.slv_fix_CO_01 <- dplyr::mutate(
  A.slv_flt_CO_01,
  int.slv.Class = replace(x = .data[[A.slv_col[3]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_CO_01[1] &
                            .data[[A.slv_col[4]]] == A.slv_vec_CO_01[2],
                          values = paste(A.slv_vec_CO_01[1],
                                         A.slv_vec_CO_01[2],
                                         sep = " ")))

A.slv_vec_CO_02 <- c("Dehalococcoidia", "S085")
A.slv_flt_CO_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[3]]] == A.slv_vec_CO_02[1],
                                 .data[[A.slv_col[4]]] == A.slv_vec_CO_02[2])
A.slv_fix_CO_02 <- dplyr::mutate(
  A.slv_flt_CO_02,
  int.slv.Class = replace(x = .data[[A.slv_col[3]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_CO_02[1] &
                            .data[[A.slv_col[4]]] == A.slv_vec_CO_02[2],
                          values = paste(A.slv_vec_CO_02[1],
                                         A.slv_vec_CO_02[2],
                                         sep = " ")))

A.slv_vec_CO_03 <- c("Elusimicrobia", "Lineage IV")
A.slv_flt_CO_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[3]]] == A.slv_vec_CO_03[1],
                                 .data[[A.slv_col[4]]] == A.slv_vec_CO_03[2])
A.slv_fix_CO_03 <- dplyr::mutate(
  A.slv_flt_CO_03,
  int.slv.Class = replace(x = .data[[A.slv_col[3]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_CO_03[1] &
                            .data[[A.slv_col[4]]] == A.slv_vec_CO_03[2],
                          values = paste(A.slv_vec_CO_03[1],
                                         A.slv_vec_CO_03[2],
                                         sep = " ")))

# Order-Family
A.slv_vec_OF_01 <- c("Bacillales", "Family XI")
A.slv_flt_OF_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_01[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_01[2])
A.slv_fix_OF_01 <- dplyr::mutate(
  A.slv_flt_OF_01,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_01[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_01[2],
                           values = paste(A.slv_vec_OF_01[1],
                                          A.slv_vec_OF_01[2],
                                          sep = " ")))

A.slv_vec_OF_02 <- c("Bacteroidales", "p-2534-18B5 gut group")
A.slv_flt_OF_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_02[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_02[2])
A.slv_fix_OF_02 <- dplyr::mutate(
  A.slv_flt_OF_02,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_02[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_02[2],
                           values = paste(A.slv_vec_OF_02[1],
                                          A.slv_vec_OF_02[2],
                                          sep = " ")))

A.slv_vec_OF_03 <- c("Clostridiales", "Family XI")
A.slv_flt_OF_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_03[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_03[2])
A.slv_fix_OF_03 <- dplyr::mutate(
  A.slv_flt_OF_03,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_03[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_03[2],
                           values = paste(A.slv_vec_OF_03[1],
                                          A.slv_vec_OF_03[2],
                                          sep = " ")))

A.slv_vec_OF_04 <- c("Flavobacteriales", "NS9 marine group")
A.slv_flt_OF_04 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_04[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_04[2])
A.slv_fix_OF_04 <- dplyr::mutate(
  A.slv_flt_OF_04,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_04[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_04[2],
                           values = paste(A.slv_vec_OF_04[1],
                                          A.slv_vec_OF_04[2],
                                          sep = " ")))

A.slv_vec_OF_05 <- c("Kryptoniales", "BSV26")
A.slv_flt_OF_05 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_05[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_05[2])
A.slv_fix_OF_05 <- dplyr::mutate(
  A.slv_flt_OF_05,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_05[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_05[2],
                           values = paste(A.slv_vec_OF_05[1],
                                          A.slv_vec_OF_05[2],
                                          sep = " ")))

A.slv_vec_OF_06 <- c("Myxococcales", "BIrii41")
A.slv_flt_OF_06 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_06[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_06[2])
A.slv_fix_OF_06 <- dplyr::mutate(
  A.slv_flt_OF_06,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_06[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_06[2],
                           values = paste(A.slv_vec_OF_06[1],
                                          A.slv_vec_OF_06[2],
                                          sep = " ")))

A.slv_vec_OF_07 <- c("Solirubrobacterales", "67-14")
A.slv_flt_OF_07 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_07[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_07[2])
A.slv_fix_OF_07 <- dplyr::mutate(
  A.slv_flt_OF_07,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_07[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_07[2],
                           values = paste(A.slv_vec_OF_07[1],
                                          A.slv_vec_OF_07[2],
                                          sep = " ")))

A.slv_vec_OF_08 <- c("Thermoanaerobacterales", "Family III")
A.slv_flt_OF_08 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_08[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_08[2])
A.slv_fix_OF_08 <- dplyr::mutate(
  A.slv_flt_OF_08,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_08[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_08[2],
                           values = paste(A.slv_vec_OF_08[1],
                                          A.slv_vec_OF_08[2],
                                          sep = " ")))

A.slv_vec_OF_09 <- c("Thermomicrobiales", "JG30-KF-CM45")
A.slv_flt_OF_09 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_09[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_09[2])
A.slv_fix_OF_09 <- dplyr::mutate(
  A.slv_flt_OF_09,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_09[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_09[2],
                           values = paste(A.slv_vec_OF_09[1],
                                          A.slv_vec_OF_09[2],
                                          sep = " ")))

A.slv_vec_OF_10 <- c("Clostridiales", "Family XIII")
A.slv_flt_OF_10 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_10[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_10[2],
                                 !grepl(pattern = A.slv_vec_OF_10[2],
                                        x = .data[[A.slv_col[6]]]))
A.slv_fix_OF_10 <- dplyr::mutate(
  A.slv_flt_OF_10,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_10[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_10[2],
                           values = paste(A.slv_vec_OF_10[1],
                                          A.slv_vec_OF_10[2],
                                          sep = " ")))

A.slv_vec_OF_11 <- c("Verrucomicrobiales", "DEV007")
A.slv_flt_OF_11 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_11[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_11[2],
                                 !grepl(pattern = A.slv_vec_OF_11[2],
                                        x = .data[[A.slv_col[6]]]))
A.slv_fix_OF_11 <- dplyr::mutate(
  A.slv_flt_OF_11,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_11[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_11[2],
                           values = paste(A.slv_vec_OF_11[1],
                                          A.slv_vec_OF_11[2],
                                          sep = " ")))

A.slv_vec_OF_12 <- c("Myxococcales", "P3OB-42")
A.slv_flt_OF_12 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_12[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_12[2],
                                 !grepl(pattern = A.slv_vec_OF_12[2],
                                        x = .data[[A.slv_col[6]]]))
A.slv_fix_OF_12 <- dplyr::mutate(
  A.slv_flt_OF_12,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_12[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_12[2],
                           values = paste(A.slv_vec_OF_12[1],
                                          A.slv_vec_OF_12[2],
                                          sep = " ")))

# Family-Genus
A.slv_vec_FG_01 <- c("Actinomycetaceae", "F0332")
A.slv_flt_FG_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_01[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_01[2])
A.slv_fix_FG_01 <- dplyr::mutate(
  A.slv_flt_FG_01,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_01[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_01[2],
                          values = paste(A.slv_vec_FG_01[1],
                                         A.slv_vec_FG_01[2],
                                         sep = " ")))

A.slv_vec_FG_02 <- c("Burkholderiaceae", "AAP99")
A.slv_flt_FG_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_02[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_02[2])
A.slv_fix_FG_02 <- dplyr::mutate(
  A.slv_flt_FG_02,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_02[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_02[2],
                          values = paste(A.slv_vec_FG_02[1],
                                         A.slv_vec_FG_02[2],
                                         sep = " ")))

A.slv_vec_FG_03 <- c("Lachnospiraceae", "A2")
A.slv_flt_FG_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_03[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_03[2])
A.slv_fix_FG_03 <- dplyr::mutate(
  A.slv_flt_FG_03,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_03[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_03[2],
                          values = paste(A.slv_vec_FG_03[1],
                                         A.slv_vec_FG_03[2],
                                         sep = " ")))

A.slv_vec_FG_04 <- c("Lachnospiraceae", "ASF356")
A.slv_flt_FG_04 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_04[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_04[2])
A.slv_fix_FG_04 <- dplyr::mutate(
  A.slv_flt_FG_04,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_04[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_04[2],
                          values = paste(A.slv_vec_FG_04[1],
                                         A.slv_vec_FG_04[2],
                                         sep = " ")))

A.slv_vec_FG_05 <- c("Lachnospiraceae", "CAG-56")
A.slv_flt_FG_05 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_05[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_05[2])
A.slv_fix_FG_05 <- dplyr::mutate(
  A.slv_flt_FG_05,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_05[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_05[2],
                          values = paste(A.slv_vec_FG_05[1],
                                         A.slv_vec_FG_05[2],
                                         sep = " ")))

A.slv_vec_FG_06 <- c("Lachnospiraceae", "GCA-900066575")
A.slv_flt_FG_06 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_06[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_06[2])
A.slv_fix_FG_06 <- dplyr::mutate(
  A.slv_flt_FG_06,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_06[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_06[2],
                          values = paste(A.slv_vec_FG_06[1],
                                         A.slv_vec_FG_06[2],
                                         sep = " ")))

A.slv_vec_FG_07 <- c("Lachnospiraceae", "GCA-900066755")
A.slv_flt_FG_07 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_07[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_07[2])
A.slv_fix_FG_07 <- dplyr::mutate(
  A.slv_flt_FG_07,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_07[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_07[2],
                          values = paste(A.slv_vec_FG_07[1],
                                         A.slv_vec_FG_07[2],
                                         sep = " ")))

A.slv_vec_FG_08 <- c("Lachnospiraceae", "UC5-1-2E3")
A.slv_flt_FG_08 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_08[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_08[2])
A.slv_fix_FG_08 <- dplyr::mutate(
  A.slv_flt_FG_08,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_08[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_08[2],
                          values = paste(A.slv_vec_FG_08[1],
                                         A.slv_vec_FG_08[2],
                                         sep = " ")))

A.slv_vec_FG_09 <- c("Pseudohongiellaceae", "BIyi10")
A.slv_flt_FG_09 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_09[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_09[2])
A.slv_fix_FG_09 <- dplyr::mutate(
  A.slv_flt_FG_09,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_09[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_09[2],
                          values = paste(A.slv_vec_FG_09[1],
                                         A.slv_vec_FG_09[2],
                                         sep = " ")))

A.slv_vec_FG_10 <- c("Ruminococcaceae", "DTU089")
A.slv_flt_FG_10 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_10[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_10[2])
A.slv_fix_FG_10 <- dplyr::mutate(
  A.slv_flt_FG_10,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_10[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_10[2],
                          values = paste(A.slv_vec_FG_10[1],
                                         A.slv_vec_FG_10[2],
                                         sep = " ")))

A.slv_vec_FG_11 <- c("Ruminococcaceae", "GCA-900066225")
A.slv_flt_FG_11 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_11[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_11[2])
A.slv_fix_FG_11 <- dplyr::mutate(
  A.slv_flt_FG_11,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_11[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_11[2],
                          values = paste(A.slv_vec_FG_11[1],
                                         A.slv_vec_FG_11[2],
                                         sep = " ")))

A.slv_vec_FG_12 <- c("Ruminococcaceae", "UBA1819")
A.slv_flt_FG_12 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_12[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_12[2])
A.slv_fix_FG_12 <- dplyr::mutate(
  A.slv_flt_FG_12,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_12[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_12[2],
                          values = paste(A.slv_vec_FG_12[1],
                                         A.slv_vec_FG_12[2],
                                         sep = " ")))

A.slv_vec_FG_13 <- c("Thermoanaerobaculaceae", "Subgroup 10")
A.slv_flt_FG_13 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_13[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_13[2])
A.slv_fix_FG_13 <- dplyr::mutate(
  A.slv_flt_FG_13,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_13[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_13[2],
                          values = paste(A.slv_vec_FG_13[1],
                                         A.slv_vec_FG_13[2],
                                         sep = " ")))

A.slv_vec_FG_14 <- c("Thermoanaerobaculaceae", "Subgroup 23")
A.slv_flt_FG_14 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_14[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_14[2])
A.slv_fix_FG_14 <- dplyr::mutate(
  A.slv_flt_FG_14,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_14[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_14[2],
                          values = paste(A.slv_vec_FG_14[1],
                                         A.slv_vec_FG_14[2],
                                         sep = " ")))

A.slv_vec_FG_15 <- c("Xanthomonadaceae", "SN8")
A.slv_flt_FG_15 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_15[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_15[2])
A.slv_fix_FG_15 <- dplyr::mutate(
  A.slv_flt_FG_15,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_15[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_15[2],
                          values = paste(A.slv_vec_FG_15[1],
                                         A.slv_vec_FG_15[2],
                                         sep = " ")))

A.slv_vec_FG_16 <- c("Ruminococcaceae", "CAG-352")
A.slv_flt_FG_16 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_16[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_16[2])
A.slv_fix_FG_16 <- dplyr::mutate(
  A.slv_flt_FG_16,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_16[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_16[2],
                          values = paste(A.slv_vec_FG_16[1],
                                         A.slv_vec_FG_16[2],
                                         sep = " ")))

# Genus-Species
A.slv_vec_GS_01 <- c("Bacteroides", "bacterium NLAE-zl-H46")
A.slv_flt_GS_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_01[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_01[2])
A.slv_fix_GS_01 <- dplyr::mutate(
  A.slv_flt_GS_01,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_01[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_01[2],
                            values = paste(A.slv_vec_GS_01[1],
                                           A.slv_vec_GS_01[2],
                                           sep = " ")))

A.slv_vec_GS_02 <- c("Christensenellaceae R-7 group", "bacterium YE57")
A.slv_flt_GS_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_02[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_02[2])
A.slv_fix_GS_02 <- dplyr::mutate(
  A.slv_flt_GS_02,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_02[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_02[2],
                            values = paste(A.slv_vec_GS_02[1],
                                           A.slv_vec_GS_02[2],
                                           sep = " ")))

A.slv_vec_GS_03 <- c("Desulfovibrio", "bacterium New Zealand D")
A.slv_flt_GS_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_03[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_03[2])
A.slv_fix_GS_03 <- dplyr::mutate(
  A.slv_flt_GS_03,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_03[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_03[2],
                            values = paste(A.slv_vec_GS_03[1],
                                           A.slv_vec_GS_03[2],
                                           sep = " ")))

A.slv_vec_GS_04 <- c("Gaiella", "actinobacterium WWH12")
A.slv_flt_GS_04 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_04[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_04[2])
A.slv_fix_GS_04 <- dplyr::mutate(
  A.slv_flt_GS_04,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_04[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_04[2],
                            values = paste(A.slv_vec_GS_04[1],
                                           A.slv_vec_GS_04[2],
                                           sep = " ")))

A.slv_vec_GS_05 <- c("Ruminiclostridium 9",
                     "bacterium enrichment culture clone M244")
A.slv_flt_GS_05 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_05[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_05[2])
A.slv_fix_GS_05 <- dplyr::mutate(
  A.slv_flt_GS_05,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_05[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_05[2],
                            values = paste(A.slv_vec_GS_05[1],
                                           A.slv_vec_GS_05[2],
                                           sep = " ")))

A.slv_vec_GS_06 <- c("Lachnospiraceae", "uncultured",
                     "Frisingicoccus caecimuris", "Frisingicoccus")
A.slv_flt_GS_06 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_06[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_06[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_06[3])
A.slv_fix_GS_06 <- dplyr::mutate(
  A.slv_flt_GS_06,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_06[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_06[2],
                          values = A.slv_vec_GS_06[4]))

A.slv_vec_GS_07 <- c("Ruminococcaceae", "uncultured",
                     "Clostridium phoceensis", "Clostridium")
A.slv_flt_GS_07 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_07[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_07[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_07[3])
A.slv_fix_GS_07 <- dplyr::mutate(
  A.slv_flt_GS_07,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_07[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_07[2],
                          values = A.slv_vec_GS_07[4]))

A.slv_vec_GS_08 <- c("Lachnospiraceae", "uncultured",
                     "Ruminococcus sp. ID1", "Ruminococcus")
A.slv_flt_GS_08 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_08[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_08[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_08[3])
A.slv_fix_GS_08 <- dplyr::mutate(
  A.slv_flt_GS_08,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_08[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_08[2],
                          values = A.slv_vec_GS_08[4]))

A.slv_vec_GS_09 <- c("Ruminococcaceae", "uncultured",
                     "Ruminococcaceae bacterium Marseille-P3738")
A.slv_flt_GS_09 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_09[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_09[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_09[3])
A.slv_fix_GS_09 <- dplyr::mutate(
  A.slv_flt_GS_09,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_09[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_09[3],
                            values = A.slv_vec_GS_09[2]))

A.slv_vec_GS_10 <- c("Lachnospiraceae", "uncultured",
                     "butyrate-producing bacterium L2-10")
A.slv_flt_GS_10 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_10[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_10[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_10[3])
A.slv_fix_GS_10 <- dplyr::mutate(
  A.slv_flt_GS_10,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_10[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_10[3],
                            values = A.slv_vec_GS_10[2]))

A.slv_vec_GS_11 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Culture-41", "Clostridium")
A.slv_flt_GS_11 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_11[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_11[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_11[3])
A.slv_fix_GS_11 <- dplyr::mutate(
  A.slv_flt_GS_11,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_11[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_11[2],
                          values = A.slv_vec_GS_11[4]))

A.slv_vec_GS_12 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Culture-54", "Clostridium")
A.slv_flt_GS_12 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_12[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_12[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_12[3])
A.slv_fix_GS_12 <- dplyr::mutate(
  A.slv_flt_GS_12,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_12[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_12[2],
                          values = A.slv_vec_GS_12[4]))

A.slv_vec_GS_13 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Clone-49", "Clostridium")
A.slv_flt_GS_13 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_13[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_13[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_13[3])
A.slv_fix_GS_13 <- dplyr::mutate(
  A.slv_flt_GS_13,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_13[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_13[2],
                          values = A.slv_vec_GS_13[4]))

A.slv_vec_GS_14 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Culture-27", "Clostridium")
A.slv_flt_GS_14 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_14[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_14[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_14[3])
A.slv_fix_GS_14 <- dplyr::mutate(
  A.slv_flt_GS_14,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_14[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_14[2],
                          values = A.slv_vec_GS_14[4]))

A.slv_vec_GS_15 <- c("Barnesiellaceae", "uncultured",
                     "Bacteroides sp. Tilapia9", "Bacteroides")
A.slv_flt_GS_15 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_15[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_15[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_15[3])
A.slv_fix_GS_15 <- dplyr::mutate(
  A.slv_flt_GS_15,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_15[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_15[2],
                          values = A.slv_vec_GS_15[4]))

A.slv_vec_GS_16 <- c("Sphingobacteriaceae", "uncultured",
                     "Sphingobacterium jejuense", "Sphingobacterium")
A.slv_flt_GS_16 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_16[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_16[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_16[3])
A.slv_fix_GS_16 <- dplyr::mutate(
  A.slv_flt_GS_16,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_16[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_16[2],
                          values = A.slv_vec_GS_16[4]))

A.slv_vec_GS_17 <- c("Burkholderiaceae", "uncultured",
                     "Burkholderiales bacterium X4")
A.slv_flt_GS_17 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_17[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_17[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_17[3])
A.slv_fix_GS_17 <- dplyr::mutate(
  A.slv_flt_GS_17,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_17[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_17[3],
                            values = A.slv_vec_GS_17[2]))

A.slv_vec_GS_18 <- c("Burkholderiaceae", "uncultured",
                     "beta proteobacterium WX163")
A.slv_flt_GS_18 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_18[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_18[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_18[3])
A.slv_fix_GS_18 <- dplyr::mutate(
  A.slv_flt_GS_18,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_18[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_18[3],
                            values = A.slv_vec_GS_18[2]))

# Order-Family-Genus
A.slv_vec_OFG_01 <- c("Clostridiales", "Family XIII",
                      "Family XIII AD3011 group")
A.slv_flt_OFG_01 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[4]]] == A.slv_vec_OFG_01[1],
                                  .data[[A.slv_col[5]]] == A.slv_vec_OFG_01[2],
                                  .data[[A.slv_col[6]]] == A.slv_vec_OFG_01[3])
A.slv_fix_OFG_01 <- dplyr::mutate(
  A.slv_flt_OFG_01,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_01[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OFG_01[2],
                           values = paste(A.slv_vec_OFG_01[1],
                                          A.slv_vec_OFG_01[2],
                                          sep = " ")),
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_01[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_OFG_01[3],
                          values = paste(A.slv_vec_OFG_01[1],
                                         A.slv_vec_OFG_01[3],
                                         sep = " ")))

A.slv_vec_OFG_02 <- c("Clostridiales", "Family XIII",
                      "Family XIII UCG-001")
A.slv_flt_OFG_02 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[4]]] == A.slv_vec_OFG_02[1],
                                  .data[[A.slv_col[5]]] == A.slv_vec_OFG_02[2],
                                  .data[[A.slv_col[6]]] == A.slv_vec_OFG_02[3])
A.slv_fix_OFG_02 <- dplyr::mutate(
  A.slv_flt_OFG_02,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_02[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OFG_02[2],
                           values = paste(A.slv_vec_OFG_02[1],
                                          A.slv_vec_OFG_02[2],
                                          sep = " ")),
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_02[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_OFG_02[3],
                          values = paste(A.slv_vec_OFG_02[1],
                                         A.slv_vec_OFG_02[3],
                                         sep = " ")))

A.slv_vec_OFG_03 <- c("Rhodospirillales", "uncultured",
                      "Azospirillum sp. 47_25", "Rhodospirillaceae",
                      "Azospirillum")
A.slv_flt_OFG_03 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[4]]] == A.slv_vec_OFG_03[1],
                                  .data[[A.slv_col[5]]] == A.slv_vec_OFG_03[2],
                                  .data[[A.slv_col[6]]] == A.slv_vec_OFG_03[3])
A.slv_fix_OFG_03 <- dplyr::mutate(
  A.slv_flt_OFG_03,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_03[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OFG_03[2],
                           values = A.slv_vec_OFG_03[4]),
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_03[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_OFG_03[3],
                          values = A.slv_vec_OFG_03[5]))

# Family-Genus-Species
A.slv_vec_FGS_01 <- c("Lachnospiraceae", "uncultured",
                      "intestinal bacterium CG19-1", "CG19-1")
A.slv_flt_FGS_01 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[5]]] == A.slv_vec_FGS_01[1],
                                  .data[[A.slv_col[6]]] == A.slv_vec_FGS_01[2],
                                  .data[[A.slv_col[7]]] == A.slv_vec_FGS_01[3])
A.slv_fix_FGS_01 <- dplyr::mutate(
  A.slv_flt_FGS_01,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FGS_01[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FGS_01[2],
                          values = paste(A.slv_vec_FGS_01[1],
                                         A.slv_vec_FGS_01[4],
                                         sep = " ")),
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list =
                              .data[[A.slv_col[5]]] == A.slv_vec_FGS_01[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_FGS_01[3],
                            values = A.slv_vec_FGS_01[2]))

# combine the _flt dfs within groups and also combine the _fix dfs within groups
A.slv_flt_PC_all <- rbind(A.slv_flt_PC_01)
A.slv_fix_PC_all <- rbind(A.slv_fix_PC_01)

A.slv_flt_CO_all <- rbind(A.slv_flt_CO_01, A.slv_flt_CO_02, A.slv_flt_CO_03)
A.slv_fix_CO_all <- rbind(A.slv_fix_CO_01, A.slv_fix_CO_02, A.slv_fix_CO_03)

A.slv_flt_OF_all <- rbind(A.slv_flt_OF_01, A.slv_flt_OF_02, A.slv_flt_OF_03,
                          A.slv_flt_OF_04, A.slv_flt_OF_05, A.slv_flt_OF_06,
                          A.slv_flt_OF_07, A.slv_flt_OF_08, A.slv_flt_OF_09,
                          A.slv_flt_OF_10, A.slv_flt_OF_11, A.slv_flt_OF_12)
A.slv_fix_OF_all <- rbind(A.slv_fix_OF_01, A.slv_fix_OF_02, A.slv_fix_OF_03,
                          A.slv_fix_OF_04, A.slv_fix_OF_05, A.slv_fix_OF_06,
                          A.slv_fix_OF_07, A.slv_fix_OF_08, A.slv_fix_OF_09,
                          A.slv_fix_OF_10, A.slv_fix_OF_11, A.slv_fix_OF_12)

A.slv_flt_FG_all <- rbind(A.slv_flt_FG_01, A.slv_flt_FG_02, A.slv_flt_FG_03,
                          A.slv_flt_FG_04, A.slv_flt_FG_05, A.slv_flt_FG_06,
                          A.slv_flt_FG_07, A.slv_flt_FG_08, A.slv_flt_FG_09,
                          A.slv_flt_FG_10, A.slv_flt_FG_11, A.slv_flt_FG_12,
                          A.slv_flt_FG_13, A.slv_flt_FG_14, A.slv_flt_FG_15,
                          A.slv_flt_FG_16)
A.slv_fix_FG_all <- rbind(A.slv_fix_FG_01, A.slv_fix_FG_02, A.slv_fix_FG_03,
                          A.slv_fix_FG_04, A.slv_fix_FG_05, A.slv_fix_FG_06,
                          A.slv_fix_FG_07, A.slv_fix_FG_08, A.slv_fix_FG_09,
                          A.slv_fix_FG_10, A.slv_fix_FG_11, A.slv_fix_FG_12,
                          A.slv_fix_FG_13, A.slv_fix_FG_14, A.slv_fix_FG_15,
                          A.slv_fix_FG_16)

A.slv_flt_GS_all <- rbind(A.slv_flt_GS_01, A.slv_flt_GS_02, A.slv_flt_GS_03,
                          A.slv_flt_GS_04, A.slv_flt_GS_05, A.slv_flt_GS_06,
                          A.slv_flt_GS_07, A.slv_flt_GS_08, A.slv_flt_GS_09,
                          A.slv_flt_GS_10, A.slv_flt_GS_11, A.slv_flt_GS_12,
                          A.slv_flt_GS_13, A.slv_flt_GS_14, A.slv_flt_GS_15,
                          A.slv_flt_GS_16, A.slv_flt_GS_17, A.slv_flt_GS_18)
A.slv_fix_GS_all <- rbind(A.slv_fix_GS_01, A.slv_fix_GS_02, A.slv_fix_GS_03,
                          A.slv_fix_GS_04, A.slv_fix_GS_05, A.slv_fix_GS_06,
                          A.slv_fix_GS_07, A.slv_fix_GS_08, A.slv_fix_GS_09,
                          A.slv_fix_GS_10, A.slv_fix_GS_11, A.slv_fix_GS_12,
                          A.slv_fix_GS_13, A.slv_fix_GS_14, A.slv_fix_GS_15,
                          A.slv_fix_GS_16, A.slv_fix_GS_17, A.slv_fix_GS_18)

A.slv_flt_OFG_all <- rbind(A.slv_flt_OFG_01, A.slv_flt_OFG_02, A.slv_flt_OFG_03)
A.slv_fix_OFG_all <- rbind(A.slv_fix_OFG_01, A.slv_fix_OFG_02, A.slv_fix_OFG_03)

A.slv_flt_FGS_all <- rbind(A.slv_flt_FGS_01)
A.slv_fix_FGS_all <- rbind(A.slv_fix_FGS_01)

# then combine all groups together
A.slv_flt_all <- rbind(A.slv_flt_PC_all, A.slv_flt_CO_all, A.slv_flt_OF_all,
                       A.slv_flt_FG_all, A.slv_flt_GS_all, A.slv_flt_OFG_all,
                       A.slv_flt_FGS_all)

A.slv_fix_all <- rbind(A.slv_fix_PC_all, A.slv_fix_CO_all, A.slv_fix_OF_all,
                       A.slv_fix_FG_all, A.slv_fix_GS_all, A.slv_fix_OFG_all,
                       A.slv_fix_FGS_all)

# filter to remove the unfixed lineages and then add in the fixed lineages
A.int_slv_trunc_flt <- dplyr::anti_join(x = A.test_slv_trunc, y = A.slv_fix_all,
                                        by = "FeatureID")

A.int_slv_trunc_fix <- rbind(A.int_slv_trunc_flt, A.slv_fix_all)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost after "fixing" truncated lineages
A.test_slv_trunc_fix <- testif_fts_lost(raw_df = A.int_slv,
                                        test_df = A.int_slv_trunc_fix,
                                        test_type = "trunc_fix",
                                        tax_type = "SILVA",
                                        col_join = "FeatureID")
A.lost_slv_trunc_fix <- gather_lost_fts(raw_df = A.int_slv,
                                        test_df = A.int_slv_trunc_fix,
                                        test_type = "trunc",
                                        tax_type = "SILVA",
                                        col_join = "FeatureID")

# check to make sure there were no NAs produced in the fixed truncated df
# ** note for KDP: possible addition to future version of testif_fts_lost() ** #
A.int_slv_trunc_test_NAs <- dplyr::filter_all(A.int_slv_trunc_fix,
                                              dplyr::any_vars(is.na(.)))

LOGICAL_A.int_slv_trunc_test_NAs <- ifelse(nrow(A.int_slv_trunc_test_NAs) == 0,
                                           yes = T, no = F)

if (!isTRUE(LOGICAL_A.int_slv_trunc_test_NAs)) {
  stop("NA values produced during fix for SILVA truncated lineages;
       run: `View(A.int_slv_trunc_test_NAs)` to see rows with NA values")
}

### ************************************
### A - STEP 1c - format taxonomy tables: determine lowest assignment level ----
### ************************************

# ** note for KDP: AS determine lowest assignment level version 0.6 ** #

# create vectors useful to determine lowest assigned level for each feature
A.int_ggs_str <- "int.ggs."
A.int_slv_str <- "int.slv."
A.val_una <- "Unassigned"
A.flt_col <- "FeatureID"

# for SILVA assignments, create vectors of unique values...
# ... using rows in all of the truncated level cols
# NOTE: these are also useful when looking for assignments that need "fixed"
A.int_slv_K_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[1]])
A.int_slv_P_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[2]])
A.int_slv_C_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[3]])
A.int_slv_O_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[4]])
A.int_slv_F_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[5]])
A.int_slv_G_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[6]])
A.int_slv_S_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[7]])

# create a vector of values we do not want in our lowest assignment column
# NOTE: uses partial string matching
A.slv_pattern <- paste(c("unass", "uncul", "unide", "metagen", "enrich"), 
                       collapse = "|")

# grep through the _unq dfs to create new vectors with the unwanted values
A.int_slv_K_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_K_unq,
                            ignore.case = T, value = T)
A.int_slv_K_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_K_unq,
                            ignore.case = T, value = T)
A.int_slv_P_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_P_unq,
                            ignore.case = T, value = T)
A.int_slv_C_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_C_unq,
                            ignore.case = T, value = T)
A.int_slv_O_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_O_unq,
                            ignore.case = T, value = T)
A.int_slv_F_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_F_unq,
                            ignore.case = T, value = T)
A.int_slv_G_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_G_unq,
                            ignore.case = T, value = T)
A.int_slv_S_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_S_unq,
                            ignore.case = T, value = T)

# combine the above vectors, pass unique() again, and collapse everything ...
# .... into a vector for use in grepl pattern matching below
A.val_slv <- paste(unique(c(A.int_slv_K_unq_unw, A.int_slv_P_unq_unw,
                            A.int_slv_C_unq_unw, A.int_slv_O_unq_unw,
                            A.int_slv_F_unq_unw, A.int_slv_G_unq_unw,
                            A.int_slv_S_unq_unw)), collapse = "|")

# now that we have our vector, process each of the truncated level cols...
# ... to find the lowest assignment level
# NOTE: any hashed lines encountered below had 0 obs. for these data... 
# ... i.e. there was no "Unassigned" for that respective level

# Unassigned Kingdom:
# A.int_ggs_U_1 <- dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col),
#                                   dplyr::all_vars(. == A.val_una))
# A.int_ggs_U_2 <- dplyr::select(A.int_ggs_U_1,
#                                dplyr::one_of(A.flt_col, A.ggs_col[1]))
# A.int_ggs_U_3 <- dplyr::rename(A.int_ggs_U_2, int.ggs.lws.txn = A.ggs_col[1])
# A.int_ggs_U <- A.int_ggs_U_3
# A.int_ggs_U$int.ggs.lws.lvl <- A.val_una

# A.int_slv_U_1 <- dplyr::filter_at(A.test_slv_trunc_fix,
#                                   dplyr::vars(A.slv_col),
#                                   dplyr::all_vars(
#                                     grepl(pattern = A.val_slv, x = .)))
# A.int_slv_U_2 <- dplyr::select(A.int_slv_U_1,
#                                dplyr::one_of(A.flt_col, A.slv_col[1]))
# A.int_slv_U_3 <- dplyr::rename(A.int_slv_U_2, int.slv.lws.txn = A.slv_col[1])
# A.int_slv_U <- A.int_slv_U_3
# A.int_slv_U$int.slv.lws.lvl <- A.val_una

# Assigned Kingdom; Unassigned Phylum:
A.int_ggs_K_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[2:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_K_2 <- dplyr::select(A.int_ggs_K_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[1]))
A.int_ggs_K_3 <- dplyr::rename(A.int_ggs_K_2, int.ggs.lws.txn = A.ggs_col[1])
A.int_ggs_K <- A.int_ggs_K_3
A.int_ggs_K$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[1], A.int_ggs_str)[[1]][2]

A.int_slv_K_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[2:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_K_2 <- dplyr::select(A.int_slv_K_1,
                               dplyr::one_of(A.flt_col, A.slv_col[1]))
A.int_slv_K_3 <- dplyr::rename(A.int_slv_K_2, int.slv.lws.txn = A.slv_col[1])
A.int_slv_K <- A.int_slv_K_3
A.int_slv_K$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[1], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum; Unassigned Class:
# A.int_ggs_P_1 <- dplyr::inner_join(
#   dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:2]),
#                    dplyr::all_vars(!. == A.val_una)),
#   dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[3:7]),
#                    dplyr::all_vars(. == A.val_una)))
# A.int_ggs_P_2 <- dplyr::select(A.int_ggs_P_1,
#                                dplyr::one_of(A.flt_col, A.ggs_col[2]))
# A.int_ggs_P_3 <- dplyr::rename(A.int_ggs_P_2, int.ggs.lws.txn = A.ggs_col[2])
# A.int_ggs_P <- A.int_ggs_P_3
# A.int_ggs_P$int.ggs.lws.lvl <- base::strsplit(
#   A.ggs_col[2], A.int_ggs_str)[[1]][2]

A.int_slv_P_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:2]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[3:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_P_2 <- dplyr::select(A.int_slv_P_1,
                               dplyr::one_of(A.flt_col, A.slv_col[2]))
A.int_slv_P_3 <- dplyr::rename(A.int_slv_P_2, int.slv.lws.txn = A.slv_col[2])
A.int_slv_P <- A.int_slv_P_3
A.int_slv_P$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[2], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class; Unassigned Order:
A.int_ggs_C_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:3]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[4:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_C_2 <- dplyr::select(A.int_ggs_C_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[3]))
A.int_ggs_C_3 <- dplyr::rename(A.int_ggs_C_2, int.ggs.lws.txn = A.ggs_col[3])
A.int_ggs_C <- A.int_ggs_C_3
A.int_ggs_C$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[3], A.int_ggs_str)[[1]][2]

A.int_slv_C_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:3]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[4:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_C_2 <- dplyr::select(A.int_slv_C_1,
                               dplyr::one_of(A.flt_col, A.slv_col[3]))
A.int_slv_C_3 <- dplyr::rename(A.int_slv_C_2, int.slv.lws.txn = A.slv_col[3])
A.int_slv_C <- A.int_slv_C_3
A.int_slv_C$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[3], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order; Unassigned Family:
A.int_ggs_O_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:4]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[5:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_O_2 <- dplyr::select(A.int_ggs_O_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[4]))
A.int_ggs_O_3 <- dplyr::rename(A.int_ggs_O_2, int.ggs.lws.txn = A.ggs_col[4])
A.int_ggs_O <- A.int_ggs_O_3
A.int_ggs_O$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[4], A.int_ggs_str)[[1]][2]

A.int_slv_O_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:4]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[5:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_O_2 <- dplyr::select(A.int_slv_O_1,
                               dplyr::one_of(A.flt_col, A.slv_col[4]))
A.int_slv_O_3 <- dplyr::rename(A.int_slv_O_2, int.slv.lws.txn = A.slv_col[4])
A.int_slv_O <- A.int_slv_O_3
A.int_slv_O$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[4], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order-Family; Unassigned Genus:
A.int_ggs_F_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:5]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[6:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_F_2 <- dplyr::select(A.int_ggs_F_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[5]))
A.int_ggs_F_3 <- dplyr::rename(A.int_ggs_F_2, int.ggs.lws.txn = A.ggs_col[5])
A.int_ggs_F <- A.int_ggs_F_3
A.int_ggs_F$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[5], A.int_ggs_str)[[1]][2]

A.int_slv_F_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:5]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[6:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_F_2 <- dplyr::select(A.int_slv_F_1,
                               dplyr::one_of(A.flt_col, A.slv_col[5]))
A.int_slv_F_3 <- dplyr::rename(A.int_slv_F_2, int.slv.lws.txn = A.slv_col[5])
A.int_slv_F <- A.int_slv_F_3
A.int_slv_F$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[5], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order-Family-Genus; Unassigned Species:
A.int_ggs_G_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:6]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_G_2 <- dplyr::select(A.int_ggs_G_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[6]))
A.int_ggs_G_3 <- dplyr::rename(A.int_ggs_G_2, int.ggs.lws.txn = A.ggs_col[6])
A.int_ggs_G <- A.int_ggs_G_3
A.int_ggs_G$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[6], A.int_ggs_str)[[1]][2]

A.int_slv_G_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:6]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_G_2 <- dplyr::select(A.int_slv_G_1,
                               dplyr::one_of(A.flt_col, A.slv_col[6]))
A.int_slv_G_3 <- dplyr::rename(A.int_slv_G_2, int.slv.lws.txn = A.slv_col[6])
A.int_slv_G <- A.int_slv_G_3
A.int_slv_G$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[6], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order-Family-Genus-Species:
A.int_ggs_S_1 <- dplyr::filter_at(A.test_ggs_trunc,
                                  dplyr::vars(A.ggs_col[1:7]),
                                  dplyr::all_vars(!. == A.val_una))
A.int_ggs_S_2 <- dplyr::select(A.int_ggs_S_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[7]))
A.int_ggs_S_3 <- dplyr::rename(A.int_ggs_S_2, int.ggs.lws.txn = A.ggs_col[7])
A.int_ggs_S <- A.int_ggs_S_3
A.int_ggs_S$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[7], A.int_ggs_str)[[1]][2]

A.int_slv_S_1 <- dplyr::filter_at(A.test_slv_trunc_fix,
                                  dplyr::vars(A.slv_col[1:7]),
                                  dplyr::all_vars(
                                    !grepl(pattern = A.val_slv, x = .)))
A.int_slv_S_2 <- dplyr::select(A.int_slv_S_1,
                               dplyr::one_of(A.flt_col, A.slv_col[7]))
A.int_slv_S_3 <- dplyr::rename(A.int_slv_S_2, int.slv.lws.txn = A.slv_col[7])
A.int_slv_S <- A.int_slv_S_3
A.int_slv_S$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[7], A.int_slv_str)[[1]][2]

# combine the above dfs
A.int_ggs_lws <- rbind(
  #A.int_ggs_U, 
  A.int_ggs_K, 
  #A.int_ggs_P,
  A.int_ggs_C, 
  A.int_ggs_O, 
  A.int_ggs_F,
  A.int_ggs_G, 
  A.int_ggs_S)

A.int_slv_lws <- rbind(
  #A.int_slv_U, 
  A.int_slv_K, 
  A.int_slv_P,
  A.int_slv_C, 
  A.int_slv_O, 
  A.int_slv_F,
  A.int_slv_G, 
  A.int_slv_S)

# run testif_fts_lost() and gather_lost_fts() functions to test for errors
# i.e. check that no features were lost in either of the lowest level dfs
A.test_ggs_lws <- testif_fts_lost(raw_df = A.int_ggs, 
                                  test_df = A.int_ggs_lws, 
                                  test_type = "lws", 
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.test_slv_lws <- testif_fts_lost(raw_df = A.int_slv, 
                                  test_df = A.int_slv_lws, 
                                  test_type = "lws",
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

A.lost_ggs_lws <- gather_lost_fts(raw_df = A.int_ggs, 
                                  test_df = A.int_ggs_lws, 
                                  test_type = "lws", 
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.lost_slv_lws <- gather_lost_fts(raw_df = A.int_slv, 
                                  test_df = A.int_slv_lws, 
                                  test_type = "lws", 
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

# if no errors were produced, merge the lws dfs with the trunc or trunc_fix dfs
A.int_ggs_frm <- merge(A.test_ggs_trunc, A.test_ggs_lws, by = "FeatureID",
                       sort = F)
A.int_slv_frm <- merge(A.int_slv_trunc_fix, A.test_slv_lws, by = "FeatureID",
                       sort = F)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost in the newly formatted dfs
A.test_ggs_frm <- testif_fts_lost(raw_df = A.int_ggs,
                                  test_df = A.int_ggs_frm,
                                  test_type = "frm",
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.test_slv_frm <- testif_fts_lost(raw_df = A.int_slv,
                                  test_df = A.int_slv_frm,
                                  test_type = "frm",
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

A.lost_ggs_frm <- gather_lost_fts(raw_df = A.int_ggs,
                                  test_df = A.int_ggs_frm,
                                  test_type = "frm",
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.lost_slv_frm <- gather_lost_fts(raw_df = A.int_slv,
                                  test_df = A.int_slv_frm,
                                  test_type = "frm",
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

# if no errors were produced, merge the formatted Greengenes and SILVA dfs
A.tax_mrg <- merge(A.test_ggs_frm, A.test_slv_frm, by = "FeatureID",
                   sort = F)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost in the merge
A.test_tax_mrg <- testif_fts_lost(raw_df = A.int_slv,
                                  test_df = A.tax_mrg,
                                  test_type = "mrg",
                                  tax_type = "merging",
                                  col_join = "FeatureID")

A.lost_tax_mrg <- gather_lost_fts(raw_df = A.int_slv,
                                  test_df = A.tax_mrg,
                                  test_type = "mrg",
                                  tax_type = "merging",
                                  col_join = "FeatureID")

### ************************************
### A - STEP  2 - create master table (EBTKS) ----
### ************************************

# ** note for KDP: AS create master table (EBTKS) version 0.7 ** #

# provide provenance for information gathering at end of STEP 2:
A.prov_secstep_AS2 <- "Section A - STEP 2"
A.prov_heading_AS2 <- "create master table (EBTKS)"
A.prov_output_obj_AS2 <- "A.EBTKS_abs_raw" # this object is output to the vault

## this step creates a master table containing:
# 'FeatureID' = MD5 sums generated in dada2 QIIME 2 feature table
# 'RepSeq' = representative DNA sequence for respective FeatureID
# 'int.ggs.tax' = taxonomic lineage with Greengenes naming convention
# 'int.ggs.cnf'= confidence in Greengenes taxonomic classification
# 'int.ggs.L' = taxon for respective level (L)
# 'int.ggs.lws.txn' = lowest assignment taxon by Greengenes
# 'int.ggs.lws.lvl' = lowest assignment level by Greengenes
# 'int.slv.tax' = taxonomic lineage with SILVA naming convention
# 'int.slv.cnf' = confidence in SILVA taxonomic classification
# 'int.slv.L' = taxon for respective level (L)
# 'int.slv.lws.txn' = lowest assignment taxon by SILVA
# 'int.slv.lws.lvl' = lowest assignment level by SILVA
# remaining columns = SampleID; row values = absolute counts of features

# read in dada2 feature table with absolute counts, and rename col #OTUID
A.raw_fts <- read.table(A.ifv_fts_tab, header = T, sep = "\t", as.is = T,
                        stringsAsFactors = F, check.names = F, skip = 1,
                        comment.char = "")

A.fts_abs <- dplyr::rename(A.raw_fts, FeatureID = "#OTU ID")

# read in the representative sequences fasta file as a data.frame (df)
A.rep_seq <- fasta_to_df(fasta_file = A.ifv_rep_seq, hdr_colname = "FeatureID",
                         seq_colname = "RepSeq")

# convert absolute feature counts to relative abundances
# NOTE: indexing with the numeric output of the which() function preserves... 
# ... column FeatureID no matter where it is located in the df
A.fts_rel <- A.fts_abs
A.fts_num <- which(names(A.fts_rel) == "FeatureID")
A.fts_rel[, -A.fts_num] <- lapply(A.fts_rel[, -A.fts_num],
                                  function(x) {x/sum(x)})

# optional, check conversion by running: colSums(A.fts_rel[, -A.fts_num])

# merge taxa, rep seqs, and fts dfs, & ensure that col FeatureID is a character
A.seq_tax <- merge(A.rep_seq, A.test_tax_mrg, by = "FeatureID", sort = F)
A.EBTKS_abs_raw <- merge(A.seq_tax, A.fts_abs, by = "FeatureID", sort = F)
A.EBTKS_rel_raw <- merge(A.seq_tax, A.fts_rel, by = "FeatureID", sort = F)

A.EBTKS_abs_raw$FeatureID <- as.character(A.EBTKS_abs_raw$FeatureID)
A.EBTKS_rel_raw$FeatureID <- as.character(A.EBTKS_rel_raw$FeatureID)

# ** note for KDP: something similar testif_fts_lost() might be needed here ** #

# gather info and provide some level of provenance outside of the script file
A.info_EBTKS_raw_AS2_a <- data.frame(
  "info" = "Number of samples",
  value = ncol(dplyr::select(A.EBTKS_abs_raw, -dplyr::one_of(com_col))),
  "object" = A.prov_output_obj_AS2, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS2, "heading" = A.prov_heading_AS2, 
  stringsAsFactors = F)
A.info_EBTKS_raw_AS2_b <- data.frame(
  "info" = "Number of features",
  value = nrow(A.EBTKS_abs_raw),
  "object" = A.prov_output_obj_AS2, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS2, "heading" = A.prov_heading_AS2, 
  stringsAsFactors = F)
A.info_EBTKS_raw_AS2_c <- data.frame(
  "info" = "Total frequency",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_raw, -dplyr::one_of(com_col)))),
  "object" = A.prov_output_obj_AS2, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS2, "heading" = A.prov_heading_AS2, 
  stringsAsFactors = F)

A.info_AS2 <- rbind(A.info_EBTKS_raw_AS2_a, A.info_EBTKS_raw_AS2_b,
                    A.info_EBTKS_raw_AS2_c)

### ************************************
### A - STEP 3a - process EBTKS: remove 'contaminants' ----
### ************************************

# ** note for KDP: AS remove 'contaminants' version 0.5 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3a <- "Section A - STEP 3a"
A.prov_heading_AS3a <- "process EBTKS: remove 'contaminants'"
A.prov_object1_AS3a <- "A.ctm_abs"
A.prov_object2_AS3a <- "A.EBTKS_abs_pro_a"

# this sub-step uses taxonomic assignments from A.tax_mrg:
# to identify and remove features likely to be contamination ...
# ... e.g. chloroplast, mitohondria, eukaryotic, Unassigned Kingdom
A.clp_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "chloropl",
                                                x = A.tax_mrg$int.ggs.tax,
                                                ignore.case = T))
A.clp_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "chloropl",
                                                x = A.tax_mrg$int.slv.tax,
                                                ignore.case = T))
A.mtc_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "mitochon",
                                                x = A.tax_mrg$int.ggs.tax,
                                                ignore.case = T))
A.mtc_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "mitochon",
                                                x = A.tax_mrg$int.slv.tax,
                                                ignore.case = T))
A.euk_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "eukaryot",
                                                x = A.tax_mrg$int.ggs.tax,
                                                ignore.case = T))
A.euk_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "eukaryot",
                                                x = A.tax_mrg$int.slv.tax,
                                                ignore.case = T))

# careful here - CANNOT be partial match
A.unk_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "Unassigned",
                                                x = A.tax_mrg$int.ggs.tax,
                                                fixed = T))
A.unk_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "Unassigned",
                                                x = A.tax_mrg$int.slv.tax,
                                                fixed = T))

# successively merge Greengenes and SILVA contaminant dfs
# unequal df lengths are handled by not specifying a 'by' variable
A.clp_int_mrg <- merge(A.clp_int_ggs, A.clp_int_slv, all = T, sort = F)
A.mtc_int_mrg <- merge(A.mtc_int_ggs, A.mtc_int_slv, all = T, sort = F)
A.euk_int_mrg <- merge(A.euk_int_ggs, A.euk_int_slv, all = T, sort = F)
A.unk_int_mrg <- merge(A.unk_int_ggs, A.unk_int_slv, all = T, sort = F)

# combine contaminant dfs
A.ctm <- rbind(A.clp_int_mrg, A.mtc_int_mrg, A.euk_int_mrg, A.unk_int_mrg)

# create ctm dfs with absolute or relative abundances
A.ctm_abs <- dplyr::semi_join(x = A.EBTKS_abs_raw, y = A.ctm, by = "FeatureID")
A.ctm_rel <- dplyr::semi_join(x = A.EBTKS_rel_raw, y = A.ctm, by = "FeatureID")

# remove contaminants from A.EBTKS_abs_raw
A.EBTKS_abs_pro_a <- dplyr::anti_join(x = A.EBTKS_abs_raw, y = A.ctm,
                                      by = "FeatureID")

### ************************************
### A - STEP 3b - process EBTKS: remove low feature count samples ----
### ************************************

# ** note for KDP: AS remove low feature count samples version 0.5 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3b <- "Section A - STEP 3b"
A.prov_heading_AS3b <- "process EBTKS: remove low feature count samples"
A.prov_object1_AS3b <- "A.low_thrsh"
A.prov_object2_AS3b <- "A.low_fts_sID"
A.prov_object3_AS3b <- "A.EBTKS_abs_pro_b"

# this sub-step identifies and removes samples with:
# low total feature count for features not considered 'contaminants'
# NOTE: this step DOES NOT remove any protocol control samples... 
# ... these are expected to have low total feature counts

# this process occurs as follows:
# (1) define the low total feature count threshold 
# (2) remove unneeded columns from A.EBTKS_abs_pro_a
# (3) sum absolute feature counts for each sample
# (4) create a vector of SampleIDs with low total feature count (to be removed)
# (5) ensure that SampleID is a character and not a factor
# (6) remove protocol controls from the sID vector
# (7) rename the sID vector
# (8) remove the identified samples from A.EBTKS_abs_pro_a

A.low_thrsh <- 999
A.EBTKS_abs_pro_frm <- dplyr::select(A.EBTKS_abs_pro_a, -dplyr::one_of(com_col))
A.EBTKS_abs_pro_fts_sum <- data.frame("SampleID" = names(A.EBTKS_abs_pro_frm),
                                      "Total" =  colSums(A.EBTKS_abs_pro_frm),
                                      row.names = NULL, stringsAsFactors = F)
A.low_fts_smp <- dplyr::filter(A.EBTKS_abs_pro_fts_sum, Total < A.low_thrsh)
A.low_fts_sID_0 <- as.character(A.low_fts_smp$SampleID)
A.low_fts_sID_1 <- A.low_fts_sID_0[!grepl(paste0("neg", collapse = "|"),
                                          A.low_fts_sID_0)]
A.low_fts_sID_2 <- A.low_fts_sID_1[!grepl(paste0("Swab", collapse = "|"),
                                          A.low_fts_sID_1)]
A.low_fts_sID_3 <- A.low_fts_sID_2[!grepl(paste0("Blank", collapse = "|"),
                                          A.low_fts_sID_2)]
A.low_fts_sID_4 <- A.low_fts_sID_3[!grepl(paste0("H2O", collapse = "|"),
                                          A.low_fts_sID_3)]
A.low_fts_sID <- A.low_fts_sID_4
A.EBTKS_abs_pro_b <- dplyr::select(A.EBTKS_abs_pro_a, 
                                   -dplyr::one_of(A.low_fts_sID))

### ************************************
### A - STEP 3c - process EBTKS: remove high 'contaminant' count samples ----
### ************************************

# ** note for KDP: AS remove high 'contaminant' count samples version 0.5 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3c <- "Section A - STEP 3c"
A.prov_heading_AS3c <- "process EBTKS: remove high 'contaminant' count samples"
A.prov_object1_AS3c <- "A.hgh_thrsh"
A.prov_object2_AS3c <- "A.hgh_ctm_sID"

# this sub-step identifies and removes samples with:
# high total relative abundance for features considered 'contaminants'
# NOTE: this can remove protocol control samples

# this process occurs as follows:
# (1) define high total relative abundance threshold 
# (2) remove unneeded columns from A.ctm_rel
# (3) sum relative abundances for each sample
# (4) create a vector of SampleIDs with high contam. rel. abund. (to be removed)
# (5) ensure that SampleID is a character and not a factor
# (6) remove the identified samples from A.EBTKS_abs_pro_b

A.hgh_thrsh <- 0.01
A.ctm_rel_frm <- dplyr::select(A.ctm_rel, -dplyr::one_of(com_col))
A.ctm_rel_frm_sum <- data.frame("SampleID" = names(A.ctm_rel_frm),
                                "Total" = colSums(A.ctm_rel_frm),
                                row.names = NULL, stringsAsFactors = F)
A.hgh_ctm_smp <- dplyr::filter(A.ctm_rel_frm_sum, Total > A.hgh_thrsh)
A.hgh_ctm_sID <- as.character(A.hgh_ctm_smp$SampleID)
A.EBTKS_abs_pro_c <- dplyr::select(A.EBTKS_abs_pro_b,
                                   -dplyr::one_of(A.hgh_ctm_sID))

### ************************************
### A - STEP 3d - process EBTKS: remove zero count features ----
### ************************************

# ** note for KDP: AS remove zero count features version 0.7 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3d <- "Section A - STEP 3d"
A.prov_heading_AS3d <- "process EBTKS: remove zero count features"
A.prov_output_obj_AS3d <- "A.EBTKS_abs_pro" # this object is output to the vault

# this sub-step identifies and removes features with:
# a total count of less than 2* across all samples^
# *threshold is in accordance with the dada2 feature table produced in QIIME 2
# ^this can & often does occur when samples are removed from a feature table ...
# ... which is a possibility after running through A - STEP 3c above
# zero count features should be removed prior to any downstream analysis

# this process occurs as follows:
# (1) create a copy of A.EBTKS_abs_pro_c to avoid overwriting data
# (2) convert column FeatureID into row.names
# (3) remove unneeded columns
# (4) sum absolute counts for each feature across all samples
# (5) create a data.frame of features that will be removed
# (6) ensure that FeatureID is a character and not a factor
# (7) anti join A.EBTKS_abs_pro_c to remove all features present in the rmv df
# NOTE: the naming convention where the sub-step's letter is appended ...
# ... to the end of the processed data.frame is dropped in this sub-step

A.EBTKS_abs_pro_1 <- A.EBTKS_abs_pro_c
row.names(A.EBTKS_abs_pro_1) <- A.EBTKS_abs_pro_1$FeatureID
A.EBTKS_abs_pro_2 <- dplyr::select(A.EBTKS_abs_pro_1, -dplyr::one_of(com_col))
A.EBTKS_abs_pro_sum <- data.frame("FeatureID" = row.names(A.EBTKS_abs_pro_2),
                                  "FeatureTotal" = rowSums(A.EBTKS_abs_pro_2),
                                  row.names = NULL, stringsAsFactors = F)
A.EBTKS_abs_pro_rmv <- dplyr::filter(A.EBTKS_abs_pro_sum, FeatureTotal < 2)
A.EBTKS_abs_pro_rmv$FeatureID <- as.character(A.EBTKS_abs_pro_rmv$FeatureID)
A.EBTKS_abs_pro <- dplyr::anti_join(x = A.EBTKS_abs_pro_c, 
                                    y = A.EBTKS_abs_pro_rmv, by = "FeatureID")

### ************************************
### A - STEP 3i - information gathering and provenance ----
### ************************************

# ** note for KDP: AS information gathering and provenance version 0.2 ** #

## A - STEP 3a:
A.info_ctm_abs_AS3a_a <- data.frame(
  "info" = "Number of features considered 'contaminants'",
  value = nrow(A.ctm_abs),
  "object" = A.prov_object1_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_ctm_abs_AS3a_b <- data.frame(
  "info" = "Total frequency of 'contaminants'",
  value = sum(colSums(dplyr::select(A.ctm_abs, -dplyr::one_of(com_col)))),
  "object" = A.prov_object1_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_ctm_abs_AS3a_c <- data.frame(
  "info" = "Number of samples 'contaminants' observed in",
  value = length(which(colSums(
    dplyr::select(A.ctm_abs, -dplyr::one_of(com_col))) > 0)),
  "object" = A.prov_object1_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_ctm_abs_AS3a_d_0 <- data.frame(
  value = which(colSums(dplyr::select(A.ctm_abs, -dplyr::one_of(com_col))) > 0))
A.info_ctm_abs_AS3a_d_1 <- A.info_ctm_abs_AS3a_d_0
A.info_ctm_abs_AS3a_d_1$info <- paste(row.names(A.info_ctm_abs_AS3a_d_1), 
                                      ": Total frequency of 'contaminants'", 
                                      sep = "")
A.info_ctm_abs_AS3a_d_2 <- A.info_ctm_abs_AS3a_d_1
A.info_ctm_abs_AS3a_d_2$object <- A.prov_object1_AS3a
A.info_ctm_abs_AS3a_d_2$script <- name_scrpt
A.info_ctm_abs_AS3a_d_2$section <- A.prov_secstep_AS3a
A.info_ctm_abs_AS3a_d_2$heading <- A.prov_heading_AS3a
A.info_ctm_abs_AS3a_d <- A.info_ctm_abs_AS3a_d_2
row.names(A.info_ctm_abs_AS3a_d) <- 1:nrow(A.info_ctm_abs_AS3a_d)

A.info_ctm_abs_AS3a <- rbind(A.info_ctm_abs_AS3a_a, A.info_ctm_abs_AS3a_b, 
                             A.info_ctm_abs_AS3a_c, A.info_ctm_abs_AS3a_d)

A.info_EBTKS_abs_pro_a_AS3a_a <- data.frame(
  "info" = "Number of samples after 'contaminant' removal",
  value = ncol(dplyr::select(A.EBTKS_abs_pro_a, -dplyr::one_of(com_col))),
  "object" = A.prov_object2_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_a_AS3a_b <- data.frame(
  "info" = "Number of features after 'contaminant' removal",
  value = nrow(A.EBTKS_abs_pro_a),
  "object" = A.prov_object2_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_a_AS3a_c <- data.frame(
  "info" = "Total frequency after 'contaminant' removal",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_pro_a, 
                                    -dplyr::one_of(com_col)))),
  "object" = A.prov_object2_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)

A.info_EBTKS_abs_pro_a_AS3a <- rbind(A.info_EBTKS_abs_pro_a_AS3a_a, 
                                     A.info_EBTKS_abs_pro_a_AS3a_b,
                                     A.info_EBTKS_abs_pro_a_AS3a_c)

A.info_AS3a <- rbind(A.info_ctm_abs_AS3a, A.info_EBTKS_abs_pro_a_AS3a)

## A - STEP 3b:
A.info_low_fts_AS3b_a <- data.frame(
  "info" = "low total feature count threshold",
  value = paste("< ", A.low_thrsh, sep = ""),
  "object" = A.prov_object1_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_low_fts_AS3b_b <- data.frame(
  "info" = "Number of samples below threshold",
  value = length(A.low_fts_sID),
  "object" = A.prov_object2_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_low_fts_AS3b_c <- data.frame(
  "info" = "SampleID below threshold",
  value = A.low_fts_sID,
  "object" = A.prov_object2_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_low_fts_smp_AS3b <- rbind(A.info_low_fts_AS3b_a, A.info_low_fts_AS3b_b,
                                 A.info_low_fts_AS3b_c)

A.info_EBTKS_abs_pro_b_AS3b_a <- data.frame(
  "info" = "Number of samples after removal of low feature count samples",
  value = ncol(dplyr::select(A.EBTKS_abs_pro_b, -dplyr::one_of(com_col))),
  "object" = A.prov_object3_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_b_AS3b_b <- data.frame(
  "info" = "Number of features after removal of low feature count samples",
  value = nrow(A.EBTKS_abs_pro_b),
  "object" = A.prov_object3_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_b_AS3b_c <- data.frame(
  "info" = "Total frequency after removal of low feature count samples",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_pro_b, 
                                    -dplyr::one_of(com_col)))),
  "object" = A.prov_object3_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)

A.info_EBTKS_abs_pro_b_AS3b <- rbind(A.info_EBTKS_abs_pro_b_AS3b_a, 
                                     A.info_EBTKS_abs_pro_b_AS3b_b,
                                     A.info_EBTKS_abs_pro_b_AS3b_c)

A.info_AS3b <- rbind(A.info_low_fts_smp_AS3b, A.info_EBTKS_abs_pro_b_AS3b)

## A - STEP 3c:
A.info_hgh_ctm_AS3c_a <- data.frame(
  "info" = "high 'contaminant' threshold",
  value = paste("> ", A.hgh_thrsh * 100, "%", sep = ""),
  "object" = A.prov_object1_AS3c, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3c, "heading" = A.prov_heading_AS3c, 
  stringsAsFactors = F)
A.info_hgh_ctm_AS3c_b <- data.frame(
  "info" = "Number of samples above threshold",
  value = length(A.hgh_ctm_sID),
  "object" = A.prov_object2_AS3c, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3c, "heading" = A.prov_heading_AS3c, 
  stringsAsFactors = F)
A.info_hgh_ctm_AS3c_c <- data.frame(
  "info" = "SampleID below threshold",
  value = A.hgh_ctm_sID,
  "object" = A.prov_object2_AS3c, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3c, "heading" = A.prov_heading_AS3c, 
  stringsAsFactors = F)
A.info_AS3c <- rbind(A.info_hgh_ctm_AS3c_a, A.info_hgh_ctm_AS3c_b,
                     A.info_hgh_ctm_AS3c_c)

# NOTE: info related to ...
# "Number of samples after removal of high 'contaminant' samples"
# "Number of features after removal of high 'contaminant' samples"
# "Total frequency after removal of high 'contaminant' samples"
# ... SHOULD NOT be calculated using the A.EBTKS_abs_pro_c df
# to understand why, read the comments at the the beginning of STEP 3d ...
# ... regarding what can happen when a sample is removed
# also, the info mentioned above is identical to the info gathered for STEP 3d

## A - STEP 3d:
A.info_EBTKS_abs_pro_AS3d_a <- data.frame(
  "info" = "Number of samples after 'process EBTKS' steps",
  value = ncol(dplyr::select(A.EBTKS_abs_pro, -dplyr::one_of(com_col))),
  "object" = A.prov_output_obj_AS3d, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3d, "heading" = A.prov_heading_AS3d, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_AS3d_b <- data.frame(
  "info" = "Number of features after 'process EBTKS' steps",
  value = nrow(A.EBTKS_abs_pro),
  "object" = A.prov_output_obj_AS3d, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3d, "heading" = A.prov_heading_AS3d, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_AS3d_c <- data.frame(
  "info" = "Total frequency after 'process EBTKS' steps",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_pro, 
                                    -dplyr::one_of(com_col)))),
  "object" = A.prov_output_obj_AS3d, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3d, "heading" = A.prov_heading_AS3d, 
  stringsAsFactors = F)

A.info_AS3d <- rbind(A.info_EBTKS_abs_pro_AS3d_a, A.info_EBTKS_abs_pro_AS3d_b,
                     A.info_EBTKS_abs_pro_AS3d_c)

# combine all A.info_ dfs created up to this point
A.info_AS3 <- rbind(A.info_AS3a, A.info_AS3b, A.info_AS3c, A.info_AS3d)

### ************************************
### A - WRITE OUTPUTS ----
### ************************************

# rbind A.info data.frames together
A.info <- rbind(A.info_AS2, A.info_AS3)

# provenance for outputs to the vault
A.prov_output1 <- data.frame("info" = "provenance for output",
                             "path" = A.ofv_EBTKS_raw,
                             "object" = A.prov_output_obj_AS2,
                             "script" = name_scrpt,
                             "section" = A.prov_secstep_AS2,
                             "heading" = A.prov_heading_AS2,
                             stringsAsFactors = F)
A.prov_output2 <- data.frame("info" = "provenance for output",
                             "path" = A.ofv_EBTKS_pro,
                             "object" = A.prov_output_obj_AS3d,
                             "script" = name_scrpt,
                             "section" = A.prov_secstep_AS3d,
                             "heading" = A.prov_heading_AS3d,
                             stringsAsFactors = F)
A.prov <- rbind(A.prov_output1, A.prov_output2)

# outputs to the vault
write.table(sep = "\t", row.names = F, 
            x = A.EBTKS_abs_raw, file = A.ofv_EBTKS_raw)
write.table(sep = "\t", row.names = F, 
            x = A.EBTKS_abs_pro, file = A.ofv_EBTKS_pro)

write.table(sep = "\t", row.names = F, x = A.info, file = A.ofv_info)
write.table(sep = "\t", row.names = F, x = A.prov, file = A.ofv_prov)

A.obj <- ls(pattern = "A.")
A.lst <- c(A.obj[grep(pattern = "A.", x = A.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS)
save(list = A.lst, file = A.ofv_wksp)

### ************************************
### S - STEP 1a - subset EBTKS: split by sample type ----
### ************************************

# ** note for KDP: code does not reflect same process used in: 
# SS subset EBTKS - version 0.1 

# this section subsets A.EBTKS_abs_pro by diet Group and Month of Study

# provide provenance for information gathering in STEP 2:
S.prov_secstep_SS1 <- "Section S - STEP 1"
S.prov_heading_SS1 <- "subset EBTKS: split by DietGroup and TimeWeeks"
S.prov_object_SecA <- "A.EBTKS_abs_pro"
S.prov_info <- paste("groups isolated from", S.prov_object_SecA, sep = " ")

# NOTE: if the environment is empty; there are some requirements:
# if section A has been run, then run the PREFACE and load section A's workspace
# if section A has not been run, then run everything above this line
# ** note for KDP:
# the choose your own adventure function that allows sections to standalone ...
# ... or stitches them together is still under construction ** #

# create new versions of the objects needed from section A ...
# ... and create a vector naming those objects (used when saving the workspace)
S.EBTKS_abs_pro <- A.EBTKS_abs_pro
S.low_fts_sID <- A.low_fts_sID
S.hgh_ctm_sID <- A.hgh_ctm_sID
S.obj_from_A <- c("A.EBTKS_abs_pro", "A.low_fts_sID", "A.hgh_ctm_sID")

# this sub-step subsets A.EBTKS_abs_pro by sample type to isolate:
# _sEPD = stool inoculum + mouse cEcum, Proximal colon, Distal colon
# _sPD = stool inoculum + mouse Proximal colon, Distal colon
# _sE = stool inoculum + mouse cEcum
# _sf = stool inoculum + mouse feces
# _s = stool inoculum
# _EPD = mouse cEcum, Proximal colon, Distal colon
# _PD = mouse Proximal colon, Distal colon
# _P = mouse Proximal colon
# _D = mouse Distal colon
# _E = mouse cEcum
# _f = mouse feces

#### ADDED fecal time series subsets
# _e1t2
# _e1t3
# _e1t4
# _e1t5
# _e1t6

# this process occurs as follows:
# (1) filter to isolate appropriate samples
# (2) check that the appropriate number of samples were retained
# (3) remove zero count feaures using same method as outlined in A - STEP 3d
S.abs_sEPD_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                              dplyr::contains("si", ignore.case = F),
                              dplyr::contains("CE", ignore.case = F),
                              dplyr::contains("CO", ignore.case = F))
S.abs_sPD_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                             dplyr::contains("si", ignore.case = F),
                             dplyr::contains("CO", ignore.case = F))
S.abs_sE_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                            dplyr::contains("si", ignore.case = F),
                            dplyr::contains("CE", ignore.case = F))
S.abs_sf_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                            dplyr::contains("si", ignore.case = F),
                            dplyr::contains("FEC", ignore.case = F))
S.abs_s_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                           dplyr::contains("si", ignore.case = F))
S.abs_EPD_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                             dplyr::contains("CE", ignore.case = F),
                             dplyr::contains("CO", ignore.case = F))
S.abs_PD_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                            dplyr::contains("CO", ignore.case = F))
S.abs_P_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                           dplyr::contains("CO_P", ignore.case = F))
S.abs_D_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                           dplyr::contains("CO_D", ignore.case = F))
S.abs_E_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                           dplyr::contains("CE", ignore.case = F))
S.abs_f_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                           dplyr::contains("FEC", ignore.case = F))

S.abs_e1t2_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                              dplyr::contains("e1t2", ignore.case = F))
S.abs_e1t3_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                              dplyr::contains("e1t3", ignore.case = F))
S.abs_e1t4_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                              dplyr::contains("e1t4", ignore.case = F))
S.abs_e1t5_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                              dplyr::contains("e1t5", ignore.case = F))
S.abs_e1t6_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                              dplyr::contains("e1t6", ignore.case = F))

S.abs_sEPD_smp_chk <- ncol(dplyr::select(S.abs_sEPD_0, -dplyr::one_of(com_col)))
S.abs_sPD_smp_chk <- ncol(dplyr::select(S.abs_sPD_0, -dplyr::one_of(com_col)))
S.abs_sE_smp_chk <- ncol(dplyr::select(S.abs_sE_0, -dplyr::one_of(com_col)))
S.abs_sf_smp_chk <- ncol(dplyr::select(S.abs_sf_0, -dplyr::one_of(com_col)))
S.abs_s_smp_chk <- ncol(dplyr::select(S.abs_s_0, -dplyr::one_of(com_col)))
S.abs_EPD_smp_chk <- ncol(dplyr::select(S.abs_EPD_0, -dplyr::one_of(com_col)))
S.abs_PD_smp_chk <- ncol(dplyr::select(S.abs_PD_0, -dplyr::one_of(com_col)))
S.abs_P_smp_chk <- ncol(dplyr::select(S.abs_P_0, -dplyr::one_of(com_col)))
S.abs_D_smp_chk <- ncol(dplyr::select(S.abs_D_0, -dplyr::one_of(com_col)))
S.abs_E_smp_chk <- ncol(dplyr::select(S.abs_E_0, -dplyr::one_of(com_col)))
S.abs_f_smp_chk <- ncol(dplyr::select(S.abs_f_0, -dplyr::one_of(com_col)))
S.abs_e1t2_smp_chk <- ncol(dplyr::select(S.abs_e1t2_0, -dplyr::one_of(com_col)))
S.abs_e1t3_smp_chk <- ncol(dplyr::select(S.abs_e1t3_0, -dplyr::one_of(com_col)))
S.abs_e1t4_smp_chk <- ncol(dplyr::select(S.abs_e1t4_0, -dplyr::one_of(com_col)))
S.abs_e1t5_smp_chk <- ncol(dplyr::select(S.abs_e1t5_0, -dplyr::one_of(com_col)))
S.abs_e1t6_smp_chk <- ncol(dplyr::select(S.abs_e1t6_0, -dplyr::one_of(com_col)))

# print(S.abs_sEPD_smp_chk) # 93
# print(S.abs_sPD_smp_chk) # 66
# print(S.abs_sE_smp_chk) # 44
# print(S.abs_sf_smp_chk) # 59
# print(S.abs_s_smp_chk) # 17
# print(S.abs_EPD_smp_chk) # 76
# print(S.abs_PD_smp_chk) # 49
# print(S.abs_P_smp_chk) # 25
# print(S.abs_D_smp_chk) # 24
# print(S.abs_E_smp_chk) # 27
# print(S.abs_f_smp_chk) # 42
# print(S.abs_e1t2_smp_chk) # 4
# print(S.abs_e1t3_smp_chk) # 8
# print(S.abs_e1t4_smp_chk) # 8
# print(S.abs_e1t5_smp_chk) # 8
# print(S.abs_e1t6_smp_chk) # 8

S.abs_sEPD_1 <- S.abs_sEPD_0
S.abs_sPD_1 <- S.abs_sPD_0
S.abs_sE_1 <- S.abs_sE_0
S.abs_sf_1 <- S.abs_sf_0
S.abs_s_1 <- S.abs_s_0
S.abs_EPD_1 <- S.abs_EPD_0
S.abs_PD_1 <- S.abs_PD_0
S.abs_P_1 <- S.abs_P_0
S.abs_D_1 <- S.abs_D_0
S.abs_E_1 <- S.abs_E_0
S.abs_f_1 <- S.abs_f_0
S.abs_e1t2_1 <- S.abs_e1t2_0
S.abs_e1t3_1 <- S.abs_e1t3_0
S.abs_e1t4_1 <- S.abs_e1t4_0
S.abs_e1t5_1 <- S.abs_e1t5_0
S.abs_e1t6_1 <- S.abs_e1t6_0

row.names(S.abs_sEPD_1) <- S.abs_sEPD_1$FeatureID
row.names(S.abs_sPD_1) <- S.abs_sPD_1$FeatureID
row.names(S.abs_sE_1) <- S.abs_sE_1$FeatureID
row.names(S.abs_sf_1) <- S.abs_sf_1$FeatureID
row.names(S.abs_s_1) <- S.abs_s_1$FeatureID
row.names(S.abs_EPD_1) <- S.abs_EPD_1$FeatureID
row.names(S.abs_PD_1) <- S.abs_PD_1$FeatureID
row.names(S.abs_P_1) <- S.abs_P_1$FeatureID
row.names(S.abs_D_1) <- S.abs_D_1$FeatureID
row.names(S.abs_E_1) <- S.abs_E_1$FeatureID
row.names(S.abs_f_1) <- S.abs_f_1$FeatureID
row.names(S.abs_e1t2_1) <- S.abs_e1t2_1$FeatureID
row.names(S.abs_e1t3_1) <- S.abs_e1t3_1$FeatureID
row.names(S.abs_e1t4_1) <- S.abs_e1t4_1$FeatureID
row.names(S.abs_e1t5_1) <- S.abs_e1t5_1$FeatureID
row.names(S.abs_e1t6_1) <- S.abs_e1t6_1$FeatureID

S.abs_sEPD_2 <- dplyr::select(S.abs_sEPD_1, -dplyr::one_of(com_col))
S.abs_sPD_2 <- dplyr::select(S.abs_sPD_1, -dplyr::one_of(com_col))
S.abs_sE_2 <- dplyr::select(S.abs_sE_1, -dplyr::one_of(com_col))
S.abs_sf_2 <- dplyr::select(S.abs_sf_1, -dplyr::one_of(com_col))
S.abs_s_2 <- dplyr::select(S.abs_s_1, -dplyr::one_of(com_col))
S.abs_EPD_2 <- dplyr::select(S.abs_EPD_1, -dplyr::one_of(com_col))
S.abs_PD_2 <- dplyr::select(S.abs_PD_1, -dplyr::one_of(com_col))
S.abs_P_2 <- dplyr::select(S.abs_P_1, -dplyr::one_of(com_col))
S.abs_D_2 <- dplyr::select(S.abs_D_1, -dplyr::one_of(com_col))
S.abs_E_2 <- dplyr::select(S.abs_E_1, -dplyr::one_of(com_col))
S.abs_f_2 <- dplyr::select(S.abs_f_1, -dplyr::one_of(com_col))
S.abs_e1t2_2 <- dplyr::select(S.abs_e1t2_1, -dplyr::one_of(com_col))
S.abs_e1t3_2 <- dplyr::select(S.abs_e1t3_1, -dplyr::one_of(com_col))
S.abs_e1t4_2 <- dplyr::select(S.abs_e1t4_1, -dplyr::one_of(com_col))
S.abs_e1t5_2 <- dplyr::select(S.abs_e1t5_1, -dplyr::one_of(com_col))
S.abs_e1t6_2 <- dplyr::select(S.abs_e1t6_1, -dplyr::one_of(com_col))

S.abs_sEPD_sum <- data.frame("FeatureID" = row.names(S.abs_sEPD_2),
                             "FeatureTotal" = rowSums(S.abs_sEPD_2),
                             row.names = NULL)
S.abs_sPD_sum <- data.frame("FeatureID" = row.names(S.abs_sPD_2),
                            "FeatureTotal" = rowSums(S.abs_sPD_2),
                            row.names = NULL)
S.abs_sE_sum <- data.frame("FeatureID" = row.names(S.abs_sE_2),
                           "FeatureTotal" = rowSums(S.abs_sE_2),
                           row.names = NULL)
S.abs_sf_sum <- data.frame("FeatureID" = row.names(S.abs_sf_2),
                           "FeatureTotal" = rowSums(S.abs_sf_2),
                           row.names = NULL)
S.abs_s_sum <- data.frame("FeatureID" = row.names(S.abs_s_2),
                          "FeatureTotal" = rowSums(S.abs_s_2),
                          row.names = NULL)
S.abs_EPD_sum <- data.frame("FeatureID" = row.names(S.abs_EPD_2),
                            "FeatureTotal" = rowSums(S.abs_EPD_2),
                            row.names = NULL)
S.abs_PD_sum <- data.frame("FeatureID" = row.names(S.abs_PD_2),
                           "FeatureTotal" = rowSums(S.abs_PD_2),
                           row.names = NULL)
S.abs_P_sum <- data.frame("FeatureID" = row.names(S.abs_P_2),
                          "FeatureTotal" = rowSums(S.abs_P_2),
                          row.names = NULL)
S.abs_D_sum <- data.frame("FeatureID" = row.names(S.abs_D_2),
                          "FeatureTotal" = rowSums(S.abs_D_2),
                          row.names = NULL)
S.abs_E_sum <- data.frame("FeatureID" = row.names(S.abs_E_2),
                          "FeatureTotal" = rowSums(S.abs_E_2),
                          row.names = NULL)
S.abs_f_sum <- data.frame("FeatureID" = row.names(S.abs_f_2),
                          "FeatureTotal" = rowSums(S.abs_f_2),
                          row.names = NULL)
S.abs_e1t2_sum <- data.frame("FeatureID" = row.names(S.abs_e1t2_2),
                             "FeatureTotal" = rowSums(S.abs_e1t2_2),
                             row.names = NULL)
S.abs_e1t3_sum <- data.frame("FeatureID" = row.names(S.abs_e1t3_2),
                             "FeatureTotal" = rowSums(S.abs_e1t3_2),
                             row.names = NULL)
S.abs_e1t4_sum <- data.frame("FeatureID" = row.names(S.abs_e1t4_2),
                             "FeatureTotal" = rowSums(S.abs_e1t4_2),
                             row.names = NULL)
S.abs_e1t5_sum <- data.frame("FeatureID" = row.names(S.abs_e1t5_2),
                             "FeatureTotal" = rowSums(S.abs_e1t5_2),
                             row.names = NULL)
S.abs_e1t6_sum <- data.frame("FeatureID" = row.names(S.abs_e1t6_2),
                             "FeatureTotal" = rowSums(S.abs_e1t6_2),
                             row.names = NULL)

S.abs_sEPD_rmv <- dplyr::filter(S.abs_sEPD_sum, FeatureTotal < 2)
S.abs_sPD_rmv <- dplyr::filter(S.abs_sPD_sum, FeatureTotal < 2)
S.abs_sE_rmv <- dplyr::filter(S.abs_sE_sum, FeatureTotal < 2)
S.abs_sf_rmv <- dplyr::filter(S.abs_sf_sum, FeatureTotal < 2)
S.abs_s_rmv <- dplyr::filter(S.abs_s_sum, FeatureTotal < 2)
S.abs_EPD_rmv <- dplyr::filter(S.abs_EPD_sum, FeatureTotal < 2)
S.abs_PD_rmv <- dplyr::filter(S.abs_PD_sum, FeatureTotal < 2)
S.abs_P_rmv <- dplyr::filter(S.abs_P_sum, FeatureTotal < 2)
S.abs_D_rmv <- dplyr::filter(S.abs_D_sum, FeatureTotal < 2)
S.abs_E_rmv <- dplyr::filter(S.abs_E_sum, FeatureTotal < 2)
S.abs_f_rmv <- dplyr::filter(S.abs_f_sum, FeatureTotal < 2)
S.abs_e1t2_rmv <- dplyr::filter(S.abs_e1t2_sum, FeatureTotal < 2)
S.abs_e1t3_rmv <- dplyr::filter(S.abs_e1t3_sum, FeatureTotal < 2)
S.abs_e1t4_rmv <- dplyr::filter(S.abs_e1t4_sum, FeatureTotal < 2)
S.abs_e1t5_rmv <- dplyr::filter(S.abs_e1t5_sum, FeatureTotal < 2)
S.abs_e1t6_rmv <- dplyr::filter(S.abs_e1t6_sum, FeatureTotal < 2)

S.abs_sEPD_rmv$FeatureID <- as.character(S.abs_sEPD_rmv$FeatureID)
S.abs_sPD_rmv$FeatureID <- as.character(S.abs_sPD_rmv$FeatureID)
S.abs_sE_rmv$FeatureID <- as.character(S.abs_sE_rmv$FeatureID)
S.abs_sf_rmv$FeatureID <- as.character(S.abs_sf_rmv$FeatureID)
S.abs_s_rmv$FeatureID <- as.character(S.abs_s_rmv$FeatureID)
S.abs_EPD_rmv$FeatureID <- as.character(S.abs_EPD_rmv$FeatureID)
S.abs_PD_rmv$FeatureID <- as.character(S.abs_PD_rmv$FeatureID)
S.abs_P_rmv$FeatureID <- as.character(S.abs_P_rmv$FeatureID)
S.abs_D_rmv$FeatureID <- as.character(S.abs_D_rmv$FeatureID)
S.abs_E_rmv$FeatureID <- as.character(S.abs_E_rmv$FeatureID)
S.abs_f_rmv$FeatureID <- as.character(S.abs_f_rmv$FeatureID)
S.abs_e1t2_rmv$FeatureID <- as.character(S.abs_e1t2_rmv$FeatureID)
S.abs_e1t3_rmv$FeatureID <- as.character(S.abs_e1t3_rmv$FeatureID)
S.abs_e1t4_rmv$FeatureID <- as.character(S.abs_e1t4_rmv$FeatureID)
S.abs_e1t5_rmv$FeatureID <- as.character(S.abs_e1t5_rmv$FeatureID)
S.abs_e1t6_rmv$FeatureID <- as.character(S.abs_e1t6_rmv$FeatureID)

S.abs_sEPD <- dplyr::anti_join(x = S.abs_sEPD_0, y = S.abs_sEPD_rmv, 
                               by = "FeatureID")
S.abs_sPD <- dplyr::anti_join(x = S.abs_sPD_0, y = S.abs_sPD_rmv, 
                              by = "FeatureID")
S.abs_sE <- dplyr::anti_join(x = S.abs_sE_0, y = S.abs_sE_rmv, 
                             by = "FeatureID")
S.abs_sf <- dplyr::anti_join(x = S.abs_sf_0, y = S.abs_sf_rmv, 
                             by = "FeatureID")
S.abs_s <- dplyr::anti_join(x = S.abs_s_0, y = S.abs_s_rmv, 
                            by = "FeatureID")
S.abs_EPD <- dplyr::anti_join(x = S.abs_EPD_0, y = S.abs_EPD_rmv, 
                              by = "FeatureID")
S.abs_PD <- dplyr::anti_join(x = S.abs_PD_0, y = S.abs_PD_rmv, 
                             by = "FeatureID")
S.abs_P <- dplyr::anti_join(x = S.abs_P_0, y = S.abs_PD_rmv, 
                            by = "FeatureID")
S.abs_D <- dplyr::anti_join(x = S.abs_D_0, y = S.abs_PD_rmv, 
                            by = "FeatureID")
S.abs_E <- dplyr::anti_join(x = S.abs_E_0, y = S.abs_E_rmv, 
                            by = "FeatureID")
S.abs_f <- dplyr::anti_join(x = S.abs_f_0, y = S.abs_f_rmv, 
                            by = "FeatureID")
S.abs_e1t2 <- dplyr::anti_join(x = S.abs_e1t2_0, y = S.abs_e1t2_rmv, 
                               by = "FeatureID")
S.abs_e1t3 <- dplyr::anti_join(x = S.abs_e1t3_0, y = S.abs_e1t3_rmv, 
                               by = "FeatureID")
S.abs_e1t4 <- dplyr::anti_join(x = S.abs_e1t4_0, y = S.abs_e1t4_rmv, 
                               by = "FeatureID")
S.abs_e1t5 <- dplyr::anti_join(x = S.abs_e1t5_0, y = S.abs_e1t5_rmv, 
                               by = "FeatureID")
S.abs_e1t6 <- dplyr::anti_join(x = S.abs_e1t6_0, y = S.abs_e1t6_rmv, 
                               by = "FeatureID")

### ************************************
### S - STEP 1b - subset EBTKS: split by sex of human donor ----
### ************************************

#### ADDED fecal time series subsets for FEMALE only
# _e1t2
# _e1t3
# _e1t4
# _e1t5
# _e1t6

# this sub-step subsets the sample type dfs by consortium group to isolate:
# _F01 = control and rice bran modified consortia from female
# _M00 = control and rice bran modified consortia from males

# this process occurs as follows:
# (1) filter to isolate appropriate samples
# (2) check that the appropriate number of samples were retained
# (3) remove zero count feaures using same method as outlined in A - STEP 3d

S.abs_sEPD_F01_0 <- dplyr::select(S.abs_sEPD, dplyr::one_of(com_col),
                                  dplyr::contains("F0", ignore.case = F))
S.abs_sPD_F01_0 <- dplyr::select(S.abs_sPD, dplyr::one_of(com_col),
                                 dplyr::contains("F0", ignore.case = F))
S.abs_sE_F01_0 <- dplyr::select(S.abs_sE, dplyr::one_of(com_col),
                                dplyr::contains("F0", ignore.case = F))
S.abs_sf_F01_0 <- dplyr::select(S.abs_sf, dplyr::one_of(com_col),
                                dplyr::contains("F0", ignore.case = F))
S.abs_s_F01_0 <- dplyr::select(S.abs_s, dplyr::one_of(com_col),
                               dplyr::contains("F0", ignore.case = F))
S.abs_EPD_F01_0 <- dplyr::select(S.abs_EPD, dplyr::one_of(com_col),
                                 dplyr::contains("F0", ignore.case = F))
S.abs_PD_F01_0 <- dplyr::select(S.abs_PD, dplyr::one_of(com_col),
                                dplyr::contains("F0", ignore.case = F))
S.abs_P_F01_0 <- dplyr::select(S.abs_P, dplyr::one_of(com_col),
                               dplyr::contains("F0", ignore.case = F))
S.abs_D_F01_0 <- dplyr::select(S.abs_D, dplyr::one_of(com_col),
                               dplyr::contains("F0", ignore.case = F))
S.abs_E_F01_0 <- dplyr::select(S.abs_E, dplyr::one_of(com_col),
                               dplyr::contains("F0", ignore.case = F))
S.abs_f_F01_0 <- dplyr::select(S.abs_f, dplyr::one_of(com_col),
                               dplyr::contains("F0", ignore.case = F))
S.abs_sEPD_M00_0 <- dplyr::select(S.abs_sEPD, dplyr::one_of(com_col),
                                  dplyr::contains("M0", ignore.case = F))
S.abs_sPD_M00_0 <- dplyr::select(S.abs_sPD, dplyr::one_of(com_col),
                                 dplyr::contains("M0", ignore.case = F))
S.abs_sE_M00_0 <- dplyr::select(S.abs_sE, dplyr::one_of(com_col),
                                dplyr::contains("M0", ignore.case = F))
S.abs_sf_M00_0 <- dplyr::select(S.abs_sf, dplyr::one_of(com_col),
                                dplyr::contains("M0", ignore.case = F))
S.abs_s_M00_0 <- dplyr::select(S.abs_s, dplyr::one_of(com_col),
                               dplyr::contains("M0", ignore.case = F))
S.abs_EPD_M00_0 <- dplyr::select(S.abs_EPD, dplyr::one_of(com_col),
                                 dplyr::contains("M0", ignore.case = F))
S.abs_PD_M00_0 <- dplyr::select(S.abs_PD, dplyr::one_of(com_col),
                                dplyr::contains("M0", ignore.case = F))
S.abs_P_M00_0 <- dplyr::select(S.abs_P, dplyr::one_of(com_col),
                               dplyr::contains("M0", ignore.case = F))
S.abs_D_M00_0 <- dplyr::select(S.abs_D, dplyr::one_of(com_col),
                               dplyr::contains("M0", ignore.case = F))
S.abs_E_M00_0 <- dplyr::select(S.abs_E, dplyr::one_of(com_col),
                               dplyr::contains("M0", ignore.case = F))
S.abs_f_M00_0 <- dplyr::select(S.abs_f, dplyr::one_of(com_col),
                               dplyr::contains("M0", ignore.case = F))

S.abs_e1t2_F01_0 <- dplyr::select(S.abs_e1t2, dplyr::one_of(com_col),
                                  dplyr::contains("F0", ignore.case = F))
S.abs_e1t3_F01_0 <- dplyr::select(S.abs_e1t3, dplyr::one_of(com_col),
                                  dplyr::contains("F0", ignore.case = F))
S.abs_e1t4_F01_0 <- dplyr::select(S.abs_e1t4, dplyr::one_of(com_col),
                                  dplyr::contains("F0", ignore.case = F))
S.abs_e1t5_F01_0 <- dplyr::select(S.abs_e1t5, dplyr::one_of(com_col),
                                  dplyr::contains("F0", ignore.case = F))
S.abs_e1t6_F01_0 <- dplyr::select(S.abs_e1t6, dplyr::one_of(com_col),
                                  dplyr::contains("F0", ignore.case = F))

S.abs_sEPD_F01_smp_chk <- ncol(dplyr::select(S.abs_sEPD_F01_0, 
                                             -dplyr::one_of(com_col)))
S.abs_sPD_F01_smp_chk <- ncol(dplyr::select(S.abs_sPD_F01_0, 
                                            -dplyr::one_of(com_col)))
S.abs_sE_F01_smp_chk <- ncol(dplyr::select(S.abs_sE_F01_0, 
                                           -dplyr::one_of(com_col)))
S.abs_sf_F01_smp_chk <- ncol(dplyr::select(S.abs_sf_F01_0, 
                                           -dplyr::one_of(com_col)))
S.abs_s_F01_smp_chk <- ncol(dplyr::select(S.abs_s_F01_0, 
                                          -dplyr::one_of(com_col)))
S.abs_EPD_F01_smp_chk <- ncol(dplyr::select(S.abs_EPD_F01_0,
                                            -dplyr::one_of(com_col)))
S.abs_PD_F01_smp_chk <- ncol(dplyr::select(S.abs_PD_F01_0,
                                           -dplyr::one_of(com_col)))
S.abs_P_F01_smp_chk <- ncol(dplyr::select(S.abs_P_F01_0,
                                          -dplyr::one_of(com_col)))
S.abs_D_F01_smp_chk <- ncol(dplyr::select(S.abs_D_F01_0,
                                          -dplyr::one_of(com_col)))
S.abs_E_F01_smp_chk <- ncol(dplyr::select(S.abs_E_F01_0, 
                                          -dplyr::one_of(com_col)))
S.abs_f_F01_smp_chk <- ncol(dplyr::select(S.abs_f_F01_0, 
                                          -dplyr::one_of(com_col)))
S.abs_sEPD_M00_smp_chk <- ncol(dplyr::select(S.abs_sEPD_M00_0, 
                                             -dplyr::one_of(com_col)))
S.abs_sPD_M00_smp_chk <- ncol(dplyr::select(S.abs_sPD_M00_0, 
                                            -dplyr::one_of(com_col)))
S.abs_sE_M00_smp_chk <- ncol(dplyr::select(S.abs_sE_M00_0, 
                                           -dplyr::one_of(com_col)))
S.abs_sf_M00_smp_chk <- ncol(dplyr::select(S.abs_sf_M00_0, 
                                           -dplyr::one_of(com_col)))
S.abs_s_M00_smp_chk <- ncol(dplyr::select(S.abs_s_M00_0, 
                                          -dplyr::one_of(com_col)))
S.abs_EPD_M00_smp_chk <- ncol(dplyr::select(S.abs_EPD_M00_0,
                                            -dplyr::one_of(com_col)))
S.abs_PD_M00_smp_chk <- ncol(dplyr::select(S.abs_PD_M00_0,
                                           -dplyr::one_of(com_col)))
S.abs_P_M00_smp_chk <- ncol(dplyr::select(S.abs_P_M00_0,
                                          -dplyr::one_of(com_col)))
S.abs_D_M00_smp_chk <- ncol(dplyr::select(S.abs_D_M00_0,
                                          -dplyr::one_of(com_col)))
S.abs_E_M00_smp_chk <- ncol(dplyr::select(S.abs_E_M00_0, 
                                          -dplyr::one_of(com_col)))
S.abs_f_M00_smp_chk <- ncol(dplyr::select(S.abs_f_M00_0, 
                                          -dplyr::one_of(com_col)))

S.abs_e1t2_F01_smp_chk <- ncol(dplyr::select(S.abs_e1t2_F01_0, 
                                             -dplyr::one_of(com_col)))
S.abs_e1t3_F01_smp_chk <- ncol(dplyr::select(S.abs_e1t3_F01_0, 
                                             -dplyr::one_of(com_col)))
S.abs_e1t4_F01_smp_chk <- ncol(dplyr::select(S.abs_e1t4_F01_0, 
                                             -dplyr::one_of(com_col)))
S.abs_e1t5_F01_smp_chk <- ncol(dplyr::select(S.abs_e1t5_F01_0, 
                                             -dplyr::one_of(com_col)))
S.abs_e1t6_F01_smp_chk <- ncol(dplyr::select(S.abs_e1t6_F01_0, 
                                             -dplyr::one_of(com_col)))

# print(S.abs_sEPD_F01_smp_chk) # 49
# print(S.abs_sPD_F01_smp_chk) # 34
# print(S.abs_sE_F01_smp_chk) # 23
# print(S.abs_sf_F01_smp_chk) # 31
# print(S.abs_s_F01_smp_chk) # 08
# print(S.abs_EPD_F01_smp_chk) # 41
# print(S.abs_PD_F01_smp_chk) # 26
# print(S.abs_P_F01_smp_chk) # 14
# print(S.abs_D_F01_smp_chk) # 12
# print(S.abs_E_F01_smp_chk) # 15
# print(S.abs_f_F01_smp_chk) # 23
# print(S.abs_sEPD_M00_smp_chk) # 44
# print(S.abs_sPD_M00_smp_chk) # 32
# print(S.abs_sE_M00_smp_chk) # 21
# print(S.abs_sf_M00_smp_chk) # 28
# print(S.abs_s_M00_smp_chk) # 09
# print(S.abs_EPD_M00_smp_chk) # 35
# print(S.abs_PD_M00_smp_chk) # 23
# print(S.abs_P_M00_smp_chk) # 11
# print(S.abs_D_M00_smp_chk) # 12
# print(S.abs_E_M00_smp_chk) # 12
# print(S.abs_f_M00_smp_chk) # 19

# print(S.abs_e1t2_F01_smp_chk) # 3
# print(S.abs_e1t3_F01_smp_chk) # 5
# print(S.abs_e1t4_F01_smp_chk) # 5
# print(S.abs_e1t5_F01_smp_chk) # 5
# print(S.abs_e1t6_F01_smp_chk) # 5

S.abs_sEPD_F01_1 <- S.abs_sEPD_F01_0
S.abs_sPD_F01_1 <- S.abs_sPD_F01_0
S.abs_sE_F01_1 <- S.abs_sE_F01_0
S.abs_sf_F01_1 <- S.abs_sf_F01_0
S.abs_s_F01_1 <- S.abs_s_F01_0
S.abs_EPD_F01_1 <- S.abs_EPD_F01_0
S.abs_PD_F01_1 <- S.abs_PD_F01_0
S.abs_P_F01_1 <- S.abs_P_F01_0
S.abs_D_F01_1 <- S.abs_D_F01_0
S.abs_E_F01_1 <- S.abs_E_F01_0
S.abs_f_F01_1 <- S.abs_f_F01_0
S.abs_sEPD_M00_1 <- S.abs_sEPD_M00_0
S.abs_sPD_M00_1 <- S.abs_sPD_M00_0
S.abs_sE_M00_1 <- S.abs_sE_M00_0
S.abs_sf_M00_1 <- S.abs_sf_M00_0
S.abs_s_M00_1 <- S.abs_s_M00_0
S.abs_EPD_M00_1 <- S.abs_EPD_M00_0
S.abs_PD_M00_1 <- S.abs_PD_M00_0
S.abs_P_M00_1 <- S.abs_P_M00_0
S.abs_D_M00_1 <- S.abs_D_M00_0
S.abs_E_M00_1 <- S.abs_E_M00_0
S.abs_f_M00_1 <- S.abs_f_M00_0

S.abs_e1t2_F01_1 <- S.abs_e1t2_F01_0
S.abs_e1t3_F01_1 <- S.abs_e1t3_F01_0
S.abs_e1t4_F01_1 <- S.abs_e1t4_F01_0
S.abs_e1t5_F01_1 <- S.abs_e1t5_F01_0
S.abs_e1t6_F01_1 <- S.abs_e1t6_F01_0

row.names(S.abs_sEPD_F01_1) <- S.abs_sEPD_F01_1$FeatureID
row.names(S.abs_sPD_F01_1) <- S.abs_sPD_F01_1$FeatureID
row.names(S.abs_sE_F01_1) <- S.abs_sE_F01_1$FeatureID
row.names(S.abs_sf_F01_1) <- S.abs_sf_F01_1$FeatureID
row.names(S.abs_s_F01_1) <- S.abs_s_F01_1$FeatureID
row.names(S.abs_EPD_F01_1) <- S.abs_EPD_F01_1$FeatureID
row.names(S.abs_PD_F01_1) <- S.abs_PD_F01_1$FeatureID
row.names(S.abs_P_F01_1) <- S.abs_P_F01_1$FeatureID
row.names(S.abs_D_F01_1) <- S.abs_D_F01_1$FeatureID
row.names(S.abs_E_F01_1) <- S.abs_E_F01_1$FeatureID
row.names(S.abs_f_F01_1) <- S.abs_f_F01_1$FeatureID
row.names(S.abs_sEPD_M00_1) <- S.abs_sEPD_M00_1$FeatureID
row.names(S.abs_sPD_M00_1) <- S.abs_sPD_M00_1$FeatureID
row.names(S.abs_sE_M00_1) <- S.abs_sE_M00_1$FeatureID
row.names(S.abs_sf_M00_1) <- S.abs_sf_M00_1$FeatureID
row.names(S.abs_s_M00_1) <- S.abs_s_M00_1$FeatureID
row.names(S.abs_EPD_M00_1) <- S.abs_EPD_M00_1$FeatureID
row.names(S.abs_PD_M00_1) <- S.abs_PD_M00_1$FeatureID
row.names(S.abs_P_M00_1) <- S.abs_P_M00_1$FeatureID
row.names(S.abs_D_M00_1) <- S.abs_D_M00_1$FeatureID
row.names(S.abs_E_M00_1) <- S.abs_E_M00_1$FeatureID
row.names(S.abs_f_M00_1) <- S.abs_f_M00_1$FeatureID

row.names(S.abs_e1t2_F01_1) <- S.abs_e1t2_F01_1$FeatureID
row.names(S.abs_e1t3_F01_1) <- S.abs_e1t3_F01_1$FeatureID
row.names(S.abs_e1t4_F01_1) <- S.abs_e1t4_F01_1$FeatureID
row.names(S.abs_e1t5_F01_1) <- S.abs_e1t5_F01_1$FeatureID
row.names(S.abs_e1t6_F01_1) <- S.abs_e1t6_F01_1$FeatureID

S.abs_sEPD_F01_2 <- dplyr::select(S.abs_sEPD_F01_1, -dplyr::one_of(com_col))
S.abs_sPD_F01_2 <- dplyr::select(S.abs_sPD_F01_1, -dplyr::one_of(com_col))
S.abs_sE_F01_2 <- dplyr::select(S.abs_sE_F01_1, -dplyr::one_of(com_col))
S.abs_sf_F01_2 <- dplyr::select(S.abs_sf_F01_1, -dplyr::one_of(com_col))
S.abs_s_F01_2 <- dplyr::select(S.abs_s_F01_1, -dplyr::one_of(com_col))
S.abs_EPD_F01_2 <- dplyr::select(S.abs_EPD_F01_1, -dplyr::one_of(com_col))
S.abs_PD_F01_2 <- dplyr::select(S.abs_PD_F01_1, -dplyr::one_of(com_col))
S.abs_P_F01_2 <- dplyr::select(S.abs_P_F01_1, -dplyr::one_of(com_col))
S.abs_D_F01_2 <- dplyr::select(S.abs_D_F01_1, -dplyr::one_of(com_col))
S.abs_E_F01_2 <- dplyr::select(S.abs_E_F01_1, -dplyr::one_of(com_col))
S.abs_f_F01_2 <- dplyr::select(S.abs_f_F01_1, -dplyr::one_of(com_col))
S.abs_sEPD_M00_2 <- dplyr::select(S.abs_sEPD_M00_1, -dplyr::one_of(com_col))
S.abs_sPD_M00_2 <- dplyr::select(S.abs_sPD_M00_1, -dplyr::one_of(com_col))
S.abs_sE_M00_2 <- dplyr::select(S.abs_sE_M00_1, -dplyr::one_of(com_col))
S.abs_sf_M00_2 <- dplyr::select(S.abs_sf_M00_1, -dplyr::one_of(com_col))
S.abs_s_M00_2 <- dplyr::select(S.abs_s_M00_1, -dplyr::one_of(com_col))
S.abs_EPD_M00_2 <- dplyr::select(S.abs_EPD_M00_1, -dplyr::one_of(com_col))
S.abs_PD_M00_2 <- dplyr::select(S.abs_PD_M00_1, -dplyr::one_of(com_col))
S.abs_P_M00_2 <- dplyr::select(S.abs_P_M00_1, -dplyr::one_of(com_col))
S.abs_D_M00_2 <- dplyr::select(S.abs_D_M00_1, -dplyr::one_of(com_col))
S.abs_E_M00_2 <- dplyr::select(S.abs_E_M00_1, -dplyr::one_of(com_col))
S.abs_f_M00_2 <- dplyr::select(S.abs_f_M00_1, -dplyr::one_of(com_col))

S.abs_e1t2_F01_2 <- dplyr::select(S.abs_e1t2_F01_1, -dplyr::one_of(com_col))
S.abs_e1t3_F01_2 <- dplyr::select(S.abs_e1t3_F01_1, -dplyr::one_of(com_col))
S.abs_e1t4_F01_2 <- dplyr::select(S.abs_e1t4_F01_1, -dplyr::one_of(com_col))
S.abs_e1t5_F01_2 <- dplyr::select(S.abs_e1t5_F01_1, -dplyr::one_of(com_col))
S.abs_e1t6_F01_2 <- dplyr::select(S.abs_e1t6_F01_1, -dplyr::one_of(com_col))

S.abs_sEPD_F01_sum <- data.frame("FeatureID" = row.names(S.abs_sEPD_F01_2),
                                 "FeatureTotal" = rowSums(S.abs_sEPD_F01_2),
                                 row.names = NULL)
S.abs_sPD_F01_sum <- data.frame("FeatureID" = row.names(S.abs_sPD_F01_2),
                                "FeatureTotal" = rowSums(S.abs_sPD_F01_2),
                                row.names = NULL)
S.abs_sE_F01_sum <- data.frame("FeatureID" = row.names(S.abs_sE_F01_2),
                               "FeatureTotal" = rowSums(S.abs_sE_F01_2),
                               row.names = NULL)
S.abs_sf_F01_sum <- data.frame("FeatureID" = row.names(S.abs_sf_F01_2),
                               "FeatureTotal" = rowSums(S.abs_sf_F01_2),
                               row.names = NULL)
S.abs_s_F01_sum <- data.frame("FeatureID" = row.names(S.abs_s_F01_2),
                              "FeatureTotal" = rowSums(S.abs_s_F01_2),
                              row.names = NULL)
S.abs_EPD_F01_sum <- data.frame("FeatureID" = row.names(S.abs_EPD_F01_2),
                                "FeatureTotal" = rowSums(S.abs_EPD_F01_2),
                                row.names = NULL)
S.abs_PD_F01_sum <- data.frame("FeatureID" = row.names(S.abs_PD_F01_2),
                               "FeatureTotal" = rowSums(S.abs_PD_F01_2),
                               row.names = NULL)
S.abs_P_F01_sum <- data.frame("FeatureID" = row.names(S.abs_P_F01_2),
                              "FeatureTotal" = rowSums(S.abs_P_F01_2),
                              row.names = NULL)
S.abs_D_F01_sum <- data.frame("FeatureID" = row.names(S.abs_D_F01_2),
                              "FeatureTotal" = rowSums(S.abs_D_F01_2),
                              row.names = NULL)
S.abs_E_F01_sum <- data.frame("FeatureID" = row.names(S.abs_E_F01_2),
                              "FeatureTotal" = rowSums(S.abs_E_F01_2),
                              row.names = NULL)
S.abs_f_F01_sum <- data.frame("FeatureID" = row.names(S.abs_f_F01_2),
                              "FeatureTotal" = rowSums(S.abs_f_F01_2),
                              row.names = NULL)
S.abs_sEPD_M00_sum <- data.frame("FeatureID" = row.names(S.abs_sEPD_M00_2),
                                 "FeatureTotal" = rowSums(S.abs_sEPD_M00_2),
                                 row.names = NULL)
S.abs_sPD_M00_sum <- data.frame("FeatureID" = row.names(S.abs_sPD_M00_2),
                                "FeatureTotal" = rowSums(S.abs_sPD_M00_2),
                                row.names = NULL)
S.abs_sE_M00_sum <- data.frame("FeatureID" = row.names(S.abs_sE_M00_2),
                               "FeatureTotal" = rowSums(S.abs_sE_M00_2),
                               row.names = NULL)
S.abs_sf_M00_sum <- data.frame("FeatureID" = row.names(S.abs_sf_M00_2),
                               "FeatureTotal" = rowSums(S.abs_sf_M00_2),
                               row.names = NULL)
S.abs_s_M00_sum <- data.frame("FeatureID" = row.names(S.abs_s_M00_2),
                              "FeatureTotal" = rowSums(S.abs_s_M00_2),
                              row.names = NULL)
S.abs_EPD_M00_sum <- data.frame("FeatureID" = row.names(S.abs_EPD_M00_2),
                                "FeatureTotal" = rowSums(S.abs_EPD_M00_2),
                                row.names = NULL)
S.abs_PD_M00_sum <- data.frame("FeatureID" = row.names(S.abs_PD_M00_2),
                               "FeatureTotal" = rowSums(S.abs_PD_M00_2),
                               row.names = NULL)
S.abs_P_M00_sum <- data.frame("FeatureID" = row.names(S.abs_P_M00_2),
                              "FeatureTotal" = rowSums(S.abs_P_M00_2),
                              row.names = NULL)
S.abs_D_M00_sum <- data.frame("FeatureID" = row.names(S.abs_D_M00_2),
                              "FeatureTotal" = rowSums(S.abs_D_M00_2),
                              row.names = NULL)
S.abs_E_M00_sum <- data.frame("FeatureID" = row.names(S.abs_E_M00_2),
                              "FeatureTotal" = rowSums(S.abs_E_M00_2),
                              row.names = NULL)
S.abs_f_M00_sum <- data.frame("FeatureID" = row.names(S.abs_f_M00_2),
                              "FeatureTotal" = rowSums(S.abs_f_M00_2),
                              row.names = NULL)

S.abs_e1t2_F01_sum <- data.frame("FeatureID" = row.names(S.abs_e1t2_F01_2),
                                 "FeatureTotal" = rowSums(S.abs_e1t2_F01_2),
                                 row.names = NULL)
S.abs_e1t3_F01_sum <- data.frame("FeatureID" = row.names(S.abs_e1t3_F01_2),
                                 "FeatureTotal" = rowSums(S.abs_e1t3_F01_2),
                                 row.names = NULL)
S.abs_e1t4_F01_sum <- data.frame("FeatureID" = row.names(S.abs_e1t4_F01_2),
                                 "FeatureTotal" = rowSums(S.abs_e1t4_F01_2),
                                 row.names = NULL)
S.abs_e1t5_F01_sum <- data.frame("FeatureID" = row.names(S.abs_e1t5_F01_2),
                                 "FeatureTotal" = rowSums(S.abs_e1t5_F01_2),
                                 row.names = NULL)
S.abs_e1t6_F01_sum <- data.frame("FeatureID" = row.names(S.abs_e1t6_F01_2),
                                 "FeatureTotal" = rowSums(S.abs_e1t6_F01_2),
                                 row.names = NULL)

S.abs_sEPD_F01_rmv <- dplyr::filter(S.abs_sEPD_F01_sum, FeatureTotal < 2)
S.abs_sPD_F01_rmv <- dplyr::filter(S.abs_sPD_F01_sum, FeatureTotal < 2)
S.abs_sE_F01_rmv <- dplyr::filter(S.abs_sE_F01_sum, FeatureTotal < 2)
S.abs_sf_F01_rmv <- dplyr::filter(S.abs_sf_F01_sum, FeatureTotal < 2)
S.abs_s_F01_rmv <- dplyr::filter(S.abs_s_F01_sum, FeatureTotal < 2)
S.abs_EPD_F01_rmv <- dplyr::filter(S.abs_EPD_F01_sum, FeatureTotal < 2)
S.abs_PD_F01_rmv <- dplyr::filter(S.abs_PD_F01_sum, FeatureTotal < 2)
S.abs_P_F01_rmv <- dplyr::filter(S.abs_P_F01_sum, FeatureTotal < 2)
S.abs_D_F01_rmv <- dplyr::filter(S.abs_D_F01_sum, FeatureTotal < 2)
S.abs_E_F01_rmv <- dplyr::filter(S.abs_E_F01_sum, FeatureTotal < 2)
S.abs_f_F01_rmv <- dplyr::filter(S.abs_f_F01_sum, FeatureTotal < 2)
S.abs_sEPD_M00_rmv <- dplyr::filter(S.abs_sEPD_M00_sum, FeatureTotal < 2)
S.abs_sPD_M00_rmv <- dplyr::filter(S.abs_sPD_M00_sum, FeatureTotal < 2)
S.abs_sE_M00_rmv <- dplyr::filter(S.abs_sE_M00_sum, FeatureTotal < 2)
S.abs_sf_M00_rmv <- dplyr::filter(S.abs_sf_M00_sum, FeatureTotal < 2)
S.abs_s_M00_rmv <- dplyr::filter(S.abs_s_M00_sum, FeatureTotal < 2)
S.abs_EPD_M00_rmv <- dplyr::filter(S.abs_EPD_M00_sum, FeatureTotal < 2)
S.abs_PD_M00_rmv <- dplyr::filter(S.abs_PD_M00_sum, FeatureTotal < 2)
S.abs_P_M00_rmv <- dplyr::filter(S.abs_P_M00_sum, FeatureTotal < 2)
S.abs_D_M00_rmv <- dplyr::filter(S.abs_D_M00_sum, FeatureTotal < 2)
S.abs_E_M00_rmv <- dplyr::filter(S.abs_E_M00_sum, FeatureTotal < 2)
S.abs_f_M00_rmv <- dplyr::filter(S.abs_f_M00_sum, FeatureTotal < 2)

S.abs_e1t2_F01_rmv <- dplyr::filter(S.abs_e1t2_F01_sum, FeatureTotal < 2)
S.abs_e1t3_F01_rmv <- dplyr::filter(S.abs_e1t3_F01_sum, FeatureTotal < 2)
S.abs_e1t4_F01_rmv <- dplyr::filter(S.abs_e1t4_F01_sum, FeatureTotal < 2)
S.abs_e1t5_F01_rmv <- dplyr::filter(S.abs_e1t5_F01_sum, FeatureTotal < 2)
S.abs_e1t6_F01_rmv <- dplyr::filter(S.abs_e1t6_F01_sum, FeatureTotal < 2)

S.abs_sEPD_F01_rmv$FeatureID <- as.character(S.abs_sEPD_F01_rmv$FeatureID)
S.abs_sPD_F01_rmv$FeatureID <- as.character(S.abs_sPD_F01_rmv$FeatureID)
S.abs_sE_F01_rmv$FeatureID <- as.character(S.abs_sE_F01_rmv$FeatureID)
S.abs_sf_F01_rmv$FeatureID <- as.character(S.abs_sf_F01_rmv$FeatureID)
S.abs_s_F01_rmv$FeatureID <- as.character(S.abs_s_F01_rmv$FeatureID)
S.abs_EPD_F01_rmv$FeatureID <- as.character(S.abs_EPD_F01_rmv$FeatureID)
S.abs_PD_F01_rmv$FeatureID <- as.character(S.abs_PD_F01_rmv$FeatureID)
S.abs_P_F01_rmv$FeatureID <- as.character(S.abs_P_F01_rmv$FeatureID)
S.abs_D_F01_rmv$FeatureID <- as.character(S.abs_D_F01_rmv$FeatureID)
S.abs_E_F01_rmv$FeatureID <- as.character(S.abs_E_F01_rmv$FeatureID)
S.abs_f_F01_rmv$FeatureID <- as.character(S.abs_f_F01_rmv$FeatureID)
S.abs_sEPD_M00_rmv$FeatureID <- as.character(S.abs_sEPD_M00_rmv$FeatureID)
S.abs_sPD_M00_rmv$FeatureID <- as.character(S.abs_sPD_M00_rmv$FeatureID)
S.abs_sE_M00_rmv$FeatureID <- as.character(S.abs_sE_M00_rmv$FeatureID)
S.abs_sf_M00_rmv$FeatureID <- as.character(S.abs_sf_M00_rmv$FeatureID)
S.abs_s_M00_rmv$FeatureID <- as.character(S.abs_s_M00_rmv$FeatureID)
S.abs_EPD_M00_rmv$FeatureID <- as.character(S.abs_EPD_M00_rmv$FeatureID)
S.abs_PD_M00_rmv$FeatureID <- as.character(S.abs_PD_M00_rmv$FeatureID)
S.abs_P_M00_rmv$FeatureID <- as.character(S.abs_P_M00_rmv$FeatureID)
S.abs_D_M00_rmv$FeatureID <- as.character(S.abs_D_M00_rmv$FeatureID)
S.abs_E_M00_rmv$FeatureID <- as.character(S.abs_E_M00_rmv$FeatureID)
S.abs_f_M00_rmv$FeatureID <- as.character(S.abs_f_M00_rmv$FeatureID)

S.abs_e1t2_F01_rmv$FeatureID <- as.character(S.abs_e1t2_F01_rmv$FeatureID)
S.abs_e1t3_F01_rmv$FeatureID <- as.character(S.abs_e1t3_F01_rmv$FeatureID)
S.abs_e1t4_F01_rmv$FeatureID <- as.character(S.abs_e1t4_F01_rmv$FeatureID)
S.abs_e1t5_F01_rmv$FeatureID <- as.character(S.abs_e1t5_F01_rmv$FeatureID)
S.abs_e1t6_F01_rmv$FeatureID <- as.character(S.abs_e1t6_F01_rmv$FeatureID)

S.abs_sEPD_F01 <- dplyr::anti_join(x = S.abs_sEPD_F01_0, y = S.abs_sEPD_F01_rmv, 
                                   by = "FeatureID")
S.abs_sPD_F01 <- dplyr::anti_join(x = S.abs_sPD_F01_0, y = S.abs_sPD_F01_rmv, 
                                  by = "FeatureID")
S.abs_sE_F01 <- dplyr::anti_join(x = S.abs_sE_F01_0, y = S.abs_sE_F01_rmv, 
                                 by = "FeatureID")
S.abs_sf_F01 <- dplyr::anti_join(x = S.abs_sf_F01_0, y = S.abs_sf_F01_rmv, 
                                 by = "FeatureID")
S.abs_s_F01 <- dplyr::anti_join(x = S.abs_s_F01_0, y = S.abs_s_F01_rmv, 
                                by = "FeatureID")
S.abs_EPD_F01 <- dplyr::anti_join(x = S.abs_EPD_F01_0, y = S.abs_EPD_F01_rmv, 
                                  by = "FeatureID")
S.abs_PD_F01 <- dplyr::anti_join(x = S.abs_PD_F01_0, y = S.abs_PD_F01_rmv, 
                                 by = "FeatureID")
S.abs_P_F01 <- dplyr::anti_join(x = S.abs_P_F01_0, y = S.abs_P_F01_rmv, 
                                by = "FeatureID")
S.abs_D_F01 <- dplyr::anti_join(x = S.abs_D_F01_0, y = S.abs_D_F01_rmv, 
                                by = "FeatureID")
S.abs_E_F01 <- dplyr::anti_join(x = S.abs_E_F01_0, y = S.abs_E_F01_rmv, 
                                by = "FeatureID")
S.abs_f_F01 <- dplyr::anti_join(x = S.abs_f_F01_0, y = S.abs_f_F01_rmv, 
                                by = "FeatureID")
S.abs_sEPD_M00 <- dplyr::anti_join(x = S.abs_sEPD_M00_0, y = S.abs_sEPD_M00_rmv, 
                                   by = "FeatureID")
S.abs_sPD_M00 <- dplyr::anti_join(x = S.abs_sPD_M00_0, y = S.abs_sPD_M00_rmv, 
                                  by = "FeatureID")
S.abs_sE_M00 <- dplyr::anti_join(x = S.abs_sE_M00_0, y = S.abs_sE_M00_rmv, 
                                 by = "FeatureID")
S.abs_sf_M00 <- dplyr::anti_join(x = S.abs_sf_M00_0, y = S.abs_sf_M00_rmv, 
                                 by = "FeatureID")
S.abs_s_M00 <- dplyr::anti_join(x = S.abs_s_M00_0, y = S.abs_s_M00_rmv, 
                                by = "FeatureID")
S.abs_EPD_M00 <- dplyr::anti_join(x = S.abs_EPD_M00_0, y = S.abs_EPD_M00_rmv, 
                                  by = "FeatureID")
S.abs_PD_M00 <- dplyr::anti_join(x = S.abs_PD_M00_0, y = S.abs_PD_M00_rmv, 
                                 by = "FeatureID")
S.abs_P_M00 <- dplyr::anti_join(x = S.abs_P_M00_0, y = S.abs_P_M00_rmv, 
                                by = "FeatureID")
S.abs_D_M00 <- dplyr::anti_join(x = S.abs_D_M00_0, y = S.abs_D_M00_rmv, 
                                by = "FeatureID")
S.abs_E_M00 <- dplyr::anti_join(x = S.abs_E_M00_0, y = S.abs_E_M00_rmv, 
                                by = "FeatureID")
S.abs_f_M00 <- dplyr::anti_join(x = S.abs_f_M00_0, y = S.abs_f_M00_rmv, 
                                by = "FeatureID")

S.abs_e1t2_F01 <- dplyr::anti_join(x = S.abs_e1t2_F01_0, y = S.abs_e1t2_F01_rmv, 
                                   by = "FeatureID")
S.abs_e1t3_F01 <- dplyr::anti_join(x = S.abs_e1t3_F01_0, y = S.abs_e1t4_F01_rmv, 
                                   by = "FeatureID")
S.abs_e1t4_F01 <- dplyr::anti_join(x = S.abs_e1t4_F01_0, y = S.abs_e1t4_F01_rmv, 
                                   by = "FeatureID")
S.abs_e1t5_F01 <- dplyr::anti_join(x = S.abs_e1t5_F01_0, y = S.abs_e1t5_F01_rmv, 
                                   by = "FeatureID")
S.abs_e1t6_F01 <- dplyr::anti_join(x = S.abs_e1t6_F01_0, y = S.abs_e1t6_F01_rmv, 
                                   by = "FeatureID")

### ************************************
### S - STEP 2 - information gathering and provenance ----
### ************************************

# ** note for KDP: code does not reflect same process used in: 
# SS subset EBTKS - version 0.1
# OR
# *SS information gathering and provenance version 0.1

# (1) bind S.info data.frames together
# (2) add info/script/section/heading columns
# (3) reorder columns

# S.info_0 <- rbind(S.info_abs_bse, S.info_abs_end, S.info_abs_CF, S.info_abs_RN)
# 
# S.info_1 <- S.info_0
# S.info_1$info <- S.prov_info
# S.info_1$script <- name_scrpt
# S.info_1$section <- S.prov_secstep_SS1
# S.info_1$heading <- S.prov_heading_SS1
# 
# S.info <- dplyr::select(S.info_1, info, value, object, script, section, heading)

### ************************************
### S - WRITE OUTPUTS ----
### ************************************

# outputs to the vault
# write.table(sep = "\t", row.names = F, x = S.info, file = S.ofv_info)

S.obj <- ls(pattern = "S.")
S.lst <- c(S.obj[grep(pattern = "S.", x = S.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS, S.obj_from_A)
save(list = S.lst, file = S.ofv_wksp)

### ************************************
### T - STEP 1 - create a rooted phylogenetic tree from processed EBTKS ----
### ************************************

# ** note for KDP: TS version 0.4 ** #

# provide provenance for information gathering at end of section:
T.prov_secstep_TS1 <- "Section T - STEP 1"
T.prov_heading_TS1 <- "create a rooted phylogenetic tree from processed EBTKS"
T.prov_output_obj_TS1 <- "T.tree_root_GTR" # this object is output to the vault

# NOTE: tree construction follows the code structure in:
# Microbiome Data Analysis: from raw reads to community analyses.
# created by Benjamin J Callahan, Kris Sankaran, Julia A Fukuyama,
# Paul Joey McMurdie, and Susan P Holmes

# NOTE: if the environment is empty; there are some requirements:
# if section A has been run, then run the PREFACE and load section A's workspace
# if section A has not been run, then run everything above this line
# ** note for KDP:
# the choose your own adventure function that allows sections to standalone ...
# ... or stitches them together is still under construction ** #

# create a new version of the object needed from section A ...
# ... and create a vector naming that object (used when saving the workspace)
T.EBTKS_abs_pro <- A.EBTKS_abs_pro
T.obj_from_A <- "A.EBTKS_abs_pro"

# the tree construction process occurs as follows:
# (1) create a vector of RepSeqs
# (2) create a copy to avoid overwriting the original
# (3) assign RepSeqs as explicit names
# (4) align RepSeqs
# (5) convert alignment to class 'phyDat'
# (6) compute pairwise distances from DNA sequences
# (7) construct a neighbor-joining tree
# (8:10) fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree
# (11) extract the unrooted tree
# (12) create a midpoint rooted tree
# NOTE: GTR+G+I = Generalized time-reversible with Gamma rate variation
# runtime (1:12): 1 to 1.5 hours

T.RepSeq_0 <- T.EBTKS_abs_pro$RepSeq
T.RepSeq_1 <- T.RepSeq_0
names(T.RepSeq_1) <- T.RepSeq_1
T.RepSeq_algn_0 <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(T.RepSeq_1),
                                       anchor = NA, verbose = F,
                                       processors = NULL)
T.RepSeq_algn_1 <- phangorn::phyDat(as(T.RepSeq_algn_0, "matrix"), type = "DNA")
T.dist_ml <- phangorn::dist.ml(T.RepSeq_algn_1, processors = NULL)
T.tree_NJ <- phangorn::NJ(T.dist_ml)
T.fit = phangorn::pml(T.tree_NJ, data = T.RepSeq_algn_1, processors = NULL)
T.fitGTR <- update(T.fit, k = 4, inv = 0.2)
T.fitGTR <- phangorn::optim.pml(T.fitGTR, model = "GTR", optInv = T,
                                optGamma = T, rearrangement = "stochastic",
                                control = phangorn::pml.control(trace = 0))
T.tree_GTR <- T.fitGTR$tree
T.tree_root_GTR <- phangorn::midpoint(tree = T.tree_GTR)

### ************************************
### T - WRITE OUTPUTS ----
### ************************************

# provenance for output to the vault
T.prov <- data.frame("info" = "provenance for output",
                     "path" = T.ofv_EBTKS_tree,
                     "object" = T.prov_output_obj_TS1,
                     "script" = name_scrpt,
                     "section" = T.prov_secstep_TS1,
                     "heading" = T.prov_heading_TS1,
                     stringsAsFactors = F)

# outputs to the vault
ape::write.tree(phy = T.tree_root_GTR, file = T.ofv_EBTKS_tree)

T.obj <- ls(pattern = "T.")
T.lst <- c(T.obj[grep(pattern = "T.", x = T.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, T.obj_from_A)
save(list = T.lst, file = T.ofv_wksp)
