### ************************************
### IMPORTANT SCRIPT NOTE ----
### ************************************

# this file will not 'source' and may be difficult to process through
# all of the code is here, but I have not yet made it 'easily' reprocudible
# following review and changes that may need to be made to figs/tabs/etc...
# following acceptance of the paper for publication...
# ... this will all be organized accordingly.
# for an older style of what the organization will look like, see:
# https://github.com/kdprkr/MerlinsManuscript

### ************************************
### BEGIN Section B ----
### ************************************

# note for kdp: taken from TransFaunation/section_B.R

#### preface ----

# R packages accessed via require:
require(dendextend, quietly = T)
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)
require(ggraph, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
ifp_EBTKS_tree <- paste(sep = "", path_vault, "/tree_EBTKS_processed.nwk")
B.tree_root <- ape::read.tree(ifp_EBTKS_tree)
load(file = S.ofv_wksp)
load(file = ofv_COSMOS_wksp)

B.ofv_wksp <- paste(sep = "", path_vault, "/WS_", "SectionB.RData")

### ************************************
### B - FUNCTIONS ----
### ************************************

# create a vector naming all Section A functions (used when saving workspace)
B.functions <- c("iso_var_exp", "iso_cord", "flt_smp_dst", "pwise_adon")

# ** note for KDP: version 0.2 ** #
# function takes an input of class "pcoa" or class "prcomp" and isolates...
# the amount of variation explained by the PC specified
# NOTE: isolates 1 value at a time as specified by pc.num
# NOTE: returns % values with 3 digits
# iso_var_exp = isolate_variation_explained
iso_var_exp <- function(data = c(pcoa, prcomp), pc.num, 
                        output = c("value", "label")) {
  # internal check to ensure correct input class
  if (!inherits(data, c("pcoa", "prcomp"))) {
    stop("data must be class 'pcoa' OR class 'prcomp'")
  }
  if (!inherits(pc.num, "numeric")) {
    stop("pc.num must be numeric")
  }
  if (!output == "value" && !output == "label") {
    stop("output must be one of: 'value' OR 'label'")
  }
  if (inherits(data, "pcoa")) {
    # isolate pcoa$values (create a data.frame)
    pcoa_vals <- data.frame(data$values)
    # isolate the Relative_eig values
    rel_eigs <- pcoa_vals$Relative_eig
    # convert rel eigs to percentage and round to 3 digits
    rel_eigs_rnd <- (round(rel_eigs, digits = 3)) * 100
    # isolate the specified PC
    pcoa_val <- rel_eigs_rnd[pc.num]
    # statements to determine what to return
    if (output == "value") {
      return(pcoa_val)
    }
    if (output == "label") {
      pcoa_lab <- paste("PC ", pc.num, " (", pcoa_val, "%)", sep = "")
      return(pcoa_lab)
    }
  }
  if (inherits(data, "prcomp")) {
    # isolate standard deviations for all pcs
    pcs_sdv <- data$sdev
    # square each sdev, sum the squares, and divide each sdev by the sum
    pcs_sdv_sqr <- pcs_sdv^2
    pcs_sdv_sqr_sum <- sum(pcs_sdv_sqr)
    pca_vals_raw <- pcs_sdv_sqr/pcs_sdv_sqr_sum
    # convert to percentage and round to 3 digits
    pca_vals_rnd <- (round(pca_vals_raw, digits = 3)) * 100
    # isolate the specified PC
    pca_val <- pca_vals_rnd[pc.num]
    # statements to determine what to return
    if (output == "value") {
      return(pca_val)
    }
    if (output == "label") {
      pca_lab <- paste("PC ", pc.num, " (", pca_val, "%)", sep = "")
      return(pca_lab)
    }
  }
}
# 
# # example usage:
# val <- iso_var_exp(data = pcoa, pc.num = 1, output = "value") # returns number
# lab <- iso_var_exp(data = pcoa, pc.num = 1, output = "label") # returns label

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.1 ** #
# function takes an input of class "pcoa" isolates...
# the x and y coordinates for each sample
# iso_cord = isolate_coordinates
iso_cord <- function(data = pcoa) {
  # internal check to ensure correct input class
  if (!inherits(data, "pcoa")) {
    stop("data must be class 'pcoa'")
  }
  # isolate pcoa$vectors (create a data.frame)
  pcoa_vec <- data.frame(data$vectors)
  # isolate cols 1 and 2 (Axis.1, Axis.2; i.e. x and y coordinates)
  pcoa_crd <- pcoa_vec[, 1:2]
  # rename cols
  pcoa_crd_rnm <- pcoa_crd
  names(pcoa_crd_rnm)[1] <- "X.Coord"
  names(pcoa_crd_rnm)[2] <- "Y.Coord"
  return(pcoa_crd_rnm)
}
#
# example usage:
# pcoa_cord_df <- iso_cord(data = pcoa)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.1 ** #
# function takes inputs of sampledata and a distance matrix and...
# ... returns filtered versions of both output as a list
# the resultant outputs contain identical samples which is useful when...
# ... running PERMANOVA/adonis or using the pwise_adon() function
# NOTE: this function is a dependency for pwise_adon()
flt_smp_dst <- function(smp_dat = data.frame, dst_obj = dist, flt_col = "", 
                        flt_row = "", match_col = "") {
  # internal checks to ensure correct input classes
  if (!inherits(smp_dat, "data.frame")) {
    stop("input smp_dat must be class 'data.frame'")
  }
  if (!inherits(dst_obj, "dist")) {
    stop("input dst_obj must be class 'dist'")
  }
  if (!inherits(flt_col, "character")) {
    stop("flt_col must be character")
  }
  if (!inherits(flt_row, "character")) {
    stop("flt_row must be character")
  }
  if (!inherits(match_col, "character")) {
    stop("match_col must be character")
  }
  # create vector of row #s in sampledata flt_col that match the flt_row
  row_num <- which(smp_dat[, flt_col] %in% flt_row)
  if (isTRUE(length(row_num) == 0)) {
    stop("`flt_row` not found in `flt_col`")
  }
  # create new sampledata df using the above vector
  new_smp_dat_0 <- smp_dat[row_num, ]
  # coerce distance matrix from class: 'dist' to class: 'matrix'
  dmx <- as.matrix(dst_obj)
  # create a vector of row #s in the dmx that match column in the sampledata
  dmx_vec_nums <- which(row.names(dmx) %in% new_smp_dat_0[, match_col])
  # filter the distance matrix to retain only the numbers in the above vector
  # NOTE: we can take advantage of matrixland as col #s and row #s are identical
  new_dmx <- dmx[dmx_vec_nums, dmx_vec_nums]
  # match the order rows in new sampledata and new distance matrix
  new_smp_dat_1 <- new_smp_dat_0[match(x = row.names(new_dmx), 
                                       table = new_smp_dat_0[, match_col]), ]
  # convert new_dmx back into class 'dist'
  new_dst <- as.dist(new_dmx)
  # create a list containing new_smp_dat and new.dist
  new_lst <- list(new_smp = new_smp_dat_1, new_dst = new_dst)
  #return(new_lst)
  return(new_lst)
}
#
# example usage:
# filtered_list <- flt_smp_dst(smp_dat = smp, dst_obj = dst, flt_col = "group",
#                              flt_row = "housecat", match_col = "SampleID")
# filtered_smp_dat <- filtered_list$new_smp
# filtered_dst_obj <- filtered_list$new_dst

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.3 ** #
# function takes inputs of sampledata and a distance matrix and...
# ... performs pairwise PERMANOVA comparisons for all possible pairs...
# ... present in the specified input 'test_col' and returns a data.frame with...
# ... the comparison, R2 value, unadjusted p-values & p values adjusted...
# ... using the specified adjustment method
# NOTE: this function is dependent on the flt_smp_dst() function
pwise_adon <- function(smp_dat = data.frame, dst_obj = dist, test_col = "", 
                       match_col = "", perms, digits,
                       adjust = c("holm", "hochberg", "hommel", "bonferroni", 
                                  "BH", "BY", "fdr", "none")) {
  # internal checks to ensure correct input classes
  if (!inherits(smp_dat, "data.frame")) {
    stop("input smp_dat must be class 'data.frame'")
  }
  if (!inherits(dst_obj, "dist")) {
    stop("input dst_obj must be class 'dist'")
  }
  if (!inherits(test_col, "character")) {
    stop("test_col must be character")
  }
  if (!inherits(match_col, "character")) {
    stop("match_col must be character")
  }
  if (!adjust == "holm" && !adjust == "hochberg" && !adjust == "hommel" &&
      !adjust == "bonferroni" && !adjust == "BH" && !adjust == "BY" &&
      !adjust == "fdr" && !adjust == "none") {
    stop("adjust must be one of the methods supported by p.adjust();
    'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'")
  }
  # create a vector of row values in sampledata flt_col
  row_val <- smp_dat[, test_col]
  # create unique factors for row values in sampledata flt_col
  factors <- as.factor(unique(row_val))
  # obtain all pairwise comparisons from the above row values
  comps <- combn(x = unique(row_val), m = 2)
  # create empty vectors to store values obtained in the for loop below
  pair_labs <- c()
  R2 <- c()
  p_raw <- c()
  # loop through the data to make adonis comparisons for each pair in comps
  for(i in 1:ncol(comps)){
    # create a vector for an individual pair
    pair <- as.vector(factors[factors %in% c(comps[1, i], comps[2, i])])
    # filter the sampledata df and the distance matrix for the paired comparison
    filt <- flt_smp_dst(flt_row = pair, smp_dat = smp_dat, dst_obj = dst_obj, 
                        flt_col = test_col, match_col = match_col)
    # create filtered data
    flt_smp <- filt$new_smp
    # create new distance matrix reflecting the filtering step above
    flt_dst <- filt$new_dst
    # perform adonis on the filtered data
    adon <- vegan::adonis(formula = flt_dst ~ flt_smp[, test_col], 
                          data = flt_smp, permutations = perms)
    # define labels for the individual pair
    pair_labs <- c(pair_labs, paste(pair[1], 'vs', pair[2]))
    # store information from the adonis data.frame
    R2 <- c(R2, adon$aov.tab[1, 5])
    p_raw <- c(p_raw, adon$aov.tab[1, 6])
  }
  # adjust raw p values for multiple comparisons
  # NOTE: number of comparisons (n) is determined by length(p_raw)
  p_adj <- p.adjust(p = p_raw, method = adjust, n = length(p_raw))
  # create new data.frame containing the relvant information
  pwise_df_0 <- data.frame("Comparison" = pair_labs, R2, p_raw, p_adj)
  # append input adjustment method to column name for p_adj
  pwise_df_1 <- pwise_df_0
  p_adj_num <- which(names(pwise_df_1) == "p_adj")
  names(pwise_df_1)[p_adj_num] <- paste("p_", adjust, sep = "")
  # round all numeric values to specified number of decimal places
  pwise_df_rnd <- dplyr::mutate_if(pwise_df_1, is.numeric, 
                                   list(~ round(., digits)))
  return(pwise_df_rnd)
}
#
# example usage:
# results_pwise_adon <- pwise_adon(smp_dat = smp, dst_obj = dst,
#                                  test_col = "group", match_col = "SampleID",
#                                  perms = 999, adjust = "BH", digits = 3)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

### ************************************
### ^^^B - STEP  1 - format data for UniFrac/PCoA & CZM/clr/PCA/clustering ----
### ************************************

# provide provenance for information gathering at end of section:
B.prov_secstep_BS1 <- "Section B - STEP 1"
B.prov_heading_BS1 <- "format data for UniFrac/PCoA & CZM/clr/PCA/clustering"
# ^^^B.prov_output_obj_BS1 <- "" # this object is output to the vault'
# ^^^info and prov have not actually been added to this script

# ^^^NOTE: if the environment is empty; there are some requirements:
# if section A has been run, then run the PREFACE and load section A's workspace
# if section A has not been run, then run everything above this line
# ** note for KDP:
# the choose your own adventure function that allows sections to standalone ...
# ... or stitches them together is still under construction ** #

# create a new version of the objects needed from sections S; T ...
# ... and create a vector naming those objects (used when saving the workspace)
T.tree_root_GTR <- B.tree_root # ^^^this line serves the B.obj_from_ST
B.abs_sEPD_F01_0 <- S.abs_sEPD_F01
B.abs_EPD_F01_0 <- S.abs_EPD_F01
B.abs_s_F01_0 <- S.abs_s_F01
B.abs_sEPD_M00_0 <- S.abs_sEPD_M00
B.abs_EPD_M00_0 <- S.abs_EPD_M00
B.abs_s_M00_0 <- S.abs_s_M00
B.obj_from_ST <- c("T.tree_root_GTR", "S.abs_sEPD_F01", "S.abs_EPD_F01", 
                   "S.abs_s_F01", "S.abs_sEPD_M00", "S.abs_EPD_M00", 
                   "S.abs_s_M00")

# format data to obtain RepSeq as columns and SampleIDs as rows
# this process occurs as follows:
# (1) create copies to avoid overwriting the originals
# (2) convert column RepSeq into row.names
# (3) remove unneeded columns
# (4) transpose to convert features to cols and samples to rows

B.abs_sEPD_F01_1 <- B.abs_sEPD_F01_0
B.abs_EPD_F01_1 <- B.abs_EPD_F01_0
B.abs_s_F01_1 <- B.abs_s_F01_0
B.abs_sEPD_M00_1 <- B.abs_sEPD_M00_0
B.abs_EPD_M00_1 <- B.abs_EPD_M00_0
B.abs_s_M00_1 <- B.abs_s_M00_0

row.names(B.abs_sEPD_F01_1) <- B.abs_sEPD_F01_1$RepSeq
row.names(B.abs_EPD_F01_1) <- B.abs_EPD_F01_1$RepSeq
row.names(B.abs_s_F01_1) <- B.abs_s_F01_1$RepSeq
row.names(B.abs_sEPD_M00_1) <- B.abs_sEPD_M00_1$RepSeq
row.names(B.abs_EPD_M00_1) <- B.abs_EPD_M00_1$RepSeq
row.names(B.abs_s_M00_1) <- B.abs_s_M00_1$RepSeq

B.abs_sEPD_F01_2 <- dplyr::select(B.abs_sEPD_F01_1, -dplyr::one_of(com_col))
B.abs_EPD_F01_2 <- dplyr::select(B.abs_EPD_F01_1, -dplyr::one_of(com_col))
B.abs_s_F01_2 <- dplyr::select(B.abs_s_F01_1, -dplyr::one_of(com_col))
B.abs_sEPD_M00_2 <- dplyr::select(B.abs_sEPD_M00_1, -dplyr::one_of(com_col))
B.abs_EPD_M00_2 <- dplyr::select(B.abs_EPD_M00_1, -dplyr::one_of(com_col))
B.abs_s_M00_2 <- dplyr::select(B.abs_s_M00_1, -dplyr::one_of(com_col))

B.abs_sEPD_F01_3 <- as.data.frame(t(B.abs_sEPD_F01_2))
B.abs_EPD_F01_3 <- as.data.frame(t(B.abs_EPD_F01_2))
B.abs_s_F01_3 <- as.data.frame(t(B.abs_s_F01_2))
B.abs_sEPD_M00_3 <- as.data.frame(t(B.abs_sEPD_M00_2))
B.abs_EPD_M00_3 <- as.data.frame(t(B.abs_EPD_M00_2))
B.abs_s_M00_3 <- as.data.frame(t(B.abs_s_M00_2))

### ************************************
### B - STEP  2 - compute UniFrac/PCoA ----
### ************************************

# this process occurs as follows:
# (1) compute UniFracs
# (2) extract relevant metrics from the unifrac lists
# (3) compute PCoA
# (4) create PCoA plot axis labels (i.e. isolate variation explained by PC1;PC2)
# (5) isolate x & y coordinates
# (6;7;8) convert row.names (i.e. SampleIDs) into a col & merge with sample data
# NOTE: for uni, input data = features as cols and samples as rows

B.uni_sEPD_F01 <- GUniFrac::GUniFrac(otu.tab = B.abs_sEPD_F01_3, 
                                     tree = B.tree_root, 
                                     alpha = c(0, 0.5, 1))$unifracs
B.uni_EPD_F01 <- GUniFrac::GUniFrac(otu.tab = B.abs_EPD_F01_3, 
                                    tree = B.tree_root, 
                                    alpha = c(0, 0.5, 1))$unifracs
B.uni_s_F01 <- GUniFrac::GUniFrac(otu.tab = B.abs_s_F01_3, 
                                  tree = B.tree_root, 
                                  alpha = c(0, 0.5, 1))$unifracs
B.uni_sEPD_M00 <- GUniFrac::GUniFrac(otu.tab = B.abs_sEPD_M00_3, 
                                     tree = B.tree_root, 
                                     alpha = c(0, 0.5, 1))$unifracs
B.uni_EPD_M00 <- GUniFrac::GUniFrac(otu.tab = B.abs_EPD_M00_3, 
                                    tree = B.tree_root, 
                                    alpha = c(0, 0.5, 1))$unifracs
B.uni_s_M00 <- GUniFrac::GUniFrac(otu.tab = B.abs_s_M00_3, 
                                  tree = B.tree_root, 
                                  alpha = c(0, 0.5, 1))$unifracs

B.gun_sEPD_F01 <- as.dist(B.uni_sEPD_F01[, , "d_0.5"])
B.uun_sEPD_F01 <- as.dist(B.uni_sEPD_F01[, , "d_UW"])
B.wun_sEPD_F01 <- as.dist(B.uni_sEPD_F01[, , "d_1"])
B.gun_EPD_F01 <- as.dist(B.uni_EPD_F01[, , "d_0.5"])
B.uun_EPD_F01 <- as.dist(B.uni_EPD_F01[, , "d_UW"])
B.wun_EPD_F01 <- as.dist(B.uni_EPD_F01[, , "d_1"])
B.gun_s_F01 <- as.dist(B.uni_s_F01[, , "d_0.5"])
B.uun_s_F01 <- as.dist(B.uni_s_F01[, , "d_UW"])
B.wun_s_F01 <- as.dist(B.uni_s_F01[, , "d_1"])
B.gun_sEPD_M00 <- as.dist(B.uni_sEPD_M00[, , "d_0.5"])
B.uun_sEPD_M00 <- as.dist(B.uni_sEPD_M00[, , "d_UW"])
B.wun_sEPD_M00 <- as.dist(B.uni_sEPD_M00[, , "d_1"])
B.gun_EPD_M00 <- as.dist(B.uni_EPD_M00[, , "d_0.5"])
B.uun_EPD_M00 <- as.dist(B.uni_EPD_M00[, , "d_UW"])
B.wun_EPD_M00 <- as.dist(B.uni_EPD_M00[, , "d_1"])
B.gun_s_M00 <- as.dist(B.uni_s_M00[, , "d_0.5"])
B.uun_s_M00 <- as.dist(B.uni_s_M00[, , "d_UW"])
B.wun_s_M00 <- as.dist(B.uni_s_M00[, , "d_1"])

B.pcoa_gun_sEPD_F01 <- ape::pcoa(B.gun_sEPD_F01)
B.pcoa_uun_sEPD_F01 <- ape::pcoa(B.uun_sEPD_F01)
B.pcoa_wun_sEPD_F01 <- ape::pcoa(B.wun_sEPD_F01)
B.pcoa_gun_sEPD_M00 <- ape::pcoa(B.gun_sEPD_M00)
B.pcoa_uun_sEPD_M00 <- ape::pcoa(B.uun_sEPD_M00)
B.pcoa_wun_sEPD_M00 <- ape::pcoa(B.wun_sEPD_M00)

B.pc1_gun_sEPD_F01 <- iso_var_exp(data = B.pcoa_gun_sEPD_F01, pc.num = 1, 
                                  output = "label")
B.pc2_gun_sEPD_F01 <- iso_var_exp(data = B.pcoa_gun_sEPD_F01, pc.num = 2, 
                                  output = "label")
B.pc1_uun_sEPD_F01 <- iso_var_exp(data = B.pcoa_uun_sEPD_F01, pc.num = 1,
                                  output = "label")
B.pc2_uun_sEPD_F01 <- iso_var_exp(data = B.pcoa_uun_sEPD_F01, pc.num = 2,
                                  output = "label")
B.pc1_wun_sEPD_F01 <- iso_var_exp(data = B.pcoa_wun_sEPD_F01, pc.num = 1,
                                  output = "label")
B.pc2_wun_sEPD_F01 <- iso_var_exp(data = B.pcoa_wun_sEPD_F01, pc.num = 2,
                                  output = "label")
B.pc1_gun_sEPD_M00 <- iso_var_exp(data = B.pcoa_gun_sEPD_M00, pc.num = 1,
                                  output = "label")
B.pc2_gun_sEPD_M00 <- iso_var_exp(data = B.pcoa_gun_sEPD_M00, pc.num = 2,
                                  output = "label")
B.pc1_uun_sEPD_M00 <- iso_var_exp(data = B.pcoa_uun_sEPD_M00, pc.num = 1,
                                  output = "label")
B.pc2_uun_sEPD_M00 <- iso_var_exp(data = B.pcoa_uun_sEPD_M00, pc.num = 2,
                                  output = "label")
B.pc1_wun_sEPD_M00 <- iso_var_exp(data = B.pcoa_wun_sEPD_M00, pc.num = 1,
                                  output = "label")
B.pc2_wun_sEPD_M00 <- iso_var_exp(data = B.pcoa_wun_sEPD_M00, pc.num = 2,
                                  output = "label")

B.cord_gun_sEPD_F01_0 <- iso_cord(data = B.pcoa_gun_sEPD_F01)
B.cord_uun_sEPD_F01_0 <- iso_cord(data = B.pcoa_uun_sEPD_F01)
B.cord_wun_sEPD_F01_0 <- iso_cord(data = B.pcoa_wun_sEPD_F01)
B.cord_gun_sEPD_M00_0 <- iso_cord(data = B.pcoa_gun_sEPD_M00)
B.cord_uun_sEPD_M00_0 <- iso_cord(data = B.pcoa_uun_sEPD_M00)
B.cord_wun_sEPD_M00_0 <- iso_cord(data = B.pcoa_wun_sEPD_M00)

B.cord_gun_sEPD_F01_1 <- B.cord_gun_sEPD_F01_0
B.cord_uun_sEPD_F01_1 <- B.cord_uun_sEPD_F01_0
B.cord_wun_sEPD_F01_1 <- B.cord_wun_sEPD_F01_0
B.cord_gun_sEPD_M00_1 <- B.cord_gun_sEPD_M00_0
B.cord_uun_sEPD_M00_1 <- B.cord_uun_sEPD_M00_0
B.cord_wun_sEPD_M00_1 <- B.cord_wun_sEPD_M00_0
B.cord_gun_sEPD_F01_1$SampleID <- row.names(B.cord_gun_sEPD_F01_1)
B.cord_uun_sEPD_F01_1$SampleID <- row.names(B.cord_uun_sEPD_F01_1)
B.cord_wun_sEPD_F01_1$SampleID <- row.names(B.cord_wun_sEPD_F01_1)
B.cord_gun_sEPD_M00_1$SampleID <- row.names(B.cord_gun_sEPD_M00_1)
B.cord_uun_sEPD_M00_1$SampleID <- row.names(B.cord_uun_sEPD_M00_1)
B.cord_wun_sEPD_M00_1$SampleID <- row.names(B.cord_wun_sEPD_M00_1)
B.smp_cord_gun_sEPD_F01 <- merge(x = B.cord_gun_sEPD_F01_1, y = smp_dat, 
                                 all = F, sort = F, by = "SampleID")
B.smp_cord_uun_sEPD_F01 <- merge(x = B.cord_uun_sEPD_F01_1, y = smp_dat, 
                                 all = F, sort = F, by = "SampleID")
B.smp_cord_wun_sEPD_F01 <- merge(x = B.cord_wun_sEPD_F01_1, y = smp_dat, 
                                 all = F, sort = F, by = "SampleID")
B.smp_cord_gun_sEPD_M00 <- merge(x = B.cord_gun_sEPD_M00_1, y = smp_dat, 
                                 all = F, sort = F, by = "SampleID")
B.smp_cord_uun_sEPD_M00 <- merge(x = B.cord_uun_sEPD_M00_1, y = smp_dat, 
                                 all = F, sort = F, by = "SampleID")
B.smp_cord_wun_sEPD_M00 <- merge(x = B.cord_wun_sEPD_M00_1, y = smp_dat, 
                                 all = F, sort = F, by = "SampleID")

### ************************************
### B - STEP  3 - compute CZM/clr/PCA/clustering ----
### ************************************

# for CZM and clr, the process occurs as follows:
# (1) replace zero counts
# (2) check proportion conversion
# (3) transpose to convert features to rows and samples to cols
# (4) clr transform
# NOTE: for cmultRepl(), input data = features as cols and samples as rows
# NOTE: for clr transform, input data = features as rows and samples as cols

B.czm_sEPD_F01_0 <- zCompositions::cmultRepl(B.abs_sEPD_F01_3, method = "CZM", 
                                             output = "prop")
B.czm_EPD_F01_0 <- zCompositions::cmultRepl(B.abs_EPD_F01_3, method = "CZM", 
                                            output = "prop")
B.czm_s_F01_0 <- zCompositions::cmultRepl(B.abs_s_F01_3, method = "CZM", 
                                          output = "prop")
B.czm_sEPD_M00_0 <- zCompositions::cmultRepl(B.abs_sEPD_M00_3, method = "CZM", 
                                             output = "prop")
B.czm_EPD_M00_0 <- zCompositions::cmultRepl(B.abs_EPD_M00_3, method = "CZM", 
                                            output = "prop")
B.czm_s_M00_0 <- zCompositions::cmultRepl(B.abs_s_M00_3, method = "CZM", 
                                          output = "prop")

print(rowSums(B.czm_sEPD_F01_0)) # all == 1
print(rowSums(B.czm_EPD_F01_0)) # all == 1
print(rowSums(B.czm_s_F01_0)) # all == 1
print(rowSums(B.czm_sEPD_M00_0)) # all == 1
print(rowSums(B.czm_EPD_M00_0)) # all == 1
print(rowSums(B.czm_s_M00_0)) # all == 1

B.czm_sEPD_F01_1 <- as.data.frame(t(B.czm_sEPD_F01_0))
B.czm_EPD_F01_1 <- as.data.frame(t(B.czm_EPD_F01_0))
B.czm_s_F01_1 <- as.data.frame(t(B.czm_s_F01_0))
B.czm_sEPD_M00_1 <- as.data.frame(t(B.czm_sEPD_M00_0))
B.czm_EPD_M00_1 <- as.data.frame(t(B.czm_EPD_M00_0))
B.czm_s_M00_1 <- as.data.frame(t(B.czm_s_M00_0))

B.clr_sEPD_F01_0 <- as.data.frame(apply(B.czm_sEPD_F01_1, 2, 
                                        function(x) {log2(x) - mean(log2(x))}))
B.clr_EPD_F01_0 <- as.data.frame(apply(B.czm_EPD_F01_1, 2, 
                                       function(x) {log2(x) - mean(log2(x))}))
B.clr_s_F01_0 <- as.data.frame(apply(B.czm_s_F01_1, 2, 
                                     function(x) {log2(x) - mean(log2(x))}))
B.clr_sEPD_M00_0 <- as.data.frame(apply(B.czm_sEPD_M00_1, 2, 
                                        function(x) {log2(x) - mean(log2(x))}))
B.clr_EPD_M00_0 <- as.data.frame(apply(B.czm_EPD_M00_1, 2, 
                                       function(x) {log2(x) - mean(log2(x))}))
B.clr_s_M00_0 <- as.data.frame(apply(B.czm_s_M00_1, 2, 
                                     function(x) {log2(x) - mean(log2(x))}))

# for PCA, the process occurs as follows:
# (1) transpose to convert features to cols and samples to rows
# (2) compute PCA
# (3) create PCA plot axis labels (i.e. isolate variation explained by PC1;PC2)
# (4;5;6) create ggbiplots, extract x & y coordinates, merge with sample data
# NOTE: for prcomp(), input data = features as cols and samples as rows

B.clr_sEPD_F01_1 <- as.data.frame(t(B.clr_sEPD_F01_0))
B.clr_sEPD_M00_1 <- as.data.frame(t(B.clr_sEPD_M00_0))

B.pca_sEPD_F01 <- prcomp(B.clr_sEPD_F01_1)
B.pca_sEPD_M00 <- prcomp(B.clr_sEPD_M00_1)

B.pc1_pca_sEPD_F01 <- iso_var_exp(data = B.pca_sEPD_F01, pc.num = 1, 
                                  output = "label")
B.pc2_pca_sEPD_F01 <- iso_var_exp(data = B.pca_sEPD_F01, pc.num = 2, 
                                  output = "label")
B.pc1_pca_sEPD_M00 <- iso_var_exp(data = B.pca_sEPD_M00, pc.num = 1, 
                                  output = "label")
B.pc2_pca_sEPD_M00 <- iso_var_exp(data = B.pca_sEPD_M00, pc.num = 2, 
                                  output = "label")

B.gbp_sEPD_F01 <- ggbiplot::ggbiplot(B.pca_sEPD_F01, scale = 0, var.axes = F,
                                     labels = row.names(B.clr_sEPD_F01_1))
B.gbp_sEPD_M00 <- ggbiplot::ggbiplot(B.pca_sEPD_M00, scale = 0, var.axes = F,
                                     labels = row.names(B.clr_sEPD_M00_1))
B.cord_pca_sEPD_F01 <- B.gbp_sEPD_F01[["data"]]
B.cord_pca_sEPD_M00 <- B.gbp_sEPD_M00[["data"]]
B.smp_cord_pca_sEPD_F01 <- merge(x = B.cord_pca_sEPD_F01, y = smp_dat, 
                                 all = F, sort = F, by.x = "labels", 
                                 by.y = "SampleID")
B.smp_cord_pca_sEPD_M00 <- merge(x = B.cord_pca_sEPD_M00, y = smp_dat, 
                                 all = F, sort = F, by.x = "labels", 
                                 by.y = "SampleID")

# for clustering, the process occurs as follows:
# (1) perform hierarchical clustering with multiscale bootstrap resampling
# (2) convert class "pvclust" to class "dendrogram" (requires dendextend)
# (3;4;5) create default ggraph; extract node data; merge with sample data
# NOTE: for pvclust(), input data = features as rows and samples as cols
B.pvc_sEPD_F01 <- pvclust::pvclust(B.clr_sEPD_F01_0, method.hclust = "ward.D2",
                                   method.dist = "euclidean", nboot = 10000, 
                                   parallel = T)
B.pvc_sEPD_M00 <- pvclust::pvclust(B.clr_sEPD_M00_0, method.hclust = "ward.D2",
                                   method.dist = "euclidean", nboot = 10000, 
                                   parallel = T)

B.dnd_sEPD_F01 <- as.dendrogram(B.pvc_sEPD_F01) # plot(B.dnd_sEPD_F01)
B.dnd_sEPD_M00 <- as.dendrogram(B.pvc_sEPD_M00) # plot(B.dnd_sEPD_M00)

B.grh_def_sEPD_F01 <- ggraph(B.dnd_sEPD_F01, 'dendrogram', circular = T) + 
  geom_node_text(aes(label = label, filter = leaf), hjust = 0.5, vjust = 0.5) +
  geom_node_point(aes(filter = leaf)) +
  geom_edge_elbow()
B.grh_def_sEPD_M00 <- ggraph(B.dnd_sEPD_M00, 'dendrogram', circular = T) + 
  geom_node_text(aes(label = label, filter = leaf), hjust = 0.5, vjust = 0.5) +
  geom_node_point(aes(filter = leaf)) +
  geom_edge_elbow()
B.node_dnd_sEPD_F01 <- B.grh_def_sEPD_F01[["data"]]
B.node_dnd_sEPD_M00 <- B.grh_def_sEPD_M00[["data"]]
B.smp_node_dnd_sEPD_F01 <- merge(x = B.node_dnd_sEPD_F01, y = smp_dat, 
                                 all = F, sort = F, by.x = "label", 
                                 by.y = "SampleID")
B.smp_node_dnd_sEPD_M00 <- merge(x = B.node_dnd_sEPD_M00, y = smp_dat, 
                                 all = F, sort = F, by.x = "label", 
                                 by.y = "SampleID")

### ************************************
### B - STEP  4 - define plot parameters/customize plot aesthetics ----
### ************************************

# NOTE: for consortium group aesthetics, order = Control - Rice.bran.modified

# vector for shape codes by sample Type; order = dist - cecum - prox - inoc
B.shp_typ <- c(shp_dtl, shp_cec, shp_prx, shp_ino)

# vector for medium colors by consortium group
B.sex_hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2])
B.sex_hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2])

# vector for shape colors by consortium group/sample type;
# order = C.D - C.E - C.P - C.s - R.D - R.E - R.P - R.s
B.shp_hex_F01 <- c(hex_F01_cmc[1], hex_F01_cmc[2], hex_F01_cmc[3], greydient[8], 
                   hex_F01_rmc[1], hex_F01_rmc[2], hex_F01_rmc[3], greydient[8])
B.shp_hex_M00 <- c(hex_M00_cmc[1], hex_M00_cmc[2], hex_M00_cmc[3], greydient[8], 
                   hex_M00_rmc[1], hex_M00_rmc[2], hex_M00_rmc[3], greydient[8])

# vector for shape fill by consortium group/sample type;
# order = C.D - C.E - C.P - C.s - R.D - R.E - R.P - R.s
B.shp_fil_F01 <- c(hex_F01_cmc[2], hex_F01_cmc[2], hex_F01_cmc[5], greydient[8],
                   hex_F01_rmc[2], hex_F01_rmc[2], hex_F01_rmc[5], greydient[8])
B.shp_fil_M00 <- c(hex_M00_cmc[2], hex_M00_cmc[2], hex_M00_cmc[5], greydient[8],
                   hex_M00_rmc[2], hex_M00_rmc[2], hex_M00_rmc[5], greydient[8])

# vector for shapes by murine sex
# plot order = Female mice - stool inoculum - Male mice
B.sex_shp <- c(shp_mur_fem, shp_ino, shp_mur_mal)

# plot labeling
B.ttl_gun <- "GUniFrac"
B.ttl_euc <- "Aitchison"
B.ttl_dnd <- " "

# sizing parameters
B.sze_smp_pts <- 1.01 # (pcoa/pca/dendrogram) size for sample points
B.sze_smp_txt <- 2.22 # (pcoa/pca/dendrogram) size for sample text
B.sze_int_lne <- 0.31 # (pcoa/pca) size for origin intercept lines
B.sze_int_crc <- 1.81 # (pcoa/pca) size for origin circle
B.sze_int_sqr <- 0.42 # (pcoa/pca) size for origin square
B.sze_smp_pts <- 1.01 # (pcoa/pca) size for sample points
B.sze_dnd_lne <- 0.42 # (dendrogram) size for branch lines
B.sze_ptl <- 08 # size for plot titles
B.sze_atl <- 08 # size for axis titles (x and y)
B.sze_ltx <- 5.5 # size for legend text
B.sze_lgn_pts <- 1.31 # size for legend points
B.sze_lgn_stk <- 0.55 # size for legend point stroke

# PCoA and PCA custom theme parameters
B.JacksonP <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(size = B.sze_atl),
  plot.title = element_text(hjust = 0.5, size = B.sze_ptl),
  legend.position = "none",
  panel.ontop = T,
  panel.background = element_rect(fill = NA))

# dendrogram custom theme parameters
B.TallTree <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  plot.title = element_text(hjust = 0.5, size = B.sze_ptl),
  legend.position = "none",
  panel.ontop = T,
  panel.background = element_rect(fill = NA))

# # legend custom theme parameters
# PurpleRain <- theme(
#   text = element_text(family = fnt, color = greydient[1]),
#   legend.position = c(0.5, 0.5),
#   legend.direction = "horizontal",
#   legend.background = element_rect(fill = greydient[8]),
#   legend.key = element_rect(fill = NA, color = NA),
#   legend.key.size = unit(0, "mm"),
#   legend.text = element_text(margin = margin(0, 0.95, 0, -0.84, "mm"),
#                              size = B.sze_ltx, color = greydient[1], hjust = 1))

B.PurpleRain <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  legend.position = c(0.5, 0.5),
  legend.direction = "vertical",
  legend.background = element_rect(fill = greydient[8]),
  legend.key = element_rect(fill = NA, color = NA),
  legend.box.margin = margin(1, 1, 2, 1, "mm"),
  legend.key.size = unit(0, "mm"),
  legend.text = element_text(
    margin = margin(0.5, 0, 0.5, 0, "mm"), size = B.sze_ltx, 
    color = greydient[1], hjust = 0))

### ************************************
### B - STEP  5 - format stool inoculum data for plotting ----
### ************************************

# isolate stool inoculum samples to be plotted as text rather than points
B.smp_cord_gun_sEPD_F01_s_0 <- dplyr::filter(B.smp_cord_gun_sEPD_F01,
                                             Type == "stool.inoculum")
B.smp_cord_uun_sEPD_F01_s_0 <- dplyr::filter(B.smp_cord_uun_sEPD_F01,
                                             Type == "stool.inoculum")
B.smp_cord_wun_sEPD_F01_s_0 <- dplyr::filter(B.smp_cord_wun_sEPD_F01,
                                             Type == "stool.inoculum")
B.smp_cord_gun_sEPD_M00_s_0 <- dplyr::filter(B.smp_cord_gun_sEPD_M00,
                                             Type == "stool.inoculum")
B.smp_cord_uun_sEPD_M00_s_0 <- dplyr::filter(B.smp_cord_uun_sEPD_M00,
                                             Type == "stool.inoculum")
B.smp_cord_wun_sEPD_M00_s_0 <- dplyr::filter(B.smp_cord_wun_sEPD_M00,
                                             Type == "stool.inoculum")
B.smp_cord_pca_sEPD_F01_s_0 <- dplyr::filter(B.smp_cord_pca_sEPD_F01,
                                             Type == "stool.inoculum")
B.smp_cord_pca_sEPD_M00_s_0 <- dplyr::filter(B.smp_cord_pca_sEPD_M00,
                                             Type == "stool.inoculum")
B.smp_node_dnd_sEPD_F01_s_0 <- dplyr::filter(B.smp_node_dnd_sEPD_F01,
                                             Type == "stool.inoculum")
B.smp_node_dnd_sEPD_M00_s_0 <- dplyr::filter(B.smp_node_dnd_sEPD_M00,
                                             Type == "stool.inoculum")

# define text colors for stool inoculums by consortium group
B.smp_cord_gun_sEPD_F01_s_1 <- B.smp_cord_gun_sEPD_F01_s_0
B.smp_cord_uun_sEPD_F01_s_1 <- B.smp_cord_uun_sEPD_F01_s_0
B.smp_cord_wun_sEPD_F01_s_1 <- B.smp_cord_wun_sEPD_F01_s_0
B.smp_cord_gun_sEPD_M00_s_1 <- B.smp_cord_gun_sEPD_M00_s_0
B.smp_cord_uun_sEPD_M00_s_1 <- B.smp_cord_uun_sEPD_M00_s_0
B.smp_cord_wun_sEPD_M00_s_1 <- B.smp_cord_wun_sEPD_M00_s_0
B.smp_cord_pca_sEPD_F01_s_1 <- B.smp_cord_pca_sEPD_F01_s_0
B.smp_cord_pca_sEPD_M00_s_1 <- B.smp_cord_pca_sEPD_M00_s_0
B.smp_node_dnd_sEPD_F01_s_1 <- B.smp_node_dnd_sEPD_F01_s_0
B.smp_node_dnd_sEPD_M00_s_1 <- B.smp_node_dnd_sEPD_M00_s_0

B.smp_cord_gun_sEPD_F01_s_1$color[B.smp_cord_gun_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_F01_cmc[2]
B.smp_cord_gun_sEPD_F01_s_1$color[B.smp_cord_gun_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_F01_rmc[2]
B.smp_cord_uun_sEPD_F01_s_1$color[B.smp_cord_uun_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_F01_cmc[2]
B.smp_cord_uun_sEPD_F01_s_1$color[B.smp_cord_uun_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_F01_rmc[2]
B.smp_cord_wun_sEPD_F01_s_1$color[B.smp_cord_wun_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_F01_cmc[2]
B.smp_cord_wun_sEPD_F01_s_1$color[B.smp_cord_wun_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_F01_rmc[2]
B.smp_cord_pca_sEPD_F01_s_1$color[B.smp_cord_pca_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_F01_cmc[2]
B.smp_cord_pca_sEPD_F01_s_1$color[B.smp_cord_pca_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_F01_rmc[2]
B.smp_node_dnd_sEPD_F01_s_1$color[B.smp_node_dnd_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_F01_cmc[2]
B.smp_node_dnd_sEPD_F01_s_1$color[B.smp_node_dnd_sEPD_F01_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_F01_rmc[2]
B.smp_cord_gun_sEPD_M00_s_1$color[B.smp_cord_gun_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_M00_cmc[2]
B.smp_cord_gun_sEPD_M00_s_1$color[B.smp_cord_gun_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_M00_rmc[2]
B.smp_cord_uun_sEPD_M00_s_1$color[B.smp_cord_uun_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_M00_cmc[2]
B.smp_cord_uun_sEPD_M00_s_1$color[B.smp_cord_uun_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_M00_rmc[2]
B.smp_cord_wun_sEPD_M00_s_1$color[B.smp_cord_wun_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_M00_cmc[2]
B.smp_cord_wun_sEPD_M00_s_1$color[B.smp_cord_wun_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_M00_rmc[2]
B.smp_cord_pca_sEPD_M00_s_1$color[B.smp_cord_pca_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_M00_cmc[2]
B.smp_cord_pca_sEPD_M00_s_1$color[B.smp_cord_pca_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_M00_rmc[2]
B.smp_node_dnd_sEPD_M00_s_1$color[B.smp_node_dnd_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "CMC"] <- hex_M00_cmc[2]
B.smp_node_dnd_sEPD_M00_s_1$color[B.smp_node_dnd_sEPD_M00_s_1$ConsortiumAbrv == 
                                    "RMC"] <- hex_M00_rmc[2]

B.smp_cord_gun_sEPD_F01_s <- B.smp_cord_gun_sEPD_F01_s_1
B.smp_cord_uun_sEPD_F01_s <- B.smp_cord_uun_sEPD_F01_s_1
B.smp_cord_wun_sEPD_F01_s <- B.smp_cord_wun_sEPD_F01_s_1
B.smp_cord_pca_sEPD_F01_s <- B.smp_cord_pca_sEPD_F01_s_1
B.smp_node_dnd_sEPD_F01_s <- B.smp_node_dnd_sEPD_F01_s_1
B.smp_cord_gun_sEPD_M00_s <- B.smp_cord_gun_sEPD_M00_s_1
B.smp_cord_uun_sEPD_M00_s <- B.smp_cord_uun_sEPD_M00_s_1
B.smp_cord_wun_sEPD_M00_s <- B.smp_cord_wun_sEPD_M00_s_1
B.smp_cord_pca_sEPD_M00_s <- B.smp_cord_pca_sEPD_M00_s_1
B.smp_node_dnd_sEPD_M00_s <- B.smp_node_dnd_sEPD_M00_s_1

### ************************************
### B - STEP 6a - plotting: PCoA ----
### ************************************

# NOTE: text drawn thrice to darken (avoids use of bold font face)

# pcoa with origin intercept line segments
B.gpr_gun_sEPD_F01_seg <- ggscatter(data = B.smp_cord_gun_sEPD_F01, 
                                    x = "X.Coord", y = "Y.Coord", 
                                    font.family = fnt, shape = 32, 
                                    color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = B.sze_int_sqr)

B.gpr_gun_sEPD_M00_seg <- ggscatter(data = B.smp_cord_gun_sEPD_M00, 
                                    x = "X.Coord", y = "Y.Coord", 
                                    font.family = fnt, shape = 32, 
                                    color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = B.sze_int_sqr)

# pcoa with samples shaped by type
B.gpr_gun_sEPD_F01 <- ggscatter(data = B.smp_cord_gun_sEPD_F01, 
                                x = "X.Coord", y = "Y.Coord", font.family = fnt,
                                color = "ConsortiumType", shape = "TypeLab", 
                                fill = "ConsortiumType", size = B.sze_smp_pts,
                                title = B.ttl_gun, xlab = B.pc1_gun_sEPD_F01, 
                                ylab = B.pc2_gun_sEPD_F01,
                                ggp = B.gpr_gun_sEPD_F01_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_F01_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, label = ConsortiumLab), 
            color = B.smp_cord_gun_sEPD_F01_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_F01_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, label = ConsortiumLab), 
            color = B.smp_cord_gun_sEPD_F01_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_F01_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, label = ConsortiumLab), 
            color = B.smp_cord_gun_sEPD_F01_s$color, size = B.sze_smp_txt) +
  scale_color_manual(values = B.shp_hex_F01) +
  scale_shape_manual(values = B.shp_typ) +
  scale_fill_manual(values = B.shp_fil_F01) +
  border(color = greydient[1]) +
  B.JacksonP

B.gpr_gun_sEPD_M00 <- ggscatter(data = B.smp_cord_gun_sEPD_M00, 
                                x = "X.Coord", y = "Y.Coord", font.family = fnt,
                                color = "ConsortiumType", shape = "TypeLab", 
                                fill = "ConsortiumType", size = B.sze_smp_pts,
                                title = B.ttl_gun, xlab = B.pc1_gun_sEPD_M00, 
                                ylab = B.pc2_gun_sEPD_M00,
                                ggp = B.gpr_gun_sEPD_M00_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_M00_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, label = ConsortiumLab), 
            color = B.smp_cord_gun_sEPD_M00_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_M00_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, label = ConsortiumLab), 
            color = B.smp_cord_gun_sEPD_M00_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_M00_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, label = ConsortiumLab), 
            color = B.smp_cord_gun_sEPD_M00_s$color, size = B.sze_smp_txt) +
  scale_color_manual(values = B.shp_hex_M00) +
  scale_shape_manual(values = B.shp_typ) +
  scale_fill_manual(values = B.shp_fil_M00) +
  border(color = greydient[1]) +
  B.JacksonP

# pcoa with samples shaped by murine sex
B.gpr_gun_sEPD_F01_sex <- ggscatter(data = B.smp_cord_gun_sEPD_F01, 
                                    x = "X.Coord", y = "Y.Coord", 
                                    font.family = fnt, color = "Consortium",
                                    shape = "AnimalSex", fill = greydient[8],
                                    title = B.ttl_gun, size = B.sze_smp_pts, 
                                    xlab = B.pc1_gun_sEPD_F01, 
                                    ylab = B.pc2_gun_sEPD_F01,
                                    ggp = B.gpr_gun_sEPD_F01_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_F01_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_F01_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_F01_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  scale_color_manual(values = B.sex_hex_F01) +
  scale_shape_manual(values = B.sex_shp) +
  border(color = greydient[1]) +
  B.JacksonP

B.gpr_gun_sEPD_M00_sex <- ggscatter(data = B.smp_cord_gun_sEPD_M00, 
                                    x = "X.Coord", y = "Y.Coord", 
                                    font.family = fnt, color = "Consortium",
                                    shape = "AnimalSex", fill = greydient[8],
                                    title = B.ttl_gun, size = B.sze_smp_pts,
                                    xlab = B.pc1_gun_sEPD_M00, 
                                    ylab = B.pc2_gun_sEPD_M00,
                                    ggp = B.gpr_gun_sEPD_M00_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_M00_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_M00_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_gun_sEPD_M00_s, family = fnt,
            aes(x = X.Coord, y = Y.Coord, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  scale_color_manual(values = B.sex_hex_M00) +
  scale_shape_manual(values = B.sex_shp) +
  border(color = greydient[1]) +
  B.JacksonP

### ************************************
### B - STEP 6b - plotting: PCA ----
### ************************************

# NOTE: text drawn thrice to darken (avoids use of bold font face)

# pca with origin intercept line segments
B.gpr_pca_sEPD_F01_seg <- ggscatter(data = B.smp_cord_pca_sEPD_F01, 
                                    x = "xvar", y = "yvar", font.family = fnt,
                                    shape = 32, color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = B.sze_int_sqr)

B.gpr_pca_sEPD_M00_seg <- ggscatter(data = B.smp_cord_pca_sEPD_M00, 
                                    x = "xvar", y = "yvar", font.family = fnt,
                                    shape = 32, color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = B.sze_int_sqr)

# pca with samples shaped by type
B.gpr_pca_sEPD_F01 <- ggscatter(data = B.smp_cord_pca_sEPD_F01, 
                                x = "xvar", y = "yvar", font.family = fnt,
                                color = "ConsortiumType", shape = "TypeLab", 
                                fill = "ConsortiumType", size = B.sze_smp_pts,
                                title = B.ttl_euc, xlab = B.pc1_pca_sEPD_F01, 
                                ylab = B.pc2_pca_sEPD_F01,
                                ggp = B.gpr_pca_sEPD_F01_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_F01_s, family = fnt,
            aes(x = xvar, y = yvar, label = ConsortiumLab), 
            color = B.smp_cord_pca_sEPD_F01_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_F01_s, family = fnt,
            aes(x = xvar, y = yvar, label = ConsortiumLab), 
            color = B.smp_cord_pca_sEPD_F01_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_F01_s, family = fnt,
            aes(x = xvar, y = yvar, label = ConsortiumLab), 
            color = B.smp_cord_pca_sEPD_F01_s$color, size = B.sze_smp_txt) +
  scale_color_manual(values = B.shp_hex_F01) +
  scale_shape_manual(values = B.shp_typ) +
  scale_fill_manual(values = B.shp_fil_F01) +
  border(color = greydient[1]) +
  B.JacksonP

B.gpr_pca_sEPD_M00 <- ggscatter(data = B.smp_cord_pca_sEPD_M00, 
                                x = "xvar", y = "yvar", font.family = fnt,
                                color = "ConsortiumType", shape = "TypeLab", 
                                fill = "ConsortiumType", size = B.sze_smp_pts,
                                title = B.ttl_euc, xlab = B.pc1_pca_sEPD_M00, 
                                ylab = B.pc2_pca_sEPD_M00,
                                ggp = B.gpr_pca_sEPD_M00_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_M00_s, family = fnt,
            aes(x = xvar, y = yvar, label = ConsortiumLab), 
            color = B.smp_cord_pca_sEPD_M00_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_M00_s, family = fnt,
            aes(x = xvar, y = yvar, label = ConsortiumLab), 
            color = B.smp_cord_pca_sEPD_M00_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_M00_s, family = fnt,
            aes(x = xvar, y = yvar, label = ConsortiumLab), 
            color = B.smp_cord_pca_sEPD_M00_s$color, size = B.sze_smp_txt) +
  scale_color_manual(values = B.shp_hex_M00) +
  scale_shape_manual(values = B.shp_typ) +
  scale_fill_manual(values = B.shp_fil_M00) +
  border(color = greydient[1]) +
  B.JacksonP

# pca with samples shaped by murine sex
B.gpr_pca_sEPD_F01_sex <- ggscatter(data = B.smp_cord_pca_sEPD_F01, 
                                    x = "xvar", y = "yvar", font.family = fnt,
                                    color = "Consortium", shape = "AnimalSex", 
                                    fill = greydient[8],
                                    title = B.ttl_euc, size = B.sze_smp_pts, 
                                    xlab = B.pc1_pca_sEPD_F01, 
                                    ylab = B.pc2_pca_sEPD_F01,
                                    ggp = B.gpr_pca_sEPD_F01_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_F01_s, family = fnt,
            aes(x = xvar, y = yvar, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_F01_s, family = fnt,
            aes(x = xvar, y = yvar, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_F01_s, family = fnt,
            aes(x = xvar, y = yvar, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  scale_color_manual(values = B.sex_hex_F01) +
  scale_shape_manual(values = B.sex_shp) +
  border(color = greydient[1]) +
  B.JacksonP

B.gpr_pca_sEPD_M00_sex <- ggscatter(data = B.smp_cord_pca_sEPD_M00, 
                                    x = "xvar", y = "yvar", font.family = fnt,
                                    color = "Consortium", shape = "AnimalSex", 
                                    fill = greydient[8],
                                    title = B.ttl_euc, size = B.sze_smp_pts, 
                                    xlab = B.pc1_pca_sEPD_M00, 
                                    ylab = B.pc2_pca_sEPD_M00,
                                    ggp = B.gpr_pca_sEPD_M00_seg) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_M00_s, family = fnt,
            aes(x = xvar, y = yvar, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_M00_s, family = fnt,
            aes(x = xvar, y = yvar, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_cord_pca_sEPD_M00_s, family = fnt,
            aes(x = xvar, y = yvar, color = Consortium, 
                label = ConsortiumLab), size = B.sze_smp_txt) +
  scale_color_manual(values = B.sex_hex_M00) +
  scale_shape_manual(values = B.sex_shp) +
  border(color = greydient[1]) +
  B.JacksonP

### ************************************
### B - STEP 6c - plotting: dendrograms ----
### ************************************

# NOTE: text drawn thrice to darken (avoids use of bold font face)

# circular dendrogram with samples shaped by type
B.grh_dnd_sEPD_F01 <- ggraph(B.dnd_sEPD_F01, 'dendrogram', circular = T) + 
  geom_edge_elbow(lineend = "square", linejoin = "mitre", linemitre = 2,
                  width = B.sze_dnd_lne) +
  labs(title = B.ttl_dnd) +
  geom_point(data = B.smp_node_dnd_sEPD_F01, size = B.sze_smp_pts,
             aes(x = x, y = y, color = ConsortiumType, shape = TypeLab, 
                 fill = ConsortiumType)) +
  geom_point(inherit.aes = F, data = B.smp_node_dnd_sEPD_F01_s, 
             aes(x = x, y = y), shape = 15, color = greydient[8], 
             size = B.sze_smp_pts * 1.5) +
  geom_text(inherit.aes = F, data = B.smp_node_dnd_sEPD_F01_s, family = fnt,
            aes(x = x, y = y, label = ConsortiumLab), 
            color = B.smp_node_dnd_sEPD_F01_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_node_dnd_sEPD_F01_s, family = fnt,
            aes(x = x, y = y, label = ConsortiumLab), 
            color = B.smp_node_dnd_sEPD_F01_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_node_dnd_sEPD_F01_s, family = fnt,
            aes(x = x, y = y, label = ConsortiumLab), 
            color = B.smp_node_dnd_sEPD_F01_s$color, size = B.sze_smp_txt) +
  scale_color_manual(values = B.shp_hex_F01) +
  scale_shape_manual(values = B.shp_typ) +
  scale_fill_manual(values = B.shp_fil_F01) +
  theme_pubr() +
  border(color = NULL) +
  B.TallTree

B.grh_dnd_sEPD_M00 <- ggraph(B.dnd_sEPD_M00, 'dendrogram', circular = T) + 
  geom_edge_elbow(lineend = "square", linejoin = "mitre", linemitre = 2,
                  width = B.sze_dnd_lne) +
  labs(title = B.ttl_dnd) +
  geom_point(data = B.smp_node_dnd_sEPD_M00, size = B.sze_smp_pts,
             aes(x = x, y = y, color = ConsortiumType, shape = TypeLab, 
                 fill = ConsortiumType)) +
  geom_point(inherit.aes = F, data = B.smp_node_dnd_sEPD_M00_s, 
             aes(x = x, y = y), shape = 15, color = greydient[8], 
             size = B.sze_smp_pts * 1.5) +
  geom_text(inherit.aes = F, data = B.smp_node_dnd_sEPD_M00_s, family = fnt,
            aes(x = x, y = y, label = ConsortiumLab), 
            color = B.smp_node_dnd_sEPD_M00_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_node_dnd_sEPD_M00_s, family = fnt,
            aes(x = x, y = y, label = ConsortiumLab), 
            color = B.smp_node_dnd_sEPD_M00_s$color, size = B.sze_smp_txt) +
  geom_text(inherit.aes = F, data = B.smp_node_dnd_sEPD_M00_s, family = fnt,
            aes(x = x, y = y, label = ConsortiumLab), 
            color = B.smp_node_dnd_sEPD_M00_s$color, size = B.sze_smp_txt) +
  scale_color_manual(values = B.shp_hex_M00) +
  scale_shape_manual(values = B.shp_typ) +
  scale_fill_manual(values = B.shp_fil_M00) +
  theme_pubr() +
  border(color = NULL) +
  B.TallTree

### ************************************
### B - STEP 6d - plotting: create custom legends ----
### ************************************

# consortium group-based legends (specific to donor groups)
# create data.frames with arbitrary x and y values that will be plotted ...
# ... and define legend breaks and labels

B.lgn_grp <- data.frame(group = c("A", "B", "C", "D", 
                                  "E",
                                  "F", "G", "H", "I", 
                                  "J",
                                  "K", "L", "M", "N"), x = 0, y = 0,
                        stringsAsFactors = F)

B.lgn_grp_brk <- c("A", "B", "C", "D", 
                   "E",
                   "F", "G", "H", "I", 
                   "J",
                   "K", "L", "M", "N")

B.lgn_grp_lab_F01 <- c("CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-f: inoculum",
                       "CMC-f: cecum", 
                       "CMC-f: proximal colon", 
                       "CMC-f: distal colon",
                       " ",
                       "RMC-f: inoculum",
                       "RMC-f: cecum", 
                       "RMC-f: proximal colon", 
                       "RMC-f: distal colon")

B.lgn_grp_lab_M00 <- c("CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-f: inoculum",
                       "CMC-f: cecum", 
                       "CMC-f: proximal colon", 
                       "CMC-f: distal colon",
                       " ",
                       "RMC-f: inoculum",
                       "RMC-f: cecum", 
                       "RMC-f: proximal colon", 
                       "RMC-f: distal colon")

# vector for shape codes by sample Type; 
# order = s - E - P - D for CMC then RMC (defined above)
B.lgn_shp <- c(32, 32, 32, 32,
               32,
               67, shp_cec, shp_prx, shp_dtl, 
               32, 
               82, shp_cec, shp_prx, shp_dtl)

# vector for shape colors by consortium group/sample type;
# vector for shape fill by consortium group/sample type;
B.lgn_hex_F01 <- c(greydient[8], greydient[8], greydient[8], greydient[8],
                   greydient[8],
                   hex_F01_cmc[2], hex_F01_cmc[2], hex_F01_cmc[3],
                   hex_F01_cmc[1], 
                   greydient[8],
                   hex_F01_rmc[2], hex_F01_rmc[2], hex_F01_rmc[3], 
                   hex_F01_rmc[1])
B.lgn_fil_F01 <- c(greydient[8], greydient[8], greydient[8], greydient[8],
                   greydient[8],
                   hex_F01_cmc[2], hex_F01_cmc[2], hex_F01_cmc[5],
                   hex_F01_cmc[2], 
                   greydient[8],
                   hex_F01_cmc[2], hex_F01_rmc[2], hex_F01_rmc[5], 
                   hex_F01_rmc[2])

B.lgn_hex_M00 <- c(greydient[8], greydient[8], greydient[8], greydient[8],
                   greydient[8],
                   hex_M00_cmc[2], hex_M00_cmc[2], hex_M00_cmc[3],
                   hex_M00_cmc[1], 
                   greydient[8],
                   hex_M00_rmc[2], hex_M00_rmc[2], hex_M00_rmc[3], 
                   hex_M00_rmc[1])
B.lgn_fil_M00 <- c(greydient[8], greydient[8], greydient[8], greydient[8],
                   greydient[8],
                   hex_M00_cmc[2], hex_M00_cmc[2], hex_M00_cmc[5],
                   hex_M00_cmc[2],
                   greydient[8],
                   hex_M00_cmc[2], hex_M00_rmc[2], hex_M00_rmc[5], 
                   hex_M00_rmc[2])

# (1) legend for CMC-f and RMC-f:
# (2) legend for CMC-m and RMC-m:

B.gtb_lgn_F01 <- ggplot(data = B.lgn_grp, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group, fill = group),
             size = B.sze_lgn_pts, stroke = B.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = B.lgn_shp, 
                     breaks = B.lgn_grp_brk, labels = B.lgn_grp_lab_F01) +
  scale_color_manual(name = NULL, values = B.lgn_hex_F01, 
                     breaks = B.lgn_grp_brk, labels = B.lgn_grp_lab_F01) +
  scale_fill_manual(name = NULL, values = B.lgn_fil_F01, 
                    breaks = B.lgn_grp_brk, labels = B.lgn_grp_lab_F01) +
  theme_void(base_family = fnt) +
  B.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

B.gtb_lgn_M00 <- ggplot(data = B.lgn_grp, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group, fill = group),
             size = B.sze_lgn_pts, stroke = B.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = B.lgn_shp, 
                     breaks = B.lgn_grp_brk, labels = B.lgn_grp_lab_M00) +
  scale_color_manual(name = NULL, values = B.lgn_hex_M00, 
                     breaks = B.lgn_grp_brk, labels = B.lgn_grp_lab_M00) +
  scale_fill_manual(name = NULL, values = B.lgn_fil_M00, 
                    breaks = B.lgn_grp_brk, labels = B.lgn_grp_lab_M00) +
  theme_void(base_family = fnt) +
  B.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# extract legends and convert from class 'gtable' to 'ggplot'
B.ggp_lgn_F01 <- get_legend(B.gtb_lgn_F01)
B.ggp_lgn_M00 <- get_legend(B.gtb_lgn_M00)
B.gpr_lgn_F01 <- as_ggplot(B.ggp_lgn_F01)
B.gpr_lgn_M00 <- as_ggplot(B.ggp_lgn_M00)

# murine sex-based legends
# create data.frames with arbitrary x and y values that will be plotted ...
# ... and define legend breaks and labels

B.lgn_sex_F01 <- data.frame(group = c("A", "B", "C", "D", "E",
                                      "F",
                                      "G", "H", 
                                      "I", 
                                      "J", "K", "L"), x = 0, y = 0,
                            stringsAsFactors = F)

B.lgn_sex_brk_F01 <- c("A", "B", "C", "D", "E",
                       "F",
                       "G", "H", 
                       "I", 
                       "J", "K", "L")
B.lgn_sex_lab_F01 <- c("  N: No inoculum", 
                       "CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-f: inoculum",
                       "CMC-f: Female mice",
                       " ",
                       "RMC-f: inoculum",
                       "RMC-f: Female mice", 
                       "RMC-f: Male mice")

B.lgn_sex_M00 <- data.frame(group = c("A", "B", "C", "D", "E",
                                      "F",
                                      "G", "H", "I", 
                                      "J", 
                                      "K", "L", "N"), x = 0, y = 0,
                            stringsAsFactors = F)

B.lgn_sex_brk_M00 <- c("A", "B", "C", "D", "E",
                       "F",
                       "G", "H", "I", 
                       "J", 
                       "K", "L", "N")
B.lgn_sex_lab_M00 <- c("  N: No inoculum", 
                       "CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-f: inoculum",
                       "CMC-f: Female mice", 
                       "CMC-f: Male mice", 
                       " ",
                       "RMC-f: inoculum",
                       "RMC-f: Female mice",
                       "RMC-f: Male mice")

# vector for colors by consortium group; vectors for shapes by murine sex
# order for F01 = C.s - C.Female - e - R.s - R.Female - R.male
# order for M00 = C.s - C.Female - C.Male - e - R.s - R.Female - R.male
B.sex_hex_F01 <- c(rep(greydient[8], times = 6),
                   rep(hex_F01_cmc[2], times = 2), 
                   greydient[8],
                   rep(hex_F01_rmc[2], times = 3))
B.sex_hex_M00 <- c(rep(greydient[8], times = 6),
                   rep(hex_M00_cmc[2], times = 3), 
                   greydient[8],
                   rep(hex_M00_rmc[2], times = 3))
B.sex_shp_F01 <- c(rep(32, times = 6),
                   67, shp_mur_fem, 
                   32, 
                   82, shp_mur_fem, shp_mur_mal)
B.sex_shp_M00 <- c(rep(32, times = 6),
                   67, shp_mur_fem, shp_mur_mal, 
                   32, 
                   82, shp_mur_fem, shp_mur_mal)

# (1) legend for CMC-f and RMC-f:
# (2) legend for CMC-m and RMC-m:

B.gtb_lgn_sex_F01 <- ggplot(data = B.lgn_sex_F01, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group), fill = greydient[8],
             size = B.sze_lgn_pts, stroke = B.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = B.sex_shp_F01, 
                     breaks = B.lgn_sex_brk_F01, labels = B.lgn_sex_lab_F01) +
  scale_color_manual(name = NULL, values = B.sex_hex_F01, 
                     breaks = B.lgn_sex_brk_F01, labels = B.lgn_sex_lab_F01) +
  theme_void(base_family = fnt) +
  B.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

B.gtb_lgn_sex_M00 <- ggplot(data = B.lgn_sex_M00, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group), fill = greydient[8],
             size = B.sze_lgn_pts, stroke = B.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = B.sex_shp_M00, 
                     breaks = B.lgn_sex_brk_M00, labels = B.lgn_sex_lab_M00) +
  scale_color_manual(name = NULL, values = B.sex_hex_M00, 
                     breaks = B.lgn_sex_brk_M00, labels = B.lgn_sex_lab_M00) +
  theme_void(base_family = fnt) +
  B.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# extract legends and convert from class 'gtable' to 'ggplot'
B.ggp_lgn_sex_F01 <- get_legend(B.gtb_lgn_sex_F01)
B.ggp_lgn_sex_M00 <- get_legend(B.gtb_lgn_sex_M00)
B.gpr_lgn_sex_F01 <- as_ggplot(B.ggp_lgn_sex_F01)
B.gpr_lgn_sex_M00 <- as_ggplot(B.ggp_lgn_sex_M00)

### ************************************
### B - STEP 6e - plotting: arrange ----
### ************************************

# arrange plots by sex of human donor
B.gga_plt_F01 <- ggarrange(B.gpr_gun_sEPD_F01, B.gpr_pca_sEPD_F01,
                           B.grh_dnd_sEPD_F01, 
                           labels = c("A", "B", ""), font.label = pan_fnt, 
                           ncol = 3, nrow = 1, align = "hv")
B.gga_plt_M00 <- ggarrange(B.gpr_gun_sEPD_M00, B.gpr_pca_sEPD_M00,
                           B.grh_dnd_sEPD_M00, 
                           labels = c("C", "D", ""), font.label = pan_fnt, 
                           ncol = 3, nrow = 1, align = "hv")

# arrange plots and legends by sex of human donor
B.gga_plt_lgn_F01 <- ggarrange(B.gga_plt_F01, B.gpr_lgn_F01, labels = NULL,
                               ncol = 2, nrow = 1, widths = c(3, 0.95))
B.gga_plt_lgn_M00 <- ggarrange(B.gga_plt_M00, B.gpr_lgn_M00, labels = NULL,
                               ncol = 2, nrow = 1, widths = c(3, 0.95))

# arrange arrangements
B.gga_sEPD <- ggarrange(B.gga_plt_lgn_F01, B.gga_plt_lgn_M00, labels = NULL,
                        ncol = 1, nrow = 2, align = "hv")

# murine sex-based plots - arrange plots by sex of human donor
B.gga_plt_sex_F01 <- ggarrange(B.gpr_gun_sEPD_F01_sex, B.gpr_pca_sEPD_F01_sex,
                               labels = c("B", ""), font.label = pan_fnt, 
                               ncol = 2, nrow = 1, align = "hv")
B.gga_plt_sex_M00 <- ggarrange(B.gpr_gun_sEPD_M00_sex, B.gpr_pca_sEPD_M00_sex,
                               labels = c("E", ""), font.label = pan_fnt, 
                               ncol = 2, nrow = 1, align = "hv")

# murine sex-based plots - arrange plots and legends by sex of human donor
B.gga_plt_sex_lgn_F01 <- ggarrange(B.gga_plt_sex_F01, B.gpr_lgn_sex_F01, 
                                   labels = NULL, ncol = 2, nrow = 1, 
                                   widths = c(2, 0.95))
B.gga_plt_sex_lgn_M00 <- ggarrange(B.gga_plt_sex_M00, B.ggp_lgn_sex_M00, 
                                   labels = NULL, ncol = 2, nrow = 1, 
                                   widths = c(2, 0.95))

### ************************************
### ^^^B - STEP  7 - perform pairwise PERMANOVA ----
### ************************************

# ^^^ = code is mostly not commented in this step
# consider breaking this up into two steps... one for perform, one for format

# for CZM/clr dfs, the process occurs as follows:
# (1) transpose to convert features to cols and samples to rows
# (2) create distance matrices using the euclidean metric
# NOTE: for dist(), input data = features as cols and samples as rows
B.clr_EPD_F01_1 <- as.data.frame(t(B.clr_EPD_F01_0))
B.clr_s_F01_1 <- as.data.frame(t(B.clr_s_F01_0))
B.clr_EPD_M00_1 <- as.data.frame(t(B.clr_EPD_M00_0))
B.clr_s_M00_1 <- as.data.frame(t(B.clr_s_M00_0))

B.euc_EPD_F01 <- dist(B.clr_EPD_F01_1, method = "euclidean")
B.euc_s_F01 <- dist(B.clr_s_F01_1, method = "euclidean")
B.euc_EPD_M00 <- dist(B.clr_EPD_M00_1, method = "euclidean")
B.euc_s_M00 <- dist(B.clr_s_M00_1, method = "euclidean")

# filter sample data & distance matrices to have identical samples
B.list_gun_EPD_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.gun_EPD_F01, 
                                  flt_col = "HumanDonorSex", flt_row = "Female", 
                                  match_col = "SampleID")
B.list_uun_EPD_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.uun_EPD_F01, 
                                  flt_col = "HumanDonorSex", flt_row = "Female", 
                                  match_col = "SampleID")
B.list_wun_EPD_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.wun_EPD_F01, 
                                  flt_col = "HumanDonorSex", flt_row = "Female", 
                                  match_col = "SampleID")
B.list_euc_EPD_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.euc_EPD_F01, 
                                  flt_col = "HumanDonorSex", flt_row = "Female", 
                                  match_col = "SampleID")
B.list_gun_s_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.gun_s_F01, 
                                flt_col = "HumanDonorSex", flt_row = "Female", 
                                match_col = "SampleID")
B.list_uun_s_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.uun_s_F01, 
                                flt_col = "HumanDonorSex", flt_row = "Female", 
                                match_col = "SampleID")
B.list_wun_s_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.wun_s_F01, 
                                flt_col = "HumanDonorSex", flt_row = "Female", 
                                match_col = "SampleID")
B.list_euc_s_F01 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.euc_s_F01, 
                                flt_col = "HumanDonorSex", flt_row = "Female", 
                                match_col = "SampleID")
B.list_gun_EPD_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.gun_EPD_M00, 
                                  flt_col = "HumanDonorSex", flt_row = "Male", 
                                  match_col = "SampleID")
B.list_uun_EPD_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.uun_EPD_M00, 
                                  flt_col = "HumanDonorSex", flt_row = "Male", 
                                  match_col = "SampleID")
B.list_wun_EPD_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.wun_EPD_M00, 
                                  flt_col = "HumanDonorSex", flt_row = "Male", 
                                  match_col = "SampleID")
B.list_euc_EPD_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.euc_EPD_M00, 
                                  flt_col = "HumanDonorSex", flt_row = "Male", 
                                  match_col = "SampleID")
B.list_gun_s_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.gun_s_M00, 
                                flt_col = "HumanDonorSex", flt_row = "Male", 
                                match_col = "SampleID")
B.list_uun_s_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.uun_s_M00, 
                                flt_col = "HumanDonorSex", flt_row = "Male", 
                                match_col = "SampleID")
B.list_wun_s_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.wun_s_M00, 
                                flt_col = "HumanDonorSex", flt_row = "Male", 
                                match_col = "SampleID")
B.list_euc_s_M00 <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.euc_s_M00, 
                                flt_col = "HumanDonorSex", flt_row = "Male", 
                                match_col = "SampleID")

B.smp_gun_EPD_F01 <- B.list_gun_EPD_F01$new_smp
B.dst_gun_EPD_F01 <- B.list_gun_EPD_F01$new_dst
B.smp_uun_EPD_F01 <- B.list_uun_EPD_F01$new_smp
B.dst_uun_EPD_F01 <- B.list_uun_EPD_F01$new_dst
B.smp_wun_EPD_F01 <- B.list_wun_EPD_F01$new_smp
B.dst_wun_EPD_F01 <- B.list_wun_EPD_F01$new_dst
B.smp_euc_EPD_F01 <- B.list_euc_EPD_F01$new_smp
B.dst_euc_EPD_F01 <- B.list_euc_EPD_F01$new_dst
B.smp_gun_s_F01 <- B.list_gun_s_F01$new_smp
B.dst_gun_s_F01 <- B.list_gun_s_F01$new_dst
B.smp_uun_s_F01 <- B.list_uun_s_F01$new_smp
B.dst_uun_s_F01 <- B.list_uun_s_F01$new_dst
B.smp_wun_s_F01 <- B.list_wun_s_F01$new_smp
B.dst_wun_s_F01 <- B.list_wun_s_F01$new_dst
B.smp_euc_s_F01 <- B.list_euc_s_F01$new_smp
B.dst_euc_s_F01 <- B.list_euc_s_F01$new_dst
B.smp_gun_EPD_M00 <- B.list_gun_EPD_M00$new_smp
B.dst_gun_EPD_M00 <- B.list_gun_EPD_M00$new_dst
B.smp_uun_EPD_M00 <- B.list_uun_EPD_M00$new_smp
B.dst_uun_EPD_M00 <- B.list_uun_EPD_M00$new_dst
B.smp_wun_EPD_M00 <- B.list_wun_EPD_M00$new_smp
B.dst_wun_EPD_M00 <- B.list_wun_EPD_M00$new_dst
B.smp_euc_EPD_M00 <- B.list_euc_EPD_M00$new_smp
B.dst_euc_EPD_M00 <- B.list_euc_EPD_M00$new_dst
B.smp_gun_s_M00 <- B.list_gun_s_M00$new_smp
B.dst_gun_s_M00 <- B.list_gun_s_M00$new_dst
B.smp_uun_s_M00 <- B.list_uun_s_M00$new_smp
B.dst_uun_s_M00 <- B.list_uun_s_M00$new_dst
B.smp_wun_s_M00 <- B.list_wun_s_M00$new_smp
B.dst_wun_s_M00 <- B.list_wun_s_M00$new_dst
B.smp_euc_s_M00 <- B.list_euc_s_M00$new_smp
B.dst_euc_s_M00 <- B.list_euc_s_M00$new_dst

# HERE for global adonis code ----
# perform PERMANOVA using adonis on the _s dfs (human stool inoculums)

# perfrom adonis on the filtered data
B.rst_gbl_gun_s_F01_0 <- vegan::adonis(
  formula = B.dst_gun_s_F01 ~ B.smp_gun_s_F01[, "ConsortiumType"],
  data = B.smp_gun_s_F01, permutations = 9999)
B.rst_gbl_uun_s_F01_0 <- vegan::adonis(
  formula = B.dst_uun_s_F01 ~ B.smp_uun_s_F01[, "ConsortiumType"],
  data = B.smp_uun_s_F01, permutations = 9999)
B.rst_gbl_wun_s_F01_0 <- vegan::adonis(
  formula = B.dst_wun_s_F01 ~ B.smp_wun_s_F01[, "ConsortiumType"],
  data = B.smp_wun_s_F01, permutations = 9999)
B.rst_gbl_euc_s_F01_0 <- vegan::adonis(
  formula = B.dst_euc_s_F01 ~ B.smp_euc_s_F01[, "ConsortiumType"],
  data = B.smp_euc_s_F01, permutations = 9999)

B.rst_gbl_gun_s_M00_0 <- vegan::adonis(
  formula = B.dst_gun_s_M00 ~ B.smp_gun_s_M00[, "ConsortiumType"],
  data = B.smp_gun_s_M00, permutations = 9999)
B.rst_gbl_uun_s_M00_0 <- vegan::adonis(
  formula = B.dst_uun_s_M00 ~ B.smp_uun_s_M00[, "ConsortiumType"],
  data = B.smp_uun_s_M00, permutations = 9999)
B.rst_gbl_wun_s_M00_0 <- vegan::adonis(
  formula = B.dst_wun_s_M00 ~ B.smp_wun_s_M00[, "ConsortiumType"],
  data = B.smp_wun_s_M00, permutations = 9999)
B.rst_gbl_euc_s_M00_0 <- vegan::adonis(
  formula = B.dst_euc_s_M00 ~ B.smp_euc_s_M00[, "ConsortiumType"],
  data = B.smp_euc_s_M00, permutations = 9999)

B.rst_gbl_gun_s_F01_0
B.rst_gbl_uun_s_F01_0
B.rst_gbl_wun_s_F01_0
B.rst_gbl_euc_s_F01_0

B.rst_gbl_gun_s_M00_0
B.rst_gbl_uun_s_M00_0
B.rst_gbl_wun_s_M00_0
B.rst_gbl_euc_s_M00_0

# flt_dst <- filt$new_dst
# # perform adonis on the filtered data
# adon <- vegan::adonis(formula = flt_dst ~ flt_smp[, test_col], 
#                       data = flt_smp, permutations = perms)
# # define labels for the individual pair
# pair_labs <- c(pair_labs, paste(pair[1], 'vs', pair[2]))
# # store information from the adonis data.frame
# R2 <- c(R2, adon$aov.tab[1, 5])
# p_raw <- c(p_raw, adon$aov.tab[1, 6])

# B.rst_prm_gun_EPD_F01 <- pwise_adon(smp_dat = B.smp_gun_EPD_F01, 
#                                     dst_obj = B.dst_gun_EPD_F01,
#                                     test_col = "ConsortiumType",
#                                     match_col = "SampleID", perms = 9999, 
#                                     adjust = "BH", digits = 3)

# code resumes ----

# perform pairwise PERMANOVA using adonis (relevant for _EPD dfs only)
B.rst_prm_gun_EPD_F01 <- pwise_adon(smp_dat = B.smp_gun_EPD_F01, 
                                    dst_obj = B.dst_gun_EPD_F01,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)
B.rst_prm_uun_EPD_F01 <- pwise_adon(smp_dat = B.smp_uun_EPD_F01,
                                    dst_obj = B.dst_uun_EPD_F01,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)
B.rst_prm_wun_EPD_F01 <- pwise_adon(smp_dat = B.smp_wun_EPD_F01,
                                    dst_obj = B.dst_wun_EPD_F01,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)
B.rst_prm_euc_EPD_F01 <- pwise_adon(smp_dat = B.smp_euc_EPD_F01, 
                                    dst_obj = B.dst_euc_EPD_F01,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)
B.rst_prm_gun_EPD_M00 <- pwise_adon(smp_dat = B.smp_gun_EPD_M00, 
                                    dst_obj = B.dst_gun_EPD_M00,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)
B.rst_prm_uun_EPD_M00 <- pwise_adon(smp_dat = B.smp_uun_EPD_M00,
                                    dst_obj = B.dst_uun_EPD_M00,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)
B.rst_prm_wun_EPD_M00 <- pwise_adon(smp_dat = B.smp_wun_EPD_M00,
                                    dst_obj = B.dst_wun_EPD_M00,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)
B.rst_prm_euc_EPD_M00 <- pwise_adon(smp_dat = B.smp_euc_EPD_M00, 
                                    dst_obj = B.dst_euc_EPD_M00,
                                    test_col = "ConsortiumType",
                                    match_col = "SampleID", perms = 9999, 
                                    adjust = "BH", digits = 3)

# extract results from between groups (type matched)
B.rst_prm_gun_EPD_F01_btw_0 <- dplyr::filter(B.rst_prm_gun_EPD_F01,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_uun_EPD_F01_btw_0 <- dplyr::filter(B.rst_prm_uun_EPD_F01,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_wun_EPD_F01_btw_0 <- dplyr::filter(B.rst_prm_wun_EPD_F01,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_euc_EPD_F01_btw_0 <- dplyr::filter(B.rst_prm_euc_EPD_F01,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_gun_EPD_M00_btw_0 <- dplyr::filter(B.rst_prm_gun_EPD_M00,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_uun_EPD_M00_btw_0 <- dplyr::filter(B.rst_prm_uun_EPD_M00,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_wun_EPD_M00_btw_0 <- dplyr::filter(B.rst_prm_wun_EPD_M00,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_euc_EPD_M00_btw_0 <- dplyr::filter(B.rst_prm_euc_EPD_M00,
                                             Comparison == "C.E vs R.E" |
                                               Comparison == "C.P vs R.P" |
                                               Comparison == "C.D vs R.D")
B.rst_prm_gun_EPD_F01_wtn_cmc_0 <- dplyr::filter(B.rst_prm_gun_EPD_F01,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_uun_EPD_F01_wtn_cmc_0 <- dplyr::filter(B.rst_prm_uun_EPD_F01,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_wun_EPD_F01_wtn_cmc_0 <- dplyr::filter(B.rst_prm_wun_EPD_F01,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_euc_EPD_F01_wtn_cmc_0 <- dplyr::filter(B.rst_prm_euc_EPD_F01,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_gun_EPD_M00_wtn_cmc_0 <- dplyr::filter(B.rst_prm_gun_EPD_M00,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_uun_EPD_M00_wtn_cmc_0 <- dplyr::filter(B.rst_prm_uun_EPD_M00,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_wun_EPD_M00_wtn_cmc_0 <- dplyr::filter(B.rst_prm_wun_EPD_M00,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_euc_EPD_M00_wtn_cmc_0 <- dplyr::filter(B.rst_prm_euc_EPD_M00,
                                                 Comparison ==  "C.E vs C.P" |
                                                   Comparison ==  "C.E vs C.D" |
                                                   Comparison ==  "C.P vs C.D")
B.rst_prm_gun_EPD_F01_wtn_rmc_0 <- dplyr::filter(B.rst_prm_gun_EPD_F01,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")
B.rst_prm_uun_EPD_F01_wtn_rmc_0 <- dplyr::filter(B.rst_prm_uun_EPD_F01,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")
B.rst_prm_wun_EPD_F01_wtn_rmc_0 <- dplyr::filter(B.rst_prm_wun_EPD_F01,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")
B.rst_prm_euc_EPD_F01_wtn_rmc_0 <- dplyr::filter(B.rst_prm_euc_EPD_F01,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")
B.rst_prm_gun_EPD_M00_wtn_rmc_0 <- dplyr::filter(B.rst_prm_gun_EPD_M00,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")
B.rst_prm_uun_EPD_M00_wtn_rmc_0 <- dplyr::filter(B.rst_prm_uun_EPD_M00,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")
B.rst_prm_wun_EPD_M00_wtn_rmc_0 <- dplyr::filter(B.rst_prm_wun_EPD_M00,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")
B.rst_prm_euc_EPD_M00_wtn_rmc_0 <- dplyr::filter(B.rst_prm_euc_EPD_M00,
                                                 Comparison == "R.E vs R.D" |
                                                   Comparison == "R.E vs R.P" |
                                                   Comparison ==  "R.D vs R.P")

B.rst_prm_gun_EPD_F01_btw_1 <- B.rst_prm_gun_EPD_F01_btw_0
B.rst_prm_uun_EPD_F01_btw_1 <- B.rst_prm_uun_EPD_F01_btw_0
B.rst_prm_wun_EPD_F01_btw_1 <- B.rst_prm_wun_EPD_F01_btw_0
B.rst_prm_euc_EPD_F01_btw_1 <- B.rst_prm_euc_EPD_F01_btw_0
B.rst_prm_gun_EPD_M00_btw_1 <- B.rst_prm_gun_EPD_M00_btw_0
B.rst_prm_uun_EPD_M00_btw_1 <- B.rst_prm_uun_EPD_M00_btw_0
B.rst_prm_wun_EPD_M00_btw_1 <- B.rst_prm_wun_EPD_M00_btw_0
B.rst_prm_euc_EPD_M00_btw_1 <- B.rst_prm_euc_EPD_M00_btw_0




B.rst_prm_gun_EPD_F01_wtn_cmc_1 <- B.rst_prm_gun_EPD_F01_wtn_cmc_0
B.rst_prm_uun_EPD_F01_wtn_cmc_1 <- B.rst_prm_uun_EPD_F01_wtn_cmc_0
B.rst_prm_wun_EPD_F01_wtn_cmc_1 <- B.rst_prm_wun_EPD_F01_wtn_cmc_0
B.rst_prm_euc_EPD_F01_wtn_cmc_1 <- B.rst_prm_euc_EPD_F01_wtn_cmc_0
B.rst_prm_gun_EPD_M00_wtn_cmc_1 <- B.rst_prm_gun_EPD_M00_wtn_cmc_0
B.rst_prm_uun_EPD_M00_wtn_cmc_1 <- B.rst_prm_uun_EPD_M00_wtn_cmc_0
B.rst_prm_wun_EPD_M00_wtn_cmc_1 <- B.rst_prm_wun_EPD_M00_wtn_cmc_0
B.rst_prm_euc_EPD_M00_wtn_cmc_1 <- B.rst_prm_euc_EPD_M00_wtn_rmc_0
B.rst_prm_gun_EPD_F01_wtn_rmc_1 <- B.rst_prm_gun_EPD_F01_wtn_rmc_0
B.rst_prm_uun_EPD_F01_wtn_rmc_1 <- B.rst_prm_uun_EPD_F01_wtn_rmc_0
B.rst_prm_wun_EPD_F01_wtn_rmc_1 <- B.rst_prm_wun_EPD_F01_wtn_rmc_0
B.rst_prm_euc_EPD_F01_wtn_rmc_1 <- B.rst_prm_euc_EPD_F01_wtn_rmc_0
B.rst_prm_gun_EPD_M00_wtn_rmc_1 <- B.rst_prm_gun_EPD_M00_wtn_rmc_0
B.rst_prm_uun_EPD_M00_wtn_rmc_1 <- B.rst_prm_uun_EPD_M00_wtn_rmc_0
B.rst_prm_wun_EPD_M00_wtn_rmc_1 <- B.rst_prm_wun_EPD_M00_wtn_rmc_0
B.rst_prm_euc_EPD_M00_wtn_rmc_1 <- B.rst_prm_euc_EPD_M00_wtn_rmc_0

B.comp_btw <- "Control Consortium vs. Rice Bran Modified Consortium"
B.comp_wtn_cmc <- "within Control Consortium"
B.comp_wtn_rmc <- "within Rice Bran Modified Consortium"
B.rst_prm_gun_EPD_F01_btw_1$comp <- B.comp_btw
B.rst_prm_uun_EPD_F01_btw_1$comp <- B.comp_btw
B.rst_prm_wun_EPD_F01_btw_1$comp <- B.comp_btw
B.rst_prm_euc_EPD_F01_btw_1$comp <- B.comp_btw
B.rst_prm_gun_EPD_M00_btw_1$comp <- B.comp_btw
B.rst_prm_uun_EPD_M00_btw_1$comp <- B.comp_btw
B.rst_prm_wun_EPD_M00_btw_1$comp <- B.comp_btw
B.rst_prm_euc_EPD_M00_btw_1$comp <- B.comp_btw
B.rst_prm_gun_EPD_F01_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_uun_EPD_F01_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_wun_EPD_F01_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_euc_EPD_F01_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_gun_EPD_M00_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_uun_EPD_M00_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_wun_EPD_M00_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_euc_EPD_M00_wtn_cmc_1$comp <- B.comp_wtn_cmc
B.rst_prm_gun_EPD_F01_wtn_rmc_1$comp <- B.comp_wtn_rmc
B.rst_prm_uun_EPD_F01_wtn_rmc_1$comp <- B.comp_wtn_rmc
B.rst_prm_wun_EPD_F01_wtn_rmc_1$comp <- B.comp_wtn_rmc
B.rst_prm_euc_EPD_F01_wtn_rmc_1$comp <- B.comp_wtn_rmc
B.rst_prm_gun_EPD_M00_wtn_rmc_1$comp <- B.comp_wtn_rmc
B.rst_prm_uun_EPD_M00_wtn_rmc_1$comp <- B.comp_wtn_rmc
B.rst_prm_wun_EPD_M00_wtn_rmc_1$comp <- B.comp_wtn_rmc
B.rst_prm_euc_EPD_M00_wtn_rmc_1$comp <- B.comp_wtn_rmc

B.stats_gun_EPD_F01_0 <- rbind(B.rst_prm_gun_EPD_F01_btw_1, 
                               B.rst_prm_gun_EPD_F01_wtn_cmc_1, 
                               B.rst_prm_gun_EPD_F01_wtn_rmc_1)
B.stats_uun_EPD_F01_0 <- rbind(B.rst_prm_uun_EPD_F01_btw_1, 
                               B.rst_prm_uun_EPD_F01_wtn_cmc_1, 
                               B.rst_prm_uun_EPD_F01_wtn_rmc_1)
B.stats_wun_EPD_F01_0 <- rbind(B.rst_prm_wun_EPD_F01_btw_1, 
                               B.rst_prm_wun_EPD_F01_wtn_cmc_1, 
                               B.rst_prm_wun_EPD_F01_wtn_rmc_1)
B.stats_euc_EPD_F01_0 <- rbind(B.rst_prm_euc_EPD_F01_btw_1, 
                               B.rst_prm_euc_EPD_F01_wtn_cmc_1, 
                               B.rst_prm_euc_EPD_F01_wtn_rmc_1)
B.stats_gun_EPD_M00_0 <- rbind(B.rst_prm_gun_EPD_M00_btw_1, 
                               B.rst_prm_gun_EPD_M00_wtn_cmc_1, 
                               B.rst_prm_gun_EPD_M00_wtn_rmc_1)
B.stats_uun_EPD_M00_0 <- rbind(B.rst_prm_uun_EPD_M00_btw_1, 
                               B.rst_prm_uun_EPD_M00_wtn_cmc_1, 
                               B.rst_prm_uun_EPD_M00_wtn_rmc_1)
B.stats_wun_EPD_M00_0 <- rbind(B.rst_prm_wun_EPD_M00_btw_1, 
                               B.rst_prm_wun_EPD_M00_wtn_cmc_1, 
                               B.rst_prm_wun_EPD_M00_wtn_rmc_1)
B.stats_euc_EPD_M00_0 <- rbind(B.rst_prm_euc_EPD_M00_btw_1, 
                               B.rst_prm_euc_EPD_M00_wtn_cmc_1, 
                               B.rst_prm_euc_EPD_M00_wtn_rmc_1)

B.stats_gun_EPD_F01_1 <- B.stats_gun_EPD_F01_0
B.stats_uun_EPD_F01_1 <- B.stats_uun_EPD_F01_0 
B.stats_wun_EPD_F01_1 <- B.stats_wun_EPD_F01_0
B.stats_euc_EPD_F01_1 <- B.stats_euc_EPD_F01_0
B.stats_gun_EPD_M00_1 <- B.stats_gun_EPD_M00_0
B.stats_uun_EPD_M00_1 <- B.stats_uun_EPD_M00_0 
B.stats_wun_EPD_M00_1 <- B.stats_wun_EPD_M00_0
B.stats_euc_EPD_M00_1 <- B.stats_euc_EPD_M00_0

B.stats_gun_EPD_F01_1$metric <- "generalized UniFrac"
B.stats_uun_EPD_F01_1$metric <- "unweighted UniFrac"
B.stats_wun_EPD_F01_1$metric <- "weighted UniFrac"
B.stats_euc_EPD_F01_1$metric <- "euclidean/Aitchison"
B.stats_gun_EPD_M00_1$metric <- "generalized UniFrac"
B.stats_uun_EPD_M00_1$metric <- "unweighted UniFrac"
B.stats_wun_EPD_M00_1$metric <- "weighted UniFrac"
B.stats_euc_EPD_M00_1$metric <- "euclidean/Aitchison"

B.stats_gun_EPD_F01_1$group <- "consortia from human female"
B.stats_uun_EPD_F01_1$group <- "consortia from human female"
B.stats_wun_EPD_F01_1$group <- "consortia from human female"
B.stats_euc_EPD_F01_1$group <- "consortia from human female"
B.stats_gun_EPD_M00_1$group <- "consortia from human male"
B.stats_uun_EPD_M00_1$group <- "consortia from human male"
B.stats_wun_EPD_M00_1$group <- "consortia from human male"
B.stats_euc_EPD_M00_1$group <- "consortia from human male"

B.stats_EPD_F01_0 <- rbind(B.stats_gun_EPD_F01_1, B.stats_uun_EPD_F01_1,
                           B.stats_wun_EPD_F01_1, B.stats_euc_EPD_F01_1)
B.stats_EPD_M00_0 <- rbind(B.stats_gun_EPD_M00_1, B.stats_uun_EPD_M00_1,
                           B.stats_wun_EPD_M00_1, B.stats_euc_EPD_M00_1)

B.stats_EPD_F01_1 <- dplyr::rename(B.stats_EPD_F01_0, type = Comparison,
                                   comparison = comp)
B.stats_EPD_M00_1 <- dplyr::rename(B.stats_EPD_M00_0, type = Comparison,
                                   comparison = comp)

B.stats_EPD_F01_2 <- B.stats_EPD_F01_1
B.stats_EPD_M00_2 <- B.stats_EPD_M00_1

B.stats_EPD_F01_2$type <- gsub(pattern = "C.E vs R.E", 
                               replacement = "cecum", 
                               x = B.stats_EPD_F01_2$type)
B.stats_EPD_F01_2$type <- gsub(pattern = "C.P vs R.P", 
                               replacement = "colon-prox", 
                               x = B.stats_EPD_F01_2$type)
B.stats_EPD_F01_2$type <- gsub(pattern = "C.D vs R.D", 
                               replacement = "colon-dist", 
                               x = B.stats_EPD_F01_2$type)
B.stats_EPD_M00_2$type <- gsub(pattern = "C.E vs R.E", 
                               replacement = "cecum", 
                               x = B.stats_EPD_M00_2$type)
B.stats_EPD_M00_2$type <- gsub(pattern = "C.P vs R.P", 
                               replacement = "colon-prox", 
                               x = B.stats_EPD_M00_2$type)
B.stats_EPD_M00_2$type <- gsub(pattern = "C.D vs R.D", 
                               replacement = "colon-dist", 
                               x = B.stats_EPD_M00_2$type)

B.stats_EPD_0 <- rbind(B.stats_EPD_F01_2, B.stats_EPD_M00_2)

# format P and adjusted p values to reflect writing in manuscript
B.stats_EPD <- dplyr::rename(B.stats_EPD_0, P = p_raw, BH_P = p_BH)

### ************************************
### B - WRITE OUTPUTS ----
### ************************************

# consortium-based plots
B.ofv_gga_plt_lgn_F01 <- "TransFaunation/vault/plot_pcoa_pca_dnd_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 45,
       filename = B.ofv_gga_plt_lgn_F01, plot = B.gga_plt_lgn_F01)

B.ofv_gga_plt_lgn_M00 <- "TransFaunation/vault/plot_pcoa_pca_dnd_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 45,
       filename = B.ofv_gga_plt_lgn_M00, plot = B.gga_plt_lgn_M00)

B.ofv_gga_sEPD <- "TransFaunation/vault/plot_pcoa_pca_dnd.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 90,
       filename = B.ofv_gga_sEPD, plot = B.gga_sEPD)

# murine sex-based plots
B.ofv_gga_F01_sex <- "TransFaunation/vault/plot_pcoa_pca_F01_animal_sex.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 125, height = 45,
       filename = B.ofv_gga_F01_sex, plot = B.gga_plt_sex_lgn_F01)

B.ofv_gga_M00_sex <- "TransFaunation/vault/plot_pcoa_pca_M00_animal_sex.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 125, height = 45,
       filename = B.ofv_gga_M00_sex, plot = B.gga_plt_sex_lgn_M00)

# stats
B.ofv_stats_EPD <- "TransFaunation/vault/table_stats_beta_div.txt"
write.table(sep = "\t", row.names = F, x = B.stats_EPD, file = B.ofv_stats_EPD)

# save workspace
B.obj <- ls(pattern = "B.")
B.lst <- c(B.obj[grep(pattern = "B.", x = B.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, B.obj_from_ST, COSMOS)
save(list = B.lst, file = B.ofv_wksp)

### ************************************
### BEGIN Section C ----
### ************************************

# note for kdp: taken from TransFaunation/section_C.R ; section_C_taxa_in_CRC.R

# NOTE: pay attention to the comments as some may not reflect this code ...
# ... but rather were taken from iBright_R/aldex2/aldex2.R

#### preface ----

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)
require(ggraph, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
load(file = A.ofv_wksp)
load(file = ofv_COSMOS_wksp)

C.ofv_wksp <- paste(sep = "", path_vault, "/WS_", "SectionC.RData")

### ************************************
### ^^^C - STEP  1 - compute CZM & collapse taxonomy ----
### ************************************

# provide provenance for information gathering at end of section:
C.prov_secstep_CS1 <- "Section C - STEP 1"
C.prov_heading_CS1 <- "compute CZM & collapse taxonomy"
# ^^^C.prov_output_obj_CS1 <- "" # this object is output to the vault'
# ^^^info and prov have not actually been added to this script

# ^^^NOTE: if the environment is empty; there are some requirements:

# create a new version of the objects needed from sections A ...
# ... and create a vector naming those objects (used when saving the workspace)
C.seq_tax <- A.seq_tax
C.EBTKS_abs_pro_0 <- A.EBTKS_abs_pro
C.obj_from_A <- c("A.seq_tax", "A.EBTKS_abs_pro")

# create data.frames to store SILVA lineage and relevant assignment levels info
# (1) retain relevant columns (renaming them in the process)
# (2) reduce data.frames to retain single assignments

C.slv_Phy_0 <- dplyr::select(C.EBTKS_abs_pro_0, taxon = int.slv.Phylum)
C.slv_Fam_0 <- dplyr::select(C.EBTKS_abs_pro_0, taxon = int.slv.Family)
C.slv_Lws_0 <- dplyr::select(C.EBTKS_abs_pro_0, taxon = int.slv.lws.txn, 
                             level = int.slv.lws.lvl)

C.slv_Phy <- dplyr::distinct(C.slv_Phy_0, .keep_all = T)
C.slv_Fam <- dplyr::distinct(C.slv_Fam_0, .keep_all = T)
C.slv_Lws <- dplyr::distinct(C.slv_Lws_0, .keep_all = T)

# format data and compute CZM, this process occurs as follows:
# (1) create copies to avoid overwriting the originals
# (2) convert column FeatureID into row.names
# (3) remove unneeded columns
# (4) transpose to convert features to cols and samples to rows
# (5) replace zero counts
# (6) check proportion conversion
# (7) transpose to convert features to rows and samples to cols
# NOTE: for cmultRepl(), input data = features as cols and samples as rows

C.EBTKS_abs_pro_1 <- C.EBTKS_abs_pro_0
row.names(C.EBTKS_abs_pro_1) <- C.EBTKS_abs_pro_1$FeatureID
C.EBTKS_abs_pro_2 <- dplyr::select(C.EBTKS_abs_pro_1, -dplyr::one_of(com_col))
C.EBTKS_abs_pro_3 <- as.data.frame(t(C.EBTKS_abs_pro_2))
C.czm_EBTKS_pro_0 <- zCompositions::cmultRepl(C.EBTKS_abs_pro_3, method = "CZM",
                                              output = "prop")
print(rowSums(C.czm_EBTKS_pro_0)) # all == 1
C.czm_EBTKS_pro_1 <- as.data.frame(t(C.czm_EBTKS_pro_0))

# collapse taxonomy at the phylum, family, and lowest assignment levels
# this process occurs as follows:
# (1) create a copy to avoid overwriting the originals
# (2;3) convert FeatureID back into a col and merge with _seq_tax data.frame
# (4) retain the columns of interest for specific taxonomic levels
# (5) collapse (sum/combine) identical taxonomic assignments

C.czm_EBTKS_pro_2 <- C.czm_EBTKS_pro_1

C.czm_EBTKS_pro_2$FeatureID <- row.names(C.czm_EBTKS_pro_2)
C.czm_EBTKS_pro_tax <- merge(x = C.seq_tax, y = C.czm_EBTKS_pro_2, sort = F, 
                             by = "FeatureID")

C.Phy_EBTKS_pro_0 <- dplyr::select(C.czm_EBTKS_pro_tax, 
                                   taxon = int.slv.Phylum,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))
C.Fam_EBTKS_pro_0 <- dplyr::select(C.czm_EBTKS_pro_tax, 
                                   taxon = int.slv.Family,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))
C.Lws_EBTKS_pro_0 <- dplyr::select(C.czm_EBTKS_pro_tax, 
                                   taxon = int.slv.lws.txn,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))

C.Phy_EBTKS_pro_1 <- dplyr::summarise_all(
  dplyr::group_by(C.Phy_EBTKS_pro_0, taxon), .fun = sum)
C.Fam_EBTKS_pro_1 <- dplyr::summarise_all(
  dplyr::group_by(C.Fam_EBTKS_pro_0, taxon), .fun = sum)
C.Lws_EBTKS_pro_1 <- dplyr::summarise_all(
  dplyr::group_by(C.Lws_EBTKS_pro_0, taxon), .fun = sum)

### ************************************
### C - STEP  2 - isolate taxa for plotting and subset data.frames ----
### ************************************

# isolate taxa that will be plotted, this process occurs as follows:
# (1) create copies to avoid overwriting the originals
# (2;3) convert column taxon into row.names and then remove column taxon
# (4) transpose to convert taxa to cols and samples to rows
# (5;6) create copies and then convert SampleID back into a col
# (7) select taxa of interest from each specific taxonomic level
# NOTE: taxa of interest are:
# Phyla Bacteroidetes; Firmicutes
# Family Lachnospiraceae
# Genera Bifidobacterium; Lactobacillus

C.Phy_EBTKS_pro_2 <- as.data.frame(C.Phy_EBTKS_pro_1)
C.Fam_EBTKS_pro_2 <- as.data.frame(C.Fam_EBTKS_pro_1)
C.Lws_EBTKS_pro_2 <- as.data.frame(C.Lws_EBTKS_pro_1)

row.names(C.Phy_EBTKS_pro_2) <- C.Phy_EBTKS_pro_2$taxon
row.names(C.Fam_EBTKS_pro_2) <- C.Fam_EBTKS_pro_2$taxon
row.names(C.Lws_EBTKS_pro_2) <- C.Lws_EBTKS_pro_2$taxon
C.Phy_EBTKS_pro_3 <-dplyr::select(C.Phy_EBTKS_pro_2, -taxon)
C.Fam_EBTKS_pro_3 <-dplyr::select(C.Fam_EBTKS_pro_2, -taxon)
C.Lws_EBTKS_pro_3 <-dplyr::select(C.Lws_EBTKS_pro_2, -taxon)

C.Phy_EBTKS_pro_4 <- as.data.frame(t(C.Phy_EBTKS_pro_3))
C.Fam_EBTKS_pro_4 <- as.data.frame(t(C.Fam_EBTKS_pro_3))
C.Lws_EBTKS_pro_4 <- as.data.frame(t(C.Lws_EBTKS_pro_3))

C.Phy_EBTKS_pro_5 <- C.Phy_EBTKS_pro_4
C.Fam_EBTKS_pro_5 <- C.Fam_EBTKS_pro_4
C.Lws_EBTKS_pro_5 <- C.Lws_EBTKS_pro_4
C.Phy_EBTKS_pro_5$SampleID <- row.names(C.Phy_EBTKS_pro_5)
C.Fam_EBTKS_pro_5$SampleID <- row.names(C.Fam_EBTKS_pro_5)
C.Lws_EBTKS_pro_5$SampleID <- row.names(C.Lws_EBTKS_pro_5)

C.PhyBact_EBTKS_0 <- dplyr::select(
  C.Phy_EBTKS_pro_5, SampleID,
  dplyr::contains("Bacteroidetes", ignore.case = F))
C.PhyFirm_EBTKS_0 <- dplyr::select(
  C.Phy_EBTKS_pro_5, SampleID,
  dplyr::contains("Firmicutes", ignore.case = F))
C.FamLach_EBTKS_0 <- dplyr::select(
  C.Fam_EBTKS_pro_5, SampleID,
  dplyr::contains("Lachnospiraceae", ignore.case = F))
C.LwsBifi_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Bifidobacterium", ignore.case = F))
C.LwsLact_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Lactobacillus", ignore.case = F))

# Lactobacillus contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsLact_EBTKS_1 <- C.LwsLact_EBTKS_0
C.smpID_num <- which(names(C.LwsLact_EBTKS_1) == "SampleID")
C.LwsLact_EBTKS_1$sum <- rowSums(C.LwsLact_EBTKS_1[, -C.smpID_num])
C.LwsLact_EBTKS_2 <- dplyr::select(C.LwsLact_EBTKS_1, SampleID, 
                                   Lactobacillus = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.PhyBact_EBTKS_1 <- reshape2::melt(C.PhyBact_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.PhyFirm_EBTKS_1 <- reshape2::melt(C.PhyFirm_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.FamLach_EBTKS_1 <- reshape2::melt(C.FamLach_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsBifi_EBTKS_1 <- reshape2::melt(C.LwsBifi_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsLact_EBTKS_3 <- reshape2::melt(C.LwsLact_EBTKS_2, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")

C.PhyBact_EBTKS <- merge(x = C.PhyBact_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")
C.PhyFirm_EBTKS <- merge(x = C.PhyFirm_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")
C.FamLach_EBTKS <- merge(x = C.FamLach_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")
C.LwsBifi_EBTKS <- merge(x = C.LwsBifi_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")
C.LwsLact_EBTKS <- merge(x = C.LwsLact_EBTKS_3, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.PhyBact_s_F01 <- dplyr::filter(C.PhyBact_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.PhyBact_f_F01_0 <- dplyr::filter(C.PhyBact_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.PhyBact_EPD_F01 <- dplyr::filter(C.PhyBact_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.PhyBact_s_M00 <- dplyr::filter(C.PhyBact_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.PhyBact_f_M00_0 <- dplyr::filter(C.PhyBact_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.PhyBact_EPD_M00 <- dplyr::filter(C.PhyBact_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

C.PhyFirm_s_F01 <- dplyr::filter(C.PhyFirm_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.PhyFirm_f_F01_0 <- dplyr::filter(C.PhyFirm_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.PhyFirm_EPD_F01 <- dplyr::filter(C.PhyFirm_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.PhyFirm_s_M00 <- dplyr::filter(C.PhyFirm_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.PhyFirm_f_M00_0 <- dplyr::filter(C.PhyFirm_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.PhyFirm_EPD_M00 <- dplyr::filter(C.PhyFirm_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

C.FamLach_s_F01 <- dplyr::filter(C.FamLach_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.FamLach_f_F01_0 <- dplyr::filter(C.FamLach_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.FamLach_EPD_F01 <- dplyr::filter(C.FamLach_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.FamLach_s_M00 <- dplyr::filter(C.FamLach_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.FamLach_f_M00_0 <- dplyr::filter(C.FamLach_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.FamLach_EPD_M00 <- dplyr::filter(C.FamLach_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

C.LwsBifi_s_F01 <- dplyr::filter(C.LwsBifi_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsBifi_f_F01_0 <- dplyr::filter(C.LwsBifi_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsBifi_EPD_F01 <- dplyr::filter(C.LwsBifi_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsBifi_s_M00 <- dplyr::filter(C.LwsBifi_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsBifi_f_M00_0 <- dplyr::filter(C.LwsBifi_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsBifi_EPD_M00 <- dplyr::filter(C.LwsBifi_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

C.LwsLact_s_F01 <- dplyr::filter(C.LwsLact_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsLact_f_F01_0 <- dplyr::filter(C.LwsLact_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsLact_EPD_F01 <- dplyr::filter(C.LwsLact_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsLact_s_M00 <- dplyr::filter(C.LwsLact_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsLact_f_M00_0 <- dplyr::filter(C.LwsLact_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsLact_EPD_M00 <- dplyr::filter(C.LwsLact_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ************************************
### ^^^^C - STEP  3 - format data.frames with fecal time series samples ----
### ************************************

# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.PhyBact_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.PhyBact_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.PhyBact_f_F01_2 <- C.PhyBact_f_F01_1
C.PhyBact_f_F01_2[is.na(C.PhyBact_f_F01_2)] <- 0
C.PhyBact_f_F01 <- C.PhyBact_f_F01_2
C.PhyBact_f_F01$TimeDays <- as.numeric(C.PhyBact_f_F01$TimeDays)

C.PhyBact_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.PhyBact_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.PhyBact_f_M00_2 <- C.PhyBact_f_M00_1
C.PhyBact_f_M00_2[is.na(C.PhyBact_f_M00_2)] <- 0
C.PhyBact_f_M00_3 <- C.PhyBact_f_M00_2
C.PhyBact_f_M00_3$TimeDays <- as.numeric(C.PhyBact_f_M00_3$TimeDays)
C.PhyBact_f_M00 <- dplyr::filter(C.PhyBact_f_M00_3, !TimeDays == 4)

C.PhyFirm_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.PhyFirm_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.PhyFirm_f_F01_2 <- C.PhyFirm_f_F01_1
C.PhyFirm_f_F01_2[is.na(C.PhyFirm_f_F01_2)] <- 0
C.PhyFirm_f_F01 <- C.PhyFirm_f_F01_2
C.PhyFirm_f_F01$TimeDays <- as.numeric(C.PhyFirm_f_F01$TimeDays)

C.PhyFirm_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.PhyFirm_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.PhyFirm_f_M00_2 <- C.PhyFirm_f_M00_1
C.PhyFirm_f_M00_2[is.na(C.PhyFirm_f_M00_2)] <- 0
C.PhyFirm_f_M00_3 <- C.PhyFirm_f_M00_2
C.PhyFirm_f_M00_3$TimeDays <- as.numeric(C.PhyFirm_f_M00_3$TimeDays)
C.PhyFirm_f_M00 <- dplyr::filter(C.PhyFirm_f_M00_3, !TimeDays == 4)

C.FamLach_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.FamLach_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.FamLach_f_F01_2 <- C.FamLach_f_F01_1
C.FamLach_f_F01_2[is.na(C.FamLach_f_F01_2)] <- 0
C.FamLach_f_F01 <- C.FamLach_f_F01_2
C.FamLach_f_F01$TimeDays <- as.numeric(C.FamLach_f_F01$TimeDays)

C.FamLach_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.FamLach_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.FamLach_f_M00_2 <- C.FamLach_f_M00_1
C.FamLach_f_M00_2[is.na(C.FamLach_f_M00_2)] <- 0
C.FamLach_f_M00_3 <- C.FamLach_f_M00_2
C.FamLach_f_M00_3$TimeDays <- as.numeric(C.FamLach_f_M00_3$TimeDays)
C.FamLach_f_M00 <- dplyr::filter(C.FamLach_f_M00_3, !TimeDays == 4)

C.LwsBifi_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBifi_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBifi_f_F01_2 <- C.LwsBifi_f_F01_1
C.LwsBifi_f_F01_2[is.na(C.LwsBifi_f_F01_2)] <- 0
C.LwsBifi_f_F01 <- C.LwsBifi_f_F01_2
C.LwsBifi_f_F01$TimeDays <- as.numeric(C.LwsBifi_f_F01$TimeDays)

C.LwsBifi_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBifi_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBifi_f_M00_2 <- C.LwsBifi_f_M00_1
C.LwsBifi_f_M00_2[is.na(C.LwsBifi_f_M00_2)] <- 0
C.LwsBifi_f_M00_3 <- C.LwsBifi_f_M00_2
C.LwsBifi_f_M00_3$TimeDays <- as.numeric(C.LwsBifi_f_M00_3$TimeDays)
C.LwsBifi_f_M00 <- dplyr::filter(C.LwsBifi_f_M00_3, !TimeDays == 4)

C.LwsLact_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsLact_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsLact_f_F01_2 <- C.LwsLact_f_F01_1
C.LwsLact_f_F01_2[is.na(C.LwsLact_f_F01_2)] <- 0
C.LwsLact_f_F01 <- C.LwsLact_f_F01_2
C.LwsLact_f_F01$TimeDays <- as.numeric(C.LwsLact_f_F01$TimeDays)

C.LwsLact_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsLact_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsLact_f_M00_2 <- C.LwsLact_f_M00_1
C.LwsLact_f_M00_2[is.na(C.LwsLact_f_M00_2)] <- 0
C.LwsLact_f_M00_3 <- C.LwsLact_f_M00_2
C.LwsLact_f_M00_3$TimeDays <- as.numeric(C.LwsLact_f_M00_3$TimeDays)
C.LwsLact_f_M00 <- dplyr::filter(C.LwsLact_f_M00_3, !TimeDays == 4)

### ************************************
### C - STEP  4 - define plot parameters/customize plot aesthetics ----
### ************************************

# NOTE: for consortium group aesthetics, order = Control - Rice.bran.modified

## universal to all three plot types:
# y axis limits, breaks, labels
C.ylim <- c(0, 1.15)
C.ybrk <- c(0, 0.25, 0.5, 0.75, 1.0)
C.ylab <- c("0", "25", "50", "75", "100")

# plot labeling
C.yttl <- "proportion (%)"
C.ttl_PhyFirm <- "Phylum: Firmicutes"
C.ttl_PhyBact <- "Phylum: Bacteroidetes"
C.ttl_FamLach <- "Family: Lachnospiraceae"
C.ttl_LwsBifi <- "Genus: Bifidobacterium"
C.ttl_LwsLact <- "Genus: Lactobacillus"
C.bars_stl <- "inoculum"
C.line_stl <- "murine feces"
C.strp_stl <- "murine tissue"

## specific to barplots and lineplots:
# vectors to fill bars, color lines and color shape outlines by consortium group
C.vec_hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2])
C.vec_hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2])

## specific to barplots:
# define x axis breaks, labels, title
C.bars_xbrk <- c("CMC", "RMC")
C.bars_xlab_F01 <- c("CMC-f", "RMC-f")
C.bars_xlab_M00 <- c("CMC-m", "RMC-m")
C.bars_xttl <- "" # empty on purpose

## specific to lineplots:
# vectors for linetpes and shape fill colors by consortium group
C.vec_lne_grp <- c(lne_typ_cmc, lne_typ_rmc)
C.vec_fil_F01 <- c(hex_F01_cmc[4], hex_F01_rmc[4])
C.vec_fil_M00 <- c(hex_M00_cmc[4], hex_M00_rmc[4])

# define x axis limits, breaks, labels, title
C.line_xlim <- c(-1, 92)
C.line_xbrk <- c(0, 7, 14, 21, 28, 42, 49, 77, 91)
C.line_xlab <- c("0", "1", "2", "3", "4", "6", "7", "11", "13")
C.line_xttl <- "week"

# create a data.frame to plot dotted lines and text labels for study timeline
C.line_time <- data.frame(x = c(0, 07, 26, 45, 91), 
                          y_lne = C.ylim[2], y_lab = C.ylim[2],
                          label = c("AOM", "DSS1", "DSS2", "DSS3", "END"),
                          stringsAsFactors = F)

## specific to stripcharts:
# define x axis order/breaks, labels, and title
C.strp_xbrk <- c("C.E", "R.E", "C.P", "R.P", "C.D", "R.D")

C.strp_xlab_F01 <- c("CMC-f", "RMC-f", "CMC-f", "RMC-f", "CMC-f", "RMC-f")
C.strp_xlab_M00 <- c("CMC-m", "RMC-m", "CMC-m", "RMC-m", "CMC-m", "RMC-m")
C.strp_xttl_F01 <- " "
C.strp_xttl_M00 <- " "

# vector for shape colors by sample type; order = dist - cecum - prox
C.shp_typ <- c(shp_dtl, shp_cec, shp_prx)

# (1) vectors for shape outline colors by sample type
# (2) vectors for shape fill colors by sample type
# (3) # vector for box outline colors by sample type
# NOTE: order = defined in 'C.xord_EPD' above

C.shp_hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2], 
                   hex_F01_cmc[3], hex_F01_rmc[3],
                   hex_F01_cmc[1], hex_F01_rmc[1])
C.shp_hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2], 
                   hex_M00_cmc[3], hex_M00_rmc[3],
                   hex_M00_cmc[1], hex_M00_rmc[1])

C.shp_fil_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2], 
                   hex_F01_cmc[5], hex_F01_rmc[5],
                   hex_F01_cmc[2], hex_F01_rmc[2])
C.shp_fil_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2], 
                   hex_M00_cmc[5], hex_M00_rmc[5],
                   hex_M00_cmc[2], hex_M00_rmc[2])

C.box_hex_F01 <- c(hex_F01_cmc[5], hex_F01_rmc[5], 
                   hex_F01_cmc[5], hex_F01_rmc[5],
                   hex_F01_cmc[5], hex_F01_rmc[5])
C.box_hex_M00 <- c(hex_M00_cmc[5], hex_M00_rmc[5], 
                   hex_M00_cmc[5], hex_M00_rmc[5],
                   hex_M00_cmc[5], hex_M00_rmc[5])

# (1) define jitter w/ set random seed to ensure reproducibility of final plots
# (2) create a data.frame to plot dotted lines
# (3) create a data.frame to plot text labels for sample type
C.strip_jitr <- position_jitter(width = 0.09, seed = 1234567)
C.strip_line <- data.frame(x = c(2.5, 4.5), y_lne = C.ylim[2],
                           stringsAsFactors = F)
C.strip_text <- data.frame(x = c(1.5, 3.5, 5.5), y_lab = C.ylim[2],
                           label = c("cecum", "proximal\ncolon",
                                     "distal\ncolon"),
                           stringsAsFactors = F)

## sizing parameters and ggplot themes for all three plot types:
C.sze_ptl <- 09 # (all) - size for plot title
C.sze_stl <- 08 # (all) - size for plot subtitle
C.sze_atl <- 08 # (all) - size for axis titles (x and y)
C.sze_txb <- 07 # (bar) - size for axis text
C.sze_txl <- 08 # (line) - size for axis text
C.sze_txs <- 07 # (strip) - size for axis text
C.sze_wid_bar <- 0.42 # (bar) - size for width of bars
C.sze_tme_lne <- 0.21 # (line) - size for timeline lines
C.sze_tme_txt <- 2.02 # (line) - size for timeline text
C.sze_pad_fec <- 0.33 # (line) - size for timeline text label padding
C.sze_lne_fec <- 0.42 # (line) - size for sample lines
C.sze_pts_fec <- 1.36 # (line) - size for sample points
C.sze_stk_fec <- 0.55 # (line) - size for sample point stroke
C.sze_typ_lne <- 0.21 # (strip) - size for sample type separator lines
C.sze_typ_txt <- 1.95 # (strip) - size for sample type text
C.sze_typ_hgt <- 0.72 # (strip) - size for sample type lineheight
C.sze_pad_EPD <- 0.42 # (strip) - size for sample type label padding
C.sze_pts_EPD <- 1.00 # (strip) - size for sample points
C.sze_bxp_wid <- 0.28 # (strip) - size for width of boxes
C.sze_bxp_lwd <- 0.30 # (strip) - size for box lines
C.sze_gga_wid <- c(1.84, 4, 3) # (all) - size for plot widths for ggarrange()

# (bars) custom theme parameters
C.Dynamics <- theme(
  axis.ticks = element_line(color = greydient[1]),
  axis.text = element_text(size = C.sze_txb),
  axis.title = element_text(size = C.sze_atl),
  legend.position = "none",
  plot.title = element_text(hjust = 0, size = C.sze_ptl),
  plot.subtitle = element_text(hjust = 0.5, size = C.sze_stl),
  panel.ontop = T,
  panel.background = element_rect(fill = NA, color = NA),
  plot.background = element_rect(fill = NA, color = NA),
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]))

# (line and strip) additional custom theme parameters; used when arranging
C.LineSave <- theme(axis.text.x = element_text(size = C.sze_txl),
                    axis.text.y = element_blank(),
                    axis.title.x = element_text(size = C.sze_atl),
                    axis.title.y = element_blank(),
                    plot.title = element_text(color = NA))
C.StrpSave <- theme(axis.text.x = element_text(size = C.sze_txs),
                    axis.text.y = element_blank(),
                    axis.title.x = element_text(size = C.sze_atl),
                    axis.title.y = element_blank(),
                    plot.title = element_text(color = NA))

### ************************************
### C - STEP 5a - plotting: barplots ----
### ************************************

# PhyBact
C.gpb_PhyBact_F01 <- ggbarplot(data = C.PhyBact_s_F01, title = C.ttl_PhyBact,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_PhyBact_M00 <- ggbarplot(data = C.PhyBact_s_M00, title = C.ttl_PhyBact,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# PhyFirm
C.gpb_PhyFirm_F01 <- ggbarplot(data = C.PhyFirm_s_F01, title = C.ttl_PhyFirm,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_PhyFirm_M00 <- ggbarplot(data = C.PhyFirm_s_M00, title = C.ttl_PhyFirm,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# FamLach
C.gpb_FamLach_F01 <- ggbarplot(data = C.FamLach_s_F01, title = C.ttl_FamLach,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_FamLach_M00 <- ggbarplot(data = C.FamLach_s_M00, title = C.ttl_FamLach,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# LwsBifi
C.gpb_LwsBifi_F01 <- ggbarplot(data = C.LwsBifi_s_F01, title = C.ttl_LwsBifi,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsBifi_M00 <- ggbarplot(data = C.LwsBifi_s_M00, title = C.ttl_LwsBifi,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# LwsLact
C.gpb_LwsLact_F01 <- ggbarplot(data = C.LwsLact_s_F01, title = C.ttl_LwsLact,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsLact_M00 <- ggbarplot(data = C.LwsLact_s_M00, title = C.ttl_LwsLact,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### ************************************
### C - STEP 5b - plotting: lineplots ----
### ************************************

# PhyBact
C.gpl_PhyBact_F01 <- ggplot(data = C.PhyBact_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_PhyBact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_PhyBact_M00 <- ggplot(data = C.PhyBact_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_PhyBact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

# PhyFirm
C.gpl_PhyFirm_F01 <- ggplot(data = C.PhyFirm_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_PhyFirm, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_PhyFirm_M00 <- ggplot(data = C.PhyFirm_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_PhyFirm, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

# PhyBact
C.gpl_PhyBact_F01 <- ggplot(data = C.PhyBact_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_PhyBact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_PhyBact_M00 <- ggplot(data = C.PhyBact_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_PhyBact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

# FamLach
C.gpl_FamLach_F01 <- ggplot(data = C.FamLach_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_FamLach, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_FamLach_M00 <- ggplot(data = C.FamLach_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_FamLach, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

# LwsBifi
C.gpl_LwsBifi_F01 <- ggplot(data = C.LwsBifi_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBifi, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsBifi_M00 <- ggplot(data = C.LwsBifi_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBifi, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

# LwsLact
C.gpl_LwsLact_F01 <- ggplot(data = C.LwsLact_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsLact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsLact_M00 <- ggplot(data = C.LwsLact_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsLact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### ************************************
### C - STEP 5c - plotting: stripcharts with boxplots ----
### ************************************

# PhyBact
C.gps_PhyBact_F01 <- ggstripchart(data = C.PhyBact_EPD_F01,
                                  title = C.ttl_PhyBact, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2], 
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_PhyBact_M00 <- ggstripchart(data = C.PhyBact_EPD_M00,
                                  title = C.ttl_PhyBact, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# PhyFirm
C.gps_PhyFirm_F01 <- ggstripchart(data = C.PhyFirm_EPD_F01,
                                  title = C.ttl_PhyFirm, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_PhyFirm_M00 <- ggstripchart(data = C.PhyFirm_EPD_M00,
                                  title = C.ttl_PhyFirm, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# FamLach
C.gps_FamLach_F01 <- ggstripchart(data = C.FamLach_EPD_F01,
                                  title = C.ttl_FamLach, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_FamLach_M00 <- ggstripchart(data = C.FamLach_EPD_M00,
                                  title = C.ttl_FamLach, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# LwsBifi
C.gps_LwsBifi_F01 <- ggstripchart(data = C.LwsBifi_EPD_F01,
                                  title = C.ttl_LwsBifi, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsBifi_M00 <- ggstripchart(data = C.LwsBifi_EPD_M00,
                                  title = C.ttl_LwsBifi, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

# LwsLact
C.gps_LwsLact_F01 <- ggstripchart(data = C.LwsLact_EPD_F01,
                                  title = C.ttl_LwsLact, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsLact_M00 <- ggstripchart(data = C.LwsLact_EPD_M00,
                                  title = C.ttl_LwsLact, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "top", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### ************************************
### C - STEP 5d - plotting: arrange ----
### ************************************

# arrange bars, lines, and strips (matched by sex of human donor)

C.gga_PhyFirm_F01 <- ggarrange((C.gpb_PhyFirm_F01),
                               (C.gpl_PhyFirm_F01 + C.LineSave),
                               (C.gps_PhyFirm_F01 + C.StrpSave),
                               labels = c("A", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid,
                               align = "h")

C.gga_PhyBact_F01 <- ggarrange(C.gpb_PhyBact_F01, 
                               (C.gpl_PhyBact_F01 + C.LineSave), 
                               (C.gps_PhyBact_F01 + C.StrpSave),
                               labels = c("B", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_FamLach_F01 <- ggarrange(C.gpb_FamLach_F01, 
                               (C.gpl_FamLach_F01 + C.LineSave), 
                               (C.gps_FamLach_F01 + C.StrpSave),
                               labels = c("C", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsBifi_F01 <- ggarrange(C.gpb_LwsBifi_F01, 
                               (C.gpl_LwsBifi_F01 + C.LineSave), 
                               (C.gps_LwsBifi_F01 + C.StrpSave),
                               labels = c("D", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsLact_F01 <- ggarrange(C.gpb_LwsLact_F01, 
                               (C.gpl_LwsLact_F01 + C.LineSave), 
                               (C.gps_LwsLact_F01 + C.StrpSave),
                               labels = c("E", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

C.gga_PhyFirm_M00 <- ggarrange(C.gpb_PhyFirm_M00, 
                               (C.gpl_PhyFirm_M00 + C.LineSave), 
                               (C.gps_PhyFirm_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_PhyBact_M00 <- ggarrange(C.gpb_PhyBact_M00, 
                               (C.gpl_PhyBact_M00 + C.LineSave), 
                               (C.gps_PhyBact_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_FamLach_M00 <- ggarrange(C.gpb_FamLach_M00, 
                               (C.gpl_FamLach_M00 + C.LineSave), 
                               (C.gps_FamLach_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsBifi_M00 <- ggarrange(C.gpb_LwsBifi_M00, 
                               (C.gpl_LwsBifi_M00 + C.LineSave), 
                               (C.gps_LwsBifi_M00 + C.StrpSave),
                               labels = c("I", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsLact_M00 <- ggarrange(C.gpb_LwsLact_M00, 
                               (C.gpl_LwsLact_M00 + C.LineSave), 
                               (C.gps_LwsLact_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

# arrang arranged taxa by sex of human donor
C.gga_taxa_F01 <- ggarrange(C.gga_PhyFirm_F01, C.gga_PhyBact_F01, 
                            C.gga_FamLach_F01, 
                            C.gga_LwsBifi_F01, C.gga_LwsLact_F01,
                            labels = NULL, ncol = 1, nrow = 5, align = "h")
C.gga_taxa_M00 <- ggarrange(C.gga_PhyFirm_M00, C.gga_PhyBact_M00, 
                            C.gga_FamLach_M00, 
                            C.gga_LwsBifi_M00, C.gga_LwsLact_M00,
                            labels = NULL, ncol = 1, nrow = 5, align = "h")

# arrange arrangements
C.gga_taxa <- ggarrange(C.gga_taxa_F01, C.gga_taxa_M00,
                        labels = NULL, ncol = 2, nrow = 1, align = "h")

### ************************************
### C - WRITE OUTPUTS ----
### ************************************

# individual taxa for female donor group:
# PhyFirmicutes
C.ofv_gga_PhyFirm_F01 <- "TransFaunation/vault/plot_PhyFirmicutes_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_PhyFirm_F01, plot = C.gga_PhyFirm_F01)

# PhyBacteroidetes
C.ofv_gga_PhyBact_F01 <- "TransFaunation/vault/plot_PhyBacteroidetes_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_PhyBact_F01, plot = C.gga_PhyBact_F01)

# FamLachnospiraceae
C.ofv_gga_FamLach_F01 <- "TransFaunation/vault/plot_FamLachnospiraceae_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_FamLach_F01, plot = C.gga_FamLach_F01)

# LwsBifidobacterium
C.ofv_gga_LwsBifi_F01 <- "TransFaunation/vault/plot_LwsBifidobacterium_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_LwsBifi_F01, plot = C.gga_LwsBifi_F01)

# LwsLactobacillus
C.ofv_gga_LwsLact_F01 <- "TransFaunation/vault/plot_LwsLactobacillus_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_LwsLact_F01, plot = C.gga_LwsLact_F01)

# arranged taxa for female donor group:
C.ofv_plot_taxa_F01 <- "TransFaunation/vault/plot_taxa_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 275,
       filename = C.ofv_plot_taxa_F01, plot = C.gga_taxa_F01)

# individual taxa for males donor group:
# PhyFirmicutes
C.ofv_gga_PhyFirm_M00 <- "TransFaunation/vault/plot_PhyFirmicutes_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_PhyFirm_M00, plot = C.gga_PhyFirm_M00)

# PhyBacteroidetes
C.ofv_gga_PhyBact_M00 <- "TransFaunation/vault/plot_PhyBacteroidetes_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_PhyBact_M00, plot = C.gga_PhyBact_M00)

# FamLachnospiraceae
C.ofv_gga_FamLach_M00 <- "TransFaunation/vault/plot_FamLachnospiraceae_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_FamLach_M00, plot = C.gga_FamLach_M00)

# LwsBifidobacterium
C.ofv_gga_LwsBifi_M00 <- "TransFaunation/vault/plot_LwsBifidobacterium_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_LwsBifi_M00, plot = C.gga_LwsBifi_M00)

# LwsLactobacillus
C.ofv_gga_LwsLact_M00 <- "TransFaunation/vault/plot_LwsLactobacillus_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_gga_LwsLact_M00, plot = C.gga_LwsLact_M00)

# arranged taxa for males donor group:
C.ofv_plot_taxa_M00 <- "TransFaunation/vault/plot_taxa_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 275,
       filename = C.ofv_plot_taxa_M00, plot = C.gga_taxa_M00)

# save workspace
C.obj <- ls(pattern = "C.")
C.lst <- c(C.obj[grep(pattern = "C.", x = C.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, C.obj_from_A, COSMOS)
save(list = C.lst, file = C.ofv_wksp)

# NOTE: pay attention to the comments as some may not reflect this code ...
# ... but rather were taken from iBright_R/aldex2/aldex2.R

#### preface ----

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)
require(ggraph, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
load(file = A.ofv_wksp)
load(file = ofv_COSMOS_wksp)

C.ofv_wksp <- paste(sep = "", path_vault, "/taxa_in_CRC", "/WS_",
                    "SectionC_taxa_in_CRC.RData")

### ************************************
### ^^^C - STEP  1 - compute CZM & collapse taxonomy ----
### ************************************

# provide provenance for information gathering at end of section:
C.prov_secstep_CS1 <- "Section C - STEP 1"
C.prov_heading_CS1 <- "compute CZM & collapse taxonomy"
# ^^^C.prov_output_obj_CS1 <- "" # this object is output to the vault'
# ^^^info and prov have not actually been added to this script

# ^^^NOTE: if the environment is empty; there are some requirements:

# create a new version of the objects needed from sections A ...
# ... and create a vector naming those objects (used when saving the workspace)
C.seq_tax <- A.seq_tax
C.EBTKS_abs_pro_0 <- A.EBTKS_abs_pro
C.obj_from_A <- c("A.seq_tax", "A.EBTKS_abs_pro")
# create data.frames to store SILVA lineage and relevant assignment levels info
# (1) retain relevant columns (renaming them in the process)
# (2) reduce data.frames to retain single assignments

C.slv_Phy_0 <- dplyr::select(C.EBTKS_abs_pro_0, taxon = int.slv.Phylum)
C.slv_Fam_0 <- dplyr::select(C.EBTKS_abs_pro_0, taxon = int.slv.Family)
C.slv_Lws_0 <- dplyr::select(C.EBTKS_abs_pro_0, taxon = int.slv.lws.txn, 
                             level = int.slv.lws.lvl)

C.slv_Phy <- dplyr::distinct(C.slv_Phy_0, .keep_all = T)
C.slv_Fam <- dplyr::distinct(C.slv_Fam_0, .keep_all = T)
C.slv_Lws <- dplyr::distinct(C.slv_Lws_0, .keep_all = T)

# format data and compute CZM, this process occurs as follows:
# (1) create copies to avoid overwriting the originals
# (2) convert column FeatureID into row.names
# (3) remove unneeded columns
# (4) transpose to convert features to cols and samples to rows
# (5) replace zero counts
# (6) check proportion conversion
# (7) transpose to convert features to rows and samples to cols
# NOTE: for cmultRepl(), input data = features as cols and samples as rows

C.EBTKS_abs_pro_1 <- C.EBTKS_abs_pro_0
row.names(C.EBTKS_abs_pro_1) <- C.EBTKS_abs_pro_1$FeatureID
C.EBTKS_abs_pro_2 <- dplyr::select(C.EBTKS_abs_pro_1, -dplyr::one_of(com_col))
C.EBTKS_abs_pro_3 <- as.data.frame(t(C.EBTKS_abs_pro_2))
C.czm_EBTKS_pro_0 <- zCompositions::cmultRepl(C.EBTKS_abs_pro_3, method = "CZM",
                                              output = "prop")
print(rowSums(C.czm_EBTKS_pro_0)) # all == 1
C.czm_EBTKS_pro_1 <- as.data.frame(t(C.czm_EBTKS_pro_0))

# collapse taxonomy at the phylum, family, and lowest assignment levels
# this process occurs as follows:
# (1) create a copy to avoid overwriting the originals
# (2;3) convert FeatureID back into a col and merge with _seq_tax data.frame
# (4) retain the columns of interest for specific taxonomic levels
# (5) collapse (sum/combine) identical taxonomic assignments

C.czm_EBTKS_pro_2 <- C.czm_EBTKS_pro_1

C.czm_EBTKS_pro_2$FeatureID <- row.names(C.czm_EBTKS_pro_2)
C.czm_EBTKS_pro_tax <- merge(x = C.seq_tax, y = C.czm_EBTKS_pro_2, sort = F, 
                             by = "FeatureID")

C.Phy_EBTKS_pro_0 <- dplyr::select(C.czm_EBTKS_pro_tax, 
                                   taxon = int.slv.Phylum,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))
C.Fam_EBTKS_pro_0 <- dplyr::select(C.czm_EBTKS_pro_tax, 
                                   taxon = int.slv.Family,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))
C.Lws_EBTKS_pro_0 <- dplyr::select(C.czm_EBTKS_pro_tax, 
                                   taxon = int.slv.lws.txn,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))

C.Phy_EBTKS_pro_1 <- dplyr::summarise_all(
  dplyr::group_by(C.Phy_EBTKS_pro_0, taxon), .fun = sum)
C.Fam_EBTKS_pro_1 <- dplyr::summarise_all(
  dplyr::group_by(C.Fam_EBTKS_pro_0, taxon), .fun = sum)
C.Lws_EBTKS_pro_1 <- dplyr::summarise_all(
  dplyr::group_by(C.Lws_EBTKS_pro_0, taxon), .fun = sum)

### ************************************
### C - STEP  2 - isolate taxa for plotting and subset data.frames ----
### ************************************

# isolate taxa that will be plotted, this process occurs as follows:
# (1) create copies to avoid overwriting the originals
# (2;3) convert column taxon into row.names and then remove column taxon
# (4) transpose to convert taxa to cols and samples to rows
# (5;6) create copies and then convert SampleID back into a col

C.Phy_EBTKS_pro_2 <- as.data.frame(C.Phy_EBTKS_pro_1)
C.Fam_EBTKS_pro_2 <- as.data.frame(C.Fam_EBTKS_pro_1)
C.Lws_EBTKS_pro_2 <- as.data.frame(C.Lws_EBTKS_pro_1)

row.names(C.Phy_EBTKS_pro_2) <- C.Phy_EBTKS_pro_2$taxon
row.names(C.Fam_EBTKS_pro_2) <- C.Fam_EBTKS_pro_2$taxon
row.names(C.Lws_EBTKS_pro_2) <- C.Lws_EBTKS_pro_2$taxon
C.Phy_EBTKS_pro_3 <-dplyr::select(C.Phy_EBTKS_pro_2, -taxon)
C.Fam_EBTKS_pro_3 <-dplyr::select(C.Fam_EBTKS_pro_2, -taxon)
C.Lws_EBTKS_pro_3 <-dplyr::select(C.Lws_EBTKS_pro_2, -taxon)

C.Phy_EBTKS_pro_4 <- as.data.frame(t(C.Phy_EBTKS_pro_3))
C.Fam_EBTKS_pro_4 <- as.data.frame(t(C.Fam_EBTKS_pro_3))
C.Lws_EBTKS_pro_4 <- as.data.frame(t(C.Lws_EBTKS_pro_3))

C.Phy_EBTKS_pro_5 <- C.Phy_EBTKS_pro_4
C.Fam_EBTKS_pro_5 <- C.Fam_EBTKS_pro_4
C.Lws_EBTKS_pro_5 <- C.Lws_EBTKS_pro_4
C.Phy_EBTKS_pro_5$SampleID <- row.names(C.Phy_EBTKS_pro_5)
C.Fam_EBTKS_pro_5$SampleID <- row.names(C.Fam_EBTKS_pro_5)
C.Lws_EBTKS_pro_5$SampleID <- row.names(C.Lws_EBTKS_pro_5)

### ************************************
### C - STEP  3 - define plot parameters/customize plot aesthetics ----
### ************************************

# NOTE: for consortium group aesthetics, order = Control - Rice.bran.modified

## universal to all three plot types:
# y axis limits, breaks, labels
C.ylim <- c(0, 1.15)
C.ybrk <- c(0, 0.25, 0.5, 0.75, 1.0)
C.ylab <- c("0", "25", "50", "75", "100")

# plot labeling
C.yttl <- "proportion (%)"
C.ttl_PhyFirm <- "Phylum: Firmicutes"
C.ttl_PhyBact <- "Phylum: Bacteroidetes"
C.ttl_FamLach <- "Family: Lachnospiraceae"
C.ttl_LwsBifi <- "Genus: Bifidobacterium"
C.ttl_LwsLact <- "Genus: Lactobacillus"
C.bars_stl <- "inoculum"
C.line_stl <- "murine feces"
C.strp_stl <- "murine tissue"

## specific to barplots and lineplots:
# vectors to fill bars, color lines and color shape outlines by consortium group
C.vec_hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2])
C.vec_hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2])

## specific to barplots:
# define x axis breaks, labels, title
C.bars_xbrk <- c("CMC", "RMC")
C.bars_xlab_F01 <- c("CMC-f", "RMC-f")
C.bars_xlab_M00 <- c("CMC-m", "RMC-m")
C.bars_xttl <- "" # empty on purpose

## specific to lineplots:
# vectors for linetpes and shape fill colors by consortium group
C.vec_lne_grp <- c(lne_typ_cmc, lne_typ_rmc)
C.vec_fil_F01 <- c(hex_F01_cmc[4], hex_F01_rmc[4])
C.vec_fil_M00 <- c(hex_M00_cmc[4], hex_M00_rmc[4])

# define x axis limits, breaks, labels, title
C.line_xlim <- c(-1, 92)
C.line_xbrk <- c(0, 7, 14, 21, 28, 42, 49, 77, 91)
C.line_xlab <- c("0", "1", "2", "3", "4", "6", "7", "11", "13")
C.line_xttl <- "week"

# create a data.frame to plot dotted lines and text labels for study timeline
C.line_time <- data.frame(x = c(0, 07, 26, 45, 91), 
                          y_lne = C.ylim[2], y_lab = C.ylim[2],
                          label = c("AOM", "DSS1", "DSS2", "DSS3", "END"),
                          stringsAsFactors = F)

## specific to stripcharts:
# define x axis order/breaks, labels, and title
C.strp_xbrk <- c("C.E", "R.E", "C.P", "R.P", "C.D", "R.D")
C.strp_xlab_F01 <- c("CMC-f", "RMC-f", "CMC-f", "RMC-f", "CMC-f", "RMC-f")
C.strp_xlab_M00 <- c("CMC-m", "RMC-m", "CMC-m", "RMC-m", "CMC-m", "RMC-m")
C.strp_xttl_F01 <- " "
C.strp_xttl_M00 <- " "

# vector for shape colors by sample type; order = dist - cecum - prox
C.shp_typ <- c(shp_dtl, shp_cec, shp_prx)

# (1) vectors for shape outline colors by sample type
# (2) vectors for shape fill colors by sample type
# (3) # vector for box outline colors by sample type
# NOTE: order = defined in 'C.xord_EPD' above

C.shp_hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2], 
                   hex_F01_cmc[3], hex_F01_rmc[3],
                   hex_F01_cmc[1], hex_F01_rmc[1])
C.shp_hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2], 
                   hex_M00_cmc[3], hex_M00_rmc[3],
                   hex_M00_cmc[1], hex_M00_rmc[1])

C.shp_fil_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2], 
                   hex_F01_cmc[5], hex_F01_rmc[5],
                   hex_F01_cmc[2], hex_F01_rmc[2])
C.shp_fil_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2], 
                   hex_M00_cmc[5], hex_M00_rmc[5],
                   hex_M00_cmc[2], hex_M00_rmc[2])

C.box_hex_F01 <- c(hex_F01_cmc[5], hex_F01_rmc[5], 
                   hex_F01_cmc[5], hex_F01_rmc[5],
                   hex_F01_cmc[5], hex_F01_rmc[5])
C.box_hex_M00 <- c(hex_M00_cmc[5], hex_M00_rmc[5], 
                   hex_M00_cmc[5], hex_M00_rmc[5],
                   hex_M00_cmc[5], hex_M00_rmc[5])

# (1) define jitter w/ set random seed to ensure reproducibility of final plots
# (2) create a data.frame to plot dotted lines
# (3) create a data.frame to plot text labels for sample type
C.strip_jitr <- position_jitter(width = 0.09, seed = 1234567)
C.strip_line <- data.frame(x = c(2.5, 4.5), y_lne = C.ylim[2],
                           stringsAsFactors = F)
C.strip_text <- data.frame(x = c(1.5, 3.5, 5.5), y_lab = C.ylim[2],
                           label = c("cecum", "proximal\ncolon",
                                     "distal\ncolon"),
                           stringsAsFactors = F)

## sizing parameters and ggplot themes for all three plot types:
C.sze_ptl <- 09 # (all) - size for plot title
C.sze_stl <- 08 # (all) - size for plot subtitle
C.sze_atl <- 08 # (all) - size for axis titles (x and y)
C.sze_txb <- 07 # (bar) - size for axis text
C.sze_txl <- 08 # (line) - size for axis text
C.sze_txs <- 07 # (strip) - size for axis text
C.sze_wid_bar <- 0.42 # (bar) - size for width of bars
C.sze_tme_lne <- 0.21 # (line) - size for timeline lines
C.sze_tme_txt <- 2.02 # (line) - size for timeline text
C.sze_pad_fec <- 0.33 # (line) - size for timeline text label padding
C.sze_lne_fec <- 0.42 # (line) - size for sample lines
C.sze_pts_fec <- 1.36 # (line) - size for sample points
C.sze_stk_fec <- 0.55 # (line) - size for sample point stroke
C.sze_typ_lne <- 0.21 # (strip) - size for sample type separator lines
C.sze_typ_txt <- 1.95 # (strip) - size for sample type text
C.sze_typ_hgt <- 0.72 # (strip) - size for sample type lineheight
C.sze_pad_EPD <- 0.42 # (strip) - size for sample type label padding
C.sze_pts_EPD <- 1.00 # (strip) - size for sample points
C.sze_bxp_wid <- 0.28 # (strip) - size for width of boxes
C.sze_bxp_lwd <- 0.30 # (strip) - size for box lines
C.sze_gga_wid <- c(1.84, 4, 3) # (all) - size for plot widths for ggarrange()

# (bars) custom theme parameters
C.Dynamics <- theme(
  axis.ticks = element_line(color = greydient[1]),
  axis.text = element_text(size = C.sze_txb),
  axis.title = element_text(size = C.sze_atl),
  legend.position = "none",
  plot.title = element_text(hjust = 0, size = C.sze_ptl),
  plot.subtitle = element_text(hjust = 0.5, size = C.sze_stl),
  panel.ontop = T,
  panel.background = element_rect(fill = NA, color = NA),
  plot.background = element_rect(fill = NA, color = NA),
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]))

# (line and strip) additional custom theme parameters; used when arranging
C.LineSave <- theme(axis.text.x = element_text(size = C.sze_txl),
                    axis.text.y = element_blank(),
                    axis.title.x = element_text(size = C.sze_atl),
                    axis.title.y = element_blank(),
                    plot.title = element_text(color = NA))
C.StrpSave <- theme(axis.text.x = element_text(size = C.sze_txs),
                    axis.text.y = element_blank(),
                    axis.title.x = element_text(size = C.sze_atl),
                    axis.title.y = element_blank(),
                    plot.title = element_text(color = NA))

### ************************************
### Blautia ----
### ************************************

# enriched in CRC - Saus 2019
# positively associated with tumor burden - Baxter 2014

# plot labeling
C.ttl_LwsBlau <- "Genus: Blautia"

# select taxa of interest from each specific taxonomic level

C.LwsBlau_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Blautia", ignore.case = F))

# contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsBlau_EBTKS_1 <- C.LwsBlau_EBTKS_0
C.smpID_num <- which(names(C.LwsBlau_EBTKS_1) == "SampleID")
C.LwsBlau_EBTKS_1$sum <- rowSums(C.LwsBlau_EBTKS_1[, -C.smpID_num])
C.LwsBlau_EBTKS_2 <- dplyr::select(C.LwsBlau_EBTKS_1, SampleID,
                                   Blautia = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsBlau_EBTKS_3 <- reshape2::melt(C.LwsBlau_EBTKS_2, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsBlau_EBTKS <- merge(x = C.LwsBlau_EBTKS_3, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsBlau_s_F01 <- dplyr::filter(C.LwsBlau_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsBlau_f_F01_0 <- dplyr::filter(C.LwsBlau_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsBlau_EPD_F01 <- dplyr::filter(C.LwsBlau_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsBlau_s_M00 <- dplyr::filter(C.LwsBlau_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsBlau_f_M00_0 <- dplyr::filter(C.LwsBlau_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsBlau_EPD_M00 <- dplyr::filter(C.LwsBlau_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsBlau_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBlau_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBlau_f_F01_2 <- C.LwsBlau_f_F01_1
C.LwsBlau_f_F01_2[is.na(C.LwsBlau_f_F01_2)] <- 0
C.LwsBlau_f_F01 <- C.LwsBlau_f_F01_2
C.LwsBlau_f_F01$TimeDays <- as.numeric(C.LwsBlau_f_F01$TimeDays)

C.LwsBlau_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBlau_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBlau_f_M00_2 <- C.LwsBlau_f_M00_1
C.LwsBlau_f_M00_2[is.na(C.LwsBlau_f_M00_2)] <- 0
C.LwsBlau_f_M00_3 <- C.LwsBlau_f_M00_2
C.LwsBlau_f_M00_3$TimeDays <- as.numeric(C.LwsBlau_f_M00_3$TimeDays)
C.LwsBlau_f_M00 <- dplyr::filter(C.LwsBlau_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsBlau_F01 <- ggbarplot(data = C.LwsBlau_s_F01, title = C.ttl_LwsBlau,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsBlau_M00 <- ggbarplot(data = C.LwsBlau_s_M00, title = C.ttl_LwsBlau,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsBlau_F01 <- ggplot(data = C.LwsBlau_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBlau, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsBlau_M00 <- ggplot(data = C.LwsBlau_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBlau, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsBlau_F01 <- ggstripchart(data = C.LwsBlau_EPD_F01,
                                  title = C.ttl_LwsBlau, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsBlau_M00 <- ggstripchart(data = C.LwsBlau_EPD_M00,
                                  title = C.ttl_LwsBlau, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsBlau_F01 <- ggarrange(C.gpb_LwsBlau_F01, 
                               (C.gpl_LwsBlau_F01 + C.LineSave), 
                               (C.gps_LwsBlau_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsBlau_M00 <- ggarrange(C.gpb_LwsBlau_M00, 
                               (C.gpl_LwsBlau_M00 + C.LineSave), 
                               (C.gps_LwsBlau_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsBlau_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBlau_F01.pdf"
C.ofv_plot_LwsBlau_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBlau_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBlau_F01, plot = C.gga_LwsBlau_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBlau_M00, plot = C.gga_LwsBlau_M00)

### ************************************
### Akkermansia and A. muciniphila ----
### ************************************

# enriched in CRC - Saus 2019
# positively associated with tumor burden - Baxter 2014

# plot labeling
C.ttl_LwsAkke <- "Genus: Akkermansia"

# select taxa of interest from each specific taxonomic level

C.LwsAkke_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Akkermansia", ignore.case = F))

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsAkke_EBTKS_1 <- reshape2::melt(C.LwsAkke_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsAkke_EBTKS <- merge(x = C.LwsAkke_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsAkke_s_F01 <- dplyr::filter(C.LwsAkke_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsAkke_f_F01_0 <- dplyr::filter(C.LwsAkke_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsAkke_EPD_F01 <- dplyr::filter(C.LwsAkke_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsAkke_s_M00 <- dplyr::filter(C.LwsAkke_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsAkke_f_M00_0 <- dplyr::filter(C.LwsAkke_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsAkke_EPD_M00 <- dplyr::filter(C.LwsAkke_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsAkke_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsAkke_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsAkke_f_F01_2 <- C.LwsAkke_f_F01_1
C.LwsAkke_f_F01_2[is.na(C.LwsAkke_f_F01_2)] <- 0
C.LwsAkke_f_F01 <- C.LwsAkke_f_F01_2
C.LwsAkke_f_F01$TimeDays <- as.numeric(C.LwsAkke_f_F01$TimeDays)

C.LwsAkke_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsAkke_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsAkke_f_M00_2 <- C.LwsAkke_f_M00_1
C.LwsAkke_f_M00_2[is.na(C.LwsAkke_f_M00_2)] <- 0
C.LwsAkke_f_M00_3 <- C.LwsAkke_f_M00_2
C.LwsAkke_f_M00_3$TimeDays <- as.numeric(C.LwsAkke_f_M00_3$TimeDays)
C.LwsAkke_f_M00 <- dplyr::filter(C.LwsAkke_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsAkke_F01 <- ggbarplot(data = C.LwsAkke_s_F01, title = C.ttl_LwsAkke,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsAkke_M00 <- ggbarplot(data = C.LwsAkke_s_M00, title = C.ttl_LwsAkke,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsAkke_F01 <- ggplot(data = C.LwsAkke_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsAkke, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsAkke_M00 <- ggplot(data = C.LwsAkke_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsAkke, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsAkke_F01 <- ggstripchart(data = C.LwsAkke_EPD_F01,
                                  title = C.ttl_LwsAkke, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsAkke_M00 <- ggstripchart(data = C.LwsAkke_EPD_M00,
                                  title = C.ttl_LwsAkke, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsAkke_F01 <- ggarrange(C.gpb_LwsAkke_F01, 
                               (C.gpl_LwsAkke_F01 + C.LineSave), 
                               (C.gps_LwsAkke_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsAkke_M00 <- ggarrange(C.gpb_LwsAkke_M00, 
                               (C.gpl_LwsAkke_M00 + C.LineSave), 
                               (C.gps_LwsAkke_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsAkke_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsAkke_F01.pdf"
C.ofv_plot_LwsAkke_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsAkke_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsAkke_F01, plot = C.gga_LwsAkke_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsAkke_M00, plot = C.gga_LwsAkke_M00)

### ************************************
### Alistipes and A. finegoldi ----
### ************************************

# (A. finegoldi) enriched in CRC - Saus 2019
# (Alistipes) positively associated with tumor burden - Baxter 2014

# NOTE: Alistipes finegoldii abundance very, very low

# plot labeling
C.ttl_LwsAlis <- "Genus: Alistipes"

# select taxa of interest from each specific taxonomic level

C.LwsAlis_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Alistipes", ignore.case = F))

# contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsAlis_EBTKS_1 <- C.LwsAlis_EBTKS_0
C.smpID_num <- which(names(C.LwsAlis_EBTKS_1) == "SampleID")
C.LwsAlis_EBTKS_1$sum <- rowSums(C.LwsAlis_EBTKS_1[, -C.smpID_num])
C.LwsAlis_EBTKS_2 <- dplyr::select(C.LwsAlis_EBTKS_1, SampleID,
                                   Alistipes = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsAlis_EBTKS_3 <- reshape2::melt(C.LwsAlis_EBTKS_2, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsAlis_EBTKS <- merge(x = C.LwsAlis_EBTKS_3, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsAlis_s_F01 <- dplyr::filter(C.LwsAlis_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsAlis_f_F01_0 <- dplyr::filter(C.LwsAlis_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsAlis_EPD_F01 <- dplyr::filter(C.LwsAlis_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsAlis_s_M00 <- dplyr::filter(C.LwsAlis_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsAlis_f_M00_0 <- dplyr::filter(C.LwsAlis_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsAlis_EPD_M00 <- dplyr::filter(C.LwsAlis_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsAlis_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsAlis_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsAlis_f_F01_2 <- C.LwsAlis_f_F01_1
C.LwsAlis_f_F01_2[is.na(C.LwsAlis_f_F01_2)] <- 0
C.LwsAlis_f_F01 <- C.LwsAlis_f_F01_2
C.LwsAlis_f_F01$TimeDays <- as.numeric(C.LwsAlis_f_F01$TimeDays)

C.LwsAlis_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsAlis_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsAlis_f_M00_2 <- C.LwsAlis_f_M00_1
C.LwsAlis_f_M00_2[is.na(C.LwsAlis_f_M00_2)] <- 0
C.LwsAlis_f_M00_3 <- C.LwsAlis_f_M00_2
C.LwsAlis_f_M00_3$TimeDays <- as.numeric(C.LwsAlis_f_M00_3$TimeDays)
C.LwsAlis_f_M00 <- dplyr::filter(C.LwsAlis_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsAlis_F01 <- ggbarplot(data = C.LwsAlis_s_F01, title = C.ttl_LwsAlis,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsAlis_M00 <- ggbarplot(data = C.LwsAlis_s_M00, title = C.ttl_LwsAlis,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsAlis_F01 <- ggplot(data = C.LwsAlis_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsAlis, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsAlis_M00 <- ggplot(data = C.LwsAlis_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsAlis, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsAlis_F01 <- ggstripchart(data = C.LwsAlis_EPD_F01,
                                  title = C.ttl_LwsAlis, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsAlis_M00 <- ggstripchart(data = C.LwsAlis_EPD_M00,
                                  title = C.ttl_LwsAlis, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsAlis_F01 <- ggarrange(C.gpb_LwsAlis_F01, 
                               (C.gpl_LwsAlis_F01 + C.LineSave), 
                               (C.gps_LwsAlis_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsAlis_M00 <- ggarrange(C.gpb_LwsAlis_M00, 
                               (C.gpl_LwsAlis_M00 + C.LineSave), 
                               (C.gps_LwsAlis_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsAlis_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsAlis_F01.pdf"
C.ofv_plot_LwsAlis_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsAlis_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsAlis_F01, plot = C.gga_LwsAlis_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsAlis_M00, plot = C.gga_LwsAlis_M00)

### ************************************
### Bacteroides fragilis ----
### ************************************

# enriched in CRC - Saus 2019

# plot labeling
C.ttl_LwsBfra <- "Genus/Species: Bacteroides fragilis"

# select taxa of interest from each specific taxonomic level

C.LwsBfra_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Bacteroides fragilis", ignore.case = F))

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsBfra_EBTKS_1 <- reshape2::melt(C.LwsBfra_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsBfra_EBTKS <- merge(x = C.LwsBfra_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsBfra_s_F01 <- dplyr::filter(C.LwsBfra_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsBfra_f_F01_0 <- dplyr::filter(C.LwsBfra_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsBfra_EPD_F01 <- dplyr::filter(C.LwsBfra_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsBfra_s_M00 <- dplyr::filter(C.LwsBfra_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsBfra_f_M00_0 <- dplyr::filter(C.LwsBfra_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsBfra_EPD_M00 <- dplyr::filter(C.LwsBfra_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsBfra_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBfra_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBfra_f_F01_2 <- C.LwsBfra_f_F01_1
C.LwsBfra_f_F01_2[is.na(C.LwsBfra_f_F01_2)] <- 0
C.LwsBfra_f_F01 <- C.LwsBfra_f_F01_2
C.LwsBfra_f_F01$TimeDays <- as.numeric(C.LwsBfra_f_F01$TimeDays)

C.LwsBfra_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBfra_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBfra_f_M00_2 <- C.LwsBfra_f_M00_1
C.LwsBfra_f_M00_2[is.na(C.LwsBfra_f_M00_2)] <- 0
C.LwsBfra_f_M00_3 <- C.LwsBfra_f_M00_2
C.LwsBfra_f_M00_3$TimeDays <- as.numeric(C.LwsBfra_f_M00_3$TimeDays)
C.LwsBfra_f_M00 <- dplyr::filter(C.LwsBfra_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsBfra_F01 <- ggbarplot(data = C.LwsBfra_s_F01, title = C.ttl_LwsBfra,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsBfra_M00 <- ggbarplot(data = C.LwsBfra_s_M00, title = C.ttl_LwsBfra,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsBfra_F01 <- ggplot(data = C.LwsBfra_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBfra, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsBfra_M00 <- ggplot(data = C.LwsBfra_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBfra, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsBfra_F01 <- ggstripchart(data = C.LwsBfra_EPD_F01,
                                  title = C.ttl_LwsBfra, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsBfra_M00 <- ggstripchart(data = C.LwsBfra_EPD_M00,
                                  title = C.ttl_LwsBfra, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsBfra_F01 <- ggarrange(C.gpb_LwsBfra_F01, 
                               (C.gpl_LwsBfra_F01 + C.LineSave), 
                               (C.gps_LwsBfra_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsBfra_M00 <- ggarrange(C.gpb_LwsBfra_M00, 
                               (C.gpl_LwsBfra_M00 + C.LineSave), 
                               (C.gps_LwsBfra_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsBfra_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBfra_F01.pdf"
C.ofv_plot_LwsBfra_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBfra_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBfra_F01, plot = C.gga_LwsBfra_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBfra_M00, plot = C.gga_LwsBfra_M00)

### ************************************
### Bacteroides (excluding B. fragilis and B. uniformis) ----
### ************************************

# positively associated with tumor burden - Baxter 2014
# enriched in CRC - Saus 2019

# plot labeling
C.ttl_LwsBact <- "Genus: Bacteroides"

# select taxa of interest from each specific taxonomic level
C.LwsBact_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Bacteroides", ignore.case = F))

# remove B. fragilis and B. uniformis
C.LwsBact_EBTKS_1 <- dplyr::select(
  C.LwsBact_EBTKS_0,
  -dplyr::contains("Bacteroides fragilis", ignore.case = F),
  -dplyr::contains("Bacteroides uniformis", ignore.case = F))

# contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsBact_EBTKS_2 <- C.LwsBact_EBTKS_1
C.smpID_num <- which(names(C.LwsBact_EBTKS_2) == "SampleID")
C.LwsBact_EBTKS_2$sum <- rowSums(C.LwsBact_EBTKS_2[, -C.smpID_num])
C.LwsBact_EBTKS_3 <- dplyr::select(C.LwsBact_EBTKS_2, SampleID,
                                   Bacteroides = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsBact_EBTKS_4 <- reshape2::melt(C.LwsBact_EBTKS_3, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsBact_EBTKS <- merge(x = C.LwsBact_EBTKS_4, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsBact_s_F01 <- dplyr::filter(C.LwsBact_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsBact_f_F01_0 <- dplyr::filter(C.LwsBact_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsBact_EPD_F01 <- dplyr::filter(C.LwsBact_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsBact_s_M00 <- dplyr::filter(C.LwsBact_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsBact_f_M00_0 <- dplyr::filter(C.LwsBact_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsBact_EPD_M00 <- dplyr::filter(C.LwsBact_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsBact_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBact_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBact_f_F01_2 <- C.LwsBact_f_F01_1
C.LwsBact_f_F01_2[is.na(C.LwsBact_f_F01_2)] <- 0
C.LwsBact_f_F01 <- C.LwsBact_f_F01_2
C.LwsBact_f_F01$TimeDays <- as.numeric(C.LwsBact_f_F01$TimeDays)

C.LwsBact_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBact_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBact_f_M00_2 <- C.LwsBact_f_M00_1
C.LwsBact_f_M00_2[is.na(C.LwsBact_f_M00_2)] <- 0
C.LwsBact_f_M00_3 <- C.LwsBact_f_M00_2
C.LwsBact_f_M00_3$TimeDays <- as.numeric(C.LwsBact_f_M00_3$TimeDays)
C.LwsBact_f_M00 <- dplyr::filter(C.LwsBact_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsBact_F01 <- ggbarplot(data = C.LwsBact_s_F01, title = C.ttl_LwsBact,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsBact_M00 <- ggbarplot(data = C.LwsBact_s_M00, title = C.ttl_LwsBact,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsBact_F01 <- ggplot(data = C.LwsBact_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsBact_M00 <- ggplot(data = C.LwsBact_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBact, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsBact_F01 <- ggstripchart(data = C.LwsBact_EPD_F01,
                                  title = C.ttl_LwsBact, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsBact_M00 <- ggstripchart(data = C.LwsBact_EPD_M00,
                                  title = C.ttl_LwsBact, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsBact_F01 <- ggarrange(C.gpb_LwsBact_F01, 
                               (C.gpl_LwsBact_F01 + C.LineSave), 
                               (C.gps_LwsBact_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsBact_M00 <- ggarrange(C.gpb_LwsBact_M00, 
                               (C.gpl_LwsBact_M00 + C.LineSave), 
                               (C.gps_LwsBact_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsBact_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBact_F01.pdf"
C.ofv_plot_LwsBact_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBact_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBact_F01, plot = C.gga_LwsBact_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBact_M00, plot = C.gga_LwsBact_M00)

### ************************************
### Escherichia-Shigella ----
### ************************************

# enriched in CRC - Saus 2019
# positively associated with tumor burden - Baxter 2014

# plot labeling
C.ttl_LwsEsSh <- "Genus: Escherichia-Shigella"

# select taxa of interest from each specific taxonomic level

C.LwsEsSh_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Escherichia", ignore.case = F))

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsEsSh_EBTKS_1 <- reshape2::melt(C.LwsEsSh_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsEsSh_EBTKS <- merge(x = C.LwsEsSh_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsEsSh_s_F01 <- dplyr::filter(C.LwsEsSh_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsEsSh_f_F01_0 <- dplyr::filter(C.LwsEsSh_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsEsSh_EPD_F01 <- dplyr::filter(C.LwsEsSh_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsEsSh_s_M00 <- dplyr::filter(C.LwsEsSh_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsEsSh_f_M00_0 <- dplyr::filter(C.LwsEsSh_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsEsSh_EPD_M00 <- dplyr::filter(C.LwsEsSh_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsEsSh_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsEsSh_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsEsSh_f_F01_2 <- C.LwsEsSh_f_F01_1
C.LwsEsSh_f_F01_2[is.na(C.LwsEsSh_f_F01_2)] <- 0
C.LwsEsSh_f_F01 <- C.LwsEsSh_f_F01_2
C.LwsEsSh_f_F01$TimeDays <- as.numeric(C.LwsEsSh_f_F01$TimeDays)

C.LwsEsSh_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsEsSh_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsEsSh_f_M00_2 <- C.LwsEsSh_f_M00_1
C.LwsEsSh_f_M00_2[is.na(C.LwsEsSh_f_M00_2)] <- 0
C.LwsEsSh_f_M00_3 <- C.LwsEsSh_f_M00_2
C.LwsEsSh_f_M00_3$TimeDays <- as.numeric(C.LwsEsSh_f_M00_3$TimeDays)
C.LwsEsSh_f_M00 <- dplyr::filter(C.LwsEsSh_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsEsSh_F01 <- ggbarplot(data = C.LwsEsSh_s_F01, title = C.ttl_LwsEsSh,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsEsSh_M00 <- ggbarplot(data = C.LwsEsSh_s_M00, title = C.ttl_LwsEsSh,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsEsSh_F01 <- ggplot(data = C.LwsEsSh_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsEsSh, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsEsSh_M00 <- ggplot(data = C.LwsEsSh_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsEsSh, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsEsSh_F01 <- ggstripchart(data = C.LwsEsSh_EPD_F01,
                                  title = C.ttl_LwsEsSh, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsEsSh_M00 <- ggstripchart(data = C.LwsEsSh_EPD_M00,
                                  title = C.ttl_LwsEsSh, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsEsSh_F01 <- ggarrange(C.gpb_LwsEsSh_F01, 
                               (C.gpl_LwsEsSh_F01 + C.LineSave), 
                               (C.gps_LwsEsSh_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsEsSh_M00 <- ggarrange(C.gpb_LwsEsSh_M00, 
                               (C.gpl_LwsEsSh_M00 + C.LineSave), 
                               (C.gps_LwsEsSh_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsEsSh_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsEsSh_F01.pdf"
C.ofv_plot_LwsEsSh_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsEsSh_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsEsSh_F01, plot = C.gga_LwsEsSh_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsEsSh_M00, plot = C.gga_LwsEsSh_M00)

### ************************************
### Parabacteroides ----
### ************************************

# positively associated with tumor burden - Baxter 2014

# plot labeling
C.ttl_LwsPara <- "Genus: Parabacteroides"

# select taxa of interest from each specific taxonomic level

C.LwsPara_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Parabacteroides", ignore.case = F))

# contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsPara_EBTKS_1 <- C.LwsPara_EBTKS_0
C.smpID_num <- which(names(C.LwsPara_EBTKS_1) == "SampleID")
C.LwsPara_EBTKS_1$sum <- rowSums(C.LwsPara_EBTKS_1[, -C.smpID_num])
C.LwsPara_EBTKS_2 <- dplyr::select(C.LwsPara_EBTKS_1, SampleID,
                                   Parabacteroides = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsPara_EBTKS_3 <- reshape2::melt(C.LwsPara_EBTKS_2, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsPara_EBTKS <- merge(x = C.LwsPara_EBTKS_3, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsPara_s_F01 <- dplyr::filter(C.LwsPara_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsPara_f_F01_0 <- dplyr::filter(C.LwsPara_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsPara_EPD_F01 <- dplyr::filter(C.LwsPara_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsPara_s_M00 <- dplyr::filter(C.LwsPara_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsPara_f_M00_0 <- dplyr::filter(C.LwsPara_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsPara_EPD_M00 <- dplyr::filter(C.LwsPara_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsPara_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsPara_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsPara_f_F01_2 <- C.LwsPara_f_F01_1
C.LwsPara_f_F01_2[is.na(C.LwsPara_f_F01_2)] <- 0
C.LwsPara_f_F01 <- C.LwsPara_f_F01_2
C.LwsPara_f_F01$TimeDays <- as.numeric(C.LwsPara_f_F01$TimeDays)

C.LwsPara_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsPara_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsPara_f_M00_2 <- C.LwsPara_f_M00_1
C.LwsPara_f_M00_2[is.na(C.LwsPara_f_M00_2)] <- 0
C.LwsPara_f_M00_3 <- C.LwsPara_f_M00_2
C.LwsPara_f_M00_3$TimeDays <- as.numeric(C.LwsPara_f_M00_3$TimeDays)
C.LwsPara_f_M00 <- dplyr::filter(C.LwsPara_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsPara_F01 <- ggbarplot(data = C.LwsPara_s_F01, title = C.ttl_LwsPara,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsPara_M00 <- ggbarplot(data = C.LwsPara_s_M00, title = C.ttl_LwsPara,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsPara_F01 <- ggplot(data = C.LwsPara_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsPara, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsPara_M00 <- ggplot(data = C.LwsPara_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsPara, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsPara_F01 <- ggstripchart(data = C.LwsPara_EPD_F01,
                                  title = C.ttl_LwsPara, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsPara_M00 <- ggstripchart(data = C.LwsPara_EPD_M00,
                                  title = C.ttl_LwsPara, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsPara_F01 <- ggarrange(C.gpb_LwsPara_F01, 
                               (C.gpl_LwsPara_F01 + C.LineSave), 
                               (C.gps_LwsPara_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsPara_M00 <- ggarrange(C.gpb_LwsPara_M00, 
                               (C.gpl_LwsPara_M00 + C.LineSave), 
                               (C.gps_LwsPara_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsPara_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsPara_F01.pdf"
C.ofv_plot_LwsPara_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsPara_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsPara_F01, plot = C.gga_LwsPara_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsPara_M00, plot = C.gga_LwsPara_M00)


### CHANGE SCALE (OVERWRITES) ----

## universal to all three plot types:
# y axis limits, breaks, labels
C.ylim <- c(0, 0.23)
C.ybrk <- c(0, 0.05, 0.10, 0.15, 0.20)
C.ylab <- c("0", "5", "10", "15", "20")

# create a data.frame to plot dotted lines and text labels for study timeline
C.line_time <- data.frame(x = c(0, 07, 26, 45, 91), 
                          y_lne = C.ylim[2], y_lab = C.ylim[2],
                          label = c("AOM", "DSS1", "DSS2", "DSS3", "END"),
                          stringsAsFactors = F)

# (1) define jitter w/ set random seed to ensure reproducibility of final plots
# (2) create a data.frame to plot dotted lines
# (3) create a data.frame to plot text labels for sample type
C.strip_jitr <- position_jitter(width = 0.09, seed = 1234567)
C.strip_line <- data.frame(x = c(2.5, 4.5), y_lne = C.ylim[2],
                           stringsAsFactors = F)
C.strip_text <- data.frame(x = c(1.5, 3.5, 5.5), y_lab = C.ylim[2],
                           label = c("cecum", "colon-prox", "colon-dist"),
                           stringsAsFactors = F)

### ************************************
### Faecalibacterium prausnitzii ----
### ************************************



# depleted in CRC - Saus 2019

# plot labeling
C.ttl_LwsFaec <- "Genus: Faecalibacterium"

# select taxa of interest from each specific taxonomic level

C.LwsFaec_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Faecalibacterium", ignore.case = F))

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsFaec_EBTKS_1 <- reshape2::melt(C.LwsFaec_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsFaec_EBTKS <- merge(x = C.LwsFaec_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsFaec_s_F01 <- dplyr::filter(C.LwsFaec_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsFaec_f_F01_0 <- dplyr::filter(C.LwsFaec_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsFaec_EPD_F01 <- dplyr::filter(C.LwsFaec_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsFaec_s_M00 <- dplyr::filter(C.LwsFaec_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsFaec_f_M00_0 <- dplyr::filter(C.LwsFaec_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsFaec_EPD_M00 <- dplyr::filter(C.LwsFaec_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsFaec_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsFaec_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsFaec_f_F01_2 <- C.LwsFaec_f_F01_1
C.LwsFaec_f_F01_2[is.na(C.LwsFaec_f_F01_2)] <- 0
C.LwsFaec_f_F01 <- C.LwsFaec_f_F01_2
C.LwsFaec_f_F01$TimeDays <- as.numeric(C.LwsFaec_f_F01$TimeDays)

C.LwsFaec_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsFaec_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsFaec_f_M00_2 <- C.LwsFaec_f_M00_1
C.LwsFaec_f_M00_2[is.na(C.LwsFaec_f_M00_2)] <- 0
C.LwsFaec_f_M00_3 <- C.LwsFaec_f_M00_2
C.LwsFaec_f_M00_3$TimeDays <- as.numeric(C.LwsFaec_f_M00_3$TimeDays)
C.LwsFaec_f_M00 <- dplyr::filter(C.LwsFaec_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsFaec_F01 <- ggbarplot(data = C.LwsFaec_s_F01, title = C.ttl_LwsFaec,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsFaec_M00 <- ggbarplot(data = C.LwsFaec_s_M00, title = C.ttl_LwsFaec,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsFaec_F01 <- ggplot(data = C.LwsFaec_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsFaec, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsFaec_M00 <- ggplot(data = C.LwsFaec_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsFaec, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsFaec_F01 <- ggstripchart(data = C.LwsFaec_EPD_F01,
                                  title = C.ttl_LwsFaec, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsFaec_M00 <- ggstripchart(data = C.LwsFaec_EPD_M00,
                                  title = C.ttl_LwsFaec, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsFaec_F01 <- ggarrange(C.gpb_LwsFaec_F01, 
                               (C.gpl_LwsFaec_F01 + C.LineSave), 
                               (C.gps_LwsFaec_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsFaec_M00 <- ggarrange(C.gpb_LwsFaec_M00, 
                               (C.gpl_LwsFaec_M00 + C.LineSave), 
                               (C.gps_LwsFaec_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsFaec_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsFaec_F01.pdf"
C.ofv_plot_LwsFaec_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsFaec_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsFaec_F01, plot = C.gga_LwsFaec_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsFaec_M00, plot = C.gga_LwsFaec_M00)

### ************************************
### Roseburia ----
### ************************************

# depleted in CRC - Saus 2019

# plot labeling
C.ttl_LwsRose <- "Genus: Roseburia"

# select taxa of interest from each specific taxonomic level

C.LwsRose_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Roseburia", ignore.case = F))

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsRose_EBTKS_1 <- reshape2::melt(C.LwsRose_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsRose_EBTKS <- merge(x = C.LwsRose_EBTKS_1, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsRose_s_F01 <- dplyr::filter(C.LwsRose_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsRose_f_F01_0 <- dplyr::filter(C.LwsRose_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsRose_EPD_F01 <- dplyr::filter(C.LwsRose_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsRose_s_M00 <- dplyr::filter(C.LwsRose_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsRose_f_M00_0 <- dplyr::filter(C.LwsRose_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsRose_EPD_M00 <- dplyr::filter(C.LwsRose_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsRose_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsRose_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsRose_f_F01_2 <- C.LwsRose_f_F01_1
C.LwsRose_f_F01_2[is.na(C.LwsRose_f_F01_2)] <- 0
C.LwsRose_f_F01 <- C.LwsRose_f_F01_2
C.LwsRose_f_F01$TimeDays <- as.numeric(C.LwsRose_f_F01$TimeDays)

C.LwsRose_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsRose_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsRose_f_M00_2 <- C.LwsRose_f_M00_1
C.LwsRose_f_M00_2[is.na(C.LwsRose_f_M00_2)] <- 0
C.LwsRose_f_M00_3 <- C.LwsRose_f_M00_2
C.LwsRose_f_M00_3$TimeDays <- as.numeric(C.LwsRose_f_M00_3$TimeDays)
C.LwsRose_f_M00 <- dplyr::filter(C.LwsRose_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsRose_F01 <- ggbarplot(data = C.LwsRose_s_F01, title = C.ttl_LwsRose,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsRose_M00 <- ggbarplot(data = C.LwsRose_s_M00, title = C.ttl_LwsRose,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsRose_F01 <- ggplot(data = C.LwsRose_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsRose, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsRose_M00 <- ggplot(data = C.LwsRose_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsRose, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsRose_F01 <- ggstripchart(data = C.LwsRose_EPD_F01,
                                  title = C.ttl_LwsRose, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsRose_M00 <- ggstripchart(data = C.LwsRose_EPD_M00,
                                  title = C.ttl_LwsRose, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsRose_F01 <- ggarrange(C.gpb_LwsRose_F01, 
                               (C.gpl_LwsRose_F01 + C.LineSave), 
                               (C.gps_LwsRose_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsRose_M00 <- ggarrange(C.gpb_LwsRose_M00, 
                               (C.gpl_LwsRose_M00 + C.LineSave), 
                               (C.gps_LwsRose_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsRose_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsRose_F01.pdf"
C.ofv_plot_LwsRose_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsRose_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsRose_F01, plot = C.gga_LwsRose_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsRose_M00, plot = C.gga_LwsRose_M00)

### ************************************
### Ruminococcus spp ----
### ************************************

# depleted in CRC - Saus 2019

# plot labeling
C.ttl_LwsRumi <- "Genus: Ruminococcus"

# select taxa of interest from each specific taxonomic level
C.LwsRumi_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Ruminococcus", ignore.case = F))

# remove [Ruminococcus]
C.LwsRumi_EBTKS_1 <- dplyr::select(
  C.LwsRumi_EBTKS_0,
  -dplyr::contains("[Ruminococcus]", ignore.case = F))

# contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsRumi_EBTKS_2 <- C.LwsRumi_EBTKS_1
C.smpID_num <- which(names(C.LwsRumi_EBTKS_2) == "SampleID")
C.LwsRumi_EBTKS_2$sum <- rowSums(C.LwsRumi_EBTKS_2[, -C.smpID_num])
C.LwsRumi_EBTKS_3 <- dplyr::select(C.LwsRumi_EBTKS_2, SampleID,
                                   Ruminococcus = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsRumi_EBTKS_4 <- reshape2::melt(C.LwsRumi_EBTKS_3, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsRumi_EBTKS <- merge(x = C.LwsRumi_EBTKS_4, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsRumi_s_F01 <- dplyr::filter(C.LwsRumi_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsRumi_f_F01_0 <- dplyr::filter(C.LwsRumi_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsRumi_EPD_F01 <- dplyr::filter(C.LwsRumi_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsRumi_s_M00 <- dplyr::filter(C.LwsRumi_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsRumi_f_M00_0 <- dplyr::filter(C.LwsRumi_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsRumi_EPD_M00 <- dplyr::filter(C.LwsRumi_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsRumi_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsRumi_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsRumi_f_F01_2 <- C.LwsRumi_f_F01_1
C.LwsRumi_f_F01_2[is.na(C.LwsRumi_f_F01_2)] <- 0
C.LwsRumi_f_F01 <- C.LwsRumi_f_F01_2
C.LwsRumi_f_F01$TimeDays <- as.numeric(C.LwsRumi_f_F01$TimeDays)

C.LwsRumi_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsRumi_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsRumi_f_M00_2 <- C.LwsRumi_f_M00_1
C.LwsRumi_f_M00_2[is.na(C.LwsRumi_f_M00_2)] <- 0
C.LwsRumi_f_M00_3 <- C.LwsRumi_f_M00_2
C.LwsRumi_f_M00_3$TimeDays <- as.numeric(C.LwsRumi_f_M00_3$TimeDays)
C.LwsRumi_f_M00 <- dplyr::filter(C.LwsRumi_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsRumi_F01 <- ggbarplot(data = C.LwsRumi_s_F01, title = C.ttl_LwsRumi,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsRumi_M00 <- ggbarplot(data = C.LwsRumi_s_M00, title = C.ttl_LwsRumi,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsRumi_F01 <- ggplot(data = C.LwsRumi_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsRumi, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsRumi_M00 <- ggplot(data = C.LwsRumi_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsRumi, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsRumi_F01 <- ggstripchart(data = C.LwsRumi_EPD_F01,
                                  title = C.ttl_LwsRumi, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsRumi_M00 <- ggstripchart(data = C.LwsRumi_EPD_M00,
                                  title = C.ttl_LwsRumi, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsRumi_F01 <- ggarrange(C.gpb_LwsRumi_F01, 
                               (C.gpl_LwsRumi_F01 + C.LineSave), 
                               (C.gps_LwsRumi_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsRumi_M00 <- ggarrange(C.gpb_LwsRumi_M00, 
                               (C.gpl_LwsRumi_M00 + C.LineSave), 
                               (C.gps_LwsRumi_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsRumi_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsRumi_F01.pdf"
C.ofv_plot_LwsRumi_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsRumi_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsRumi_F01, plot = C.gga_LwsRumi_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsRumi_M00, plot = C.gga_LwsRumi_M00)

### [Ruminococcus] spp ----
### ************************************

# depleted in CRC - Saus 2019

# plot labeling
C.ttl_LwsRcon <- "Genus: [Ruminococcus]"

# select taxa of interest from each specific taxonomic level

C.LwsRcon_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("[Ruminococcus]", ignore.case = F))

# contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsRcon_EBTKS_1 <- C.LwsRcon_EBTKS_0
C.smpID_num <- which(names(C.LwsRcon_EBTKS_1) == "SampleID")
C.LwsRcon_EBTKS_1$sum <- rowSums(C.LwsRcon_EBTKS_1[, -C.smpID_num])
C.LwsRcon_EBTKS_2 <- dplyr::select(C.LwsRcon_EBTKS_1, SampleID,
                                   `[Ruminococcus]` = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsRcon_EBTKS_3 <- reshape2::melt(C.LwsRcon_EBTKS_2, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsRcon_EBTKS <- merge(x = C.LwsRcon_EBTKS_3, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsRcon_s_F01 <- dplyr::filter(C.LwsRcon_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsRcon_f_F01_0 <- dplyr::filter(C.LwsRcon_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsRcon_EPD_F01 <- dplyr::filter(C.LwsRcon_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsRcon_s_M00 <- dplyr::filter(C.LwsRcon_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsRcon_f_M00_0 <- dplyr::filter(C.LwsRcon_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsRcon_EPD_M00 <- dplyr::filter(C.LwsRcon_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsRcon_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsRcon_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsRcon_f_F01_2 <- C.LwsRcon_f_F01_1
C.LwsRcon_f_F01_2[is.na(C.LwsRcon_f_F01_2)] <- 0
C.LwsRcon_f_F01 <- C.LwsRcon_f_F01_2
C.LwsRcon_f_F01$TimeDays <- as.numeric(C.LwsRcon_f_F01$TimeDays)

C.LwsRcon_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsRcon_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsRcon_f_M00_2 <- C.LwsRcon_f_M00_1
C.LwsRcon_f_M00_2[is.na(C.LwsRcon_f_M00_2)] <- 0
C.LwsRcon_f_M00_3 <- C.LwsRcon_f_M00_2
C.LwsRcon_f_M00_3$TimeDays <- as.numeric(C.LwsRcon_f_M00_3$TimeDays)
C.LwsRcon_f_M00 <- dplyr::filter(C.LwsRcon_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsRcon_F01 <- ggbarplot(data = C.LwsRcon_s_F01, title = C.ttl_LwsRcon,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsRcon_M00 <- ggbarplot(data = C.LwsRcon_s_M00, title = C.ttl_LwsRcon,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsRcon_F01 <- ggplot(data = C.LwsRcon_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsRcon, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsRcon_M00 <- ggplot(data = C.LwsRcon_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsRcon, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsRcon_F01 <- ggstripchart(data = C.LwsRcon_EPD_F01,
                                  title = C.ttl_LwsRcon, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsRcon_M00 <- ggstripchart(data = C.LwsRcon_EPD_M00,
                                  title = C.ttl_LwsRcon, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsRcon_F01 <- ggarrange(C.gpb_LwsRcon_F01, 
                               (C.gpl_LwsRcon_F01 + C.LineSave), 
                               (C.gps_LwsRcon_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsRcon_M00 <- ggarrange(C.gpb_LwsRcon_M00, 
                               (C.gpl_LwsRcon_M00 + C.LineSave), 
                               (C.gps_LwsRcon_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsRcon_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsRcon_F01.pdf"
C.ofv_plot_LwsRcon_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsRcon_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsRcon_F01, plot = C.gga_LwsRcon_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsRcon_M00, plot = C.gga_LwsRcon_M00)

### ************************************
### Bacteroides uniformis ----
### ************************************

# depleted in CRC - Saus 2019

# plot labeling
C.ttl_LwsBuni <- "Genus/Species: Bacteroides uniformis"

# select taxa of interest from each specific taxonomic level

C.LwsBuni_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Bacteroides uniformis", ignore.case = F))

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsBuni_EBTKS_1 <- reshape2::melt(C.LwsBuni_EBTKS_0, id.vars = "SampleID",
                                    variable.name = "taxon",
                                    value.name = "prop")
C.LwsBuni_EBTKS <- merge(x = C.LwsBuni_EBTKS_1, y = smp_dat, sort = F,
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ...
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsBuni_s_F01 <- dplyr::filter(C.LwsBuni_EBTKS,
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsBuni_f_F01_0 <- dplyr::filter(C.LwsBuni_EBTKS,
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsBuni_EPD_F01 <- dplyr::filter(C.LwsBuni_EBTKS,
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal",
                                   HumanDonorSex == "Female")
C.LwsBuni_s_M00 <- dplyr::filter(C.LwsBuni_EBTKS,
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsBuni_f_M00_0 <- dplyr::filter(C.LwsBuni_EBTKS,
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsBuni_EPD_M00 <- dplyr::filter(C.LwsBuni_EBTKS,
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal",
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample...
# ... therefore calculation of standard erroe will produce an NA...
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsBuni_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBuni_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBuni_f_F01_2 <- C.LwsBuni_f_F01_1
C.LwsBuni_f_F01_2[is.na(C.LwsBuni_f_F01_2)] <- 0
C.LwsBuni_f_F01 <- C.LwsBuni_f_F01_2
C.LwsBuni_f_F01$TimeDays <- as.numeric(C.LwsBuni_f_F01$TimeDays)

C.LwsBuni_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsBuni_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsBuni_f_M00_2 <- C.LwsBuni_f_M00_1
C.LwsBuni_f_M00_2[is.na(C.LwsBuni_f_M00_2)] <- 0
C.LwsBuni_f_M00_3 <- C.LwsBuni_f_M00_2
C.LwsBuni_f_M00_3$TimeDays <- as.numeric(C.LwsBuni_f_M00_3$TimeDays)
C.LwsBuni_f_M00 <- dplyr::filter(C.LwsBuni_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsBuni_F01 <- ggbarplot(data = C.LwsBuni_s_F01, title = C.ttl_LwsBuni,
                               palette = C.vec_hex_F01, color = greydient[1],
                               font.family = fnt, width = C.sze_wid_bar,
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl,
                               fill = "ConsortiumAbrv", add = "mean_se",
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsBuni_M00 <- ggbarplot(data = C.LwsBuni_s_M00, title = C.ttl_LwsBuni,
                               palette = C.vec_hex_M00, color = greydient[1],
                               font.family = fnt, width = C.sze_wid_bar,
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl,
                               fill = "ConsortiumAbrv", add = "mean_se",
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsBuni_F01 <- ggplot(data = C.LwsBuni_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv,
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBuni, y = C.yttl,
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk,
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsBuni_M00 <- ggplot(data = C.LwsBuni_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv,
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsBuni, y = C.yttl,
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk,
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsBuni_F01 <- ggstripchart(data = C.LwsBuni_EPD_F01,
                                  title = C.ttl_LwsBuni, font.family = fnt,
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType",
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType",
                                  position = C.strip_jitr, order = C.strp_xbrk,
                                  size = C.sze_pts_EPD, add = "boxplot",
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsBuni_M00 <- ggstripchart(data = C.LwsBuni_EPD_M00,
                                  title = C.ttl_LwsBuni, font.family = fnt,
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType",
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType",
                                  position = C.strip_jitr, order = C.strp_xbrk,
                                  size = C.sze_pts_EPD, add = "boxplot",
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsBuni_F01 <- ggarrange(C.gpb_LwsBuni_F01,
                               (C.gpl_LwsBuni_F01 + C.LineSave),
                               (C.gps_LwsBuni_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid,
                               align = "h")
C.gga_LwsBuni_M00 <- ggarrange(C.gpb_LwsBuni_M00,
                               (C.gpl_LwsBuni_M00 + C.LineSave),
                               (C.gps_LwsBuni_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid,
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsBuni_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBuni_F01.pdf"
C.ofv_plot_LwsBuni_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsBuni_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBuni_F01, plot = C.gga_LwsBuni_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsBuni_M00, plot = C.gga_LwsBuni_M00)

### ************************************
### Anaerostipes ----
### ************************************

# depleted in CRC - Saus 2019

# plot labeling
C.ttl_LwsAnae <- "Genus: Anaerostipes"

# select taxa of interest from each specific taxonomic level

C.LwsAnae_EBTKS_0 <- dplyr::select(
  C.Lws_EBTKS_pro_5, SampleID,
  dplyr::contains("Anaerostipes", ignore.case = F))

# contains multiple columns (i.e. multiple lws lvl assignments)
# add them together to obtain a data.frame w/ a single column for Lactobacillus
# (1) copy to avoid overwriting the original
# (2) create an index with a numeric output to preserve column SampleID

C.LwsAnae_EBTKS_1 <- C.LwsAnae_EBTKS_0
C.smpID_num <- which(names(C.LwsAnae_EBTKS_1) == "SampleID")
C.LwsAnae_EBTKS_1$sum <- rowSums(C.LwsAnae_EBTKS_1[, -C.smpID_num])
C.LwsAnae_EBTKS_2 <- dplyr::select(C.LwsAnae_EBTKS_1, SampleID,
                                   Anaerostipes = sum)

# reshape data.frames for plotting, this process occurs as follows:
# (1) melt the data to add a new column specifying the taxon and its proportion
# (2) merge with sample data

C.LwsAnae_EBTKS_3 <- reshape2::melt(C.LwsAnae_EBTKS_2, id.vars = "SampleID",
                                    variable.name = "taxon", 
                                    value.name = "prop")
C.LwsAnae_EBTKS <- merge(x = C.LwsAnae_EBTKS_3, y = smp_dat, sort = F, 
                         by = "SampleID")

# subset by sample type and human donor (inoculum) group to obtain ... 
# ... the appropriate groups of samples that will ultimately be plotted
# NOTE: samples will be split by the sex of human donor into the following:
# _s = stool inoculum
# _EPD = mouse cEcum; mouse Proximal colon; mouse Distal colon
# _f = mouse feces

C.LwsAnae_s_F01 <- dplyr::filter(C.LwsAnae_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Female")
C.LwsAnae_f_F01_0 <- dplyr::filter(C.LwsAnae_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Female")
C.LwsAnae_EPD_F01 <- dplyr::filter(C.LwsAnae_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Female")
C.LwsAnae_s_M00 <- dplyr::filter(C.LwsAnae_EBTKS, 
                                 Type == "stool.inoculum",
                                 HumanDonorSex == "Male")
C.LwsAnae_f_M00_0 <- dplyr::filter(C.LwsAnae_EBTKS, 
                                   Type == "feces",
                                   HumanDonorSex == "Male")
C.LwsAnae_EPD_M00 <- dplyr::filter(C.LwsAnae_EBTKS, 
                                   Type == "cEcum.material" |
                                     Type == "colon.Proximal" |
                                     Type == "colon.Distal", 
                                   HumanDonorSex == "Male")

### ^^^^C - STEP  3 - format data.frames with fecal time series samples
# ^^^^the comment thread below is terrible

# for fecal samples, calculate mean and standard error
# NOTE: some timepoints have only one sample... 
# ... therefore calculation of standard erroe will produce an NA... 
# ... to combat this, replace any NA's with 0's
# also, coerce TimeDays into class numeric
# for _M00 data.frames, remove TimeDays = 4

C.LwsAnae_f_F01_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsAnae_f_F01_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsAnae_f_F01_2 <- C.LwsAnae_f_F01_1
C.LwsAnae_f_F01_2[is.na(C.LwsAnae_f_F01_2)] <- 0
C.LwsAnae_f_F01 <- C.LwsAnae_f_F01_2
C.LwsAnae_f_F01$TimeDays <- as.numeric(C.LwsAnae_f_F01$TimeDays)

C.LwsAnae_f_M00_1 <- dplyr::summarise(
  dplyr::group_by_at(C.LwsAnae_f_M00_0, dplyr::vars(TimeDays, ConsortiumAbrv)),
  mean_prop = mean(prop), SE = sd(prop) / sqrt(dplyr::n()))
C.LwsAnae_f_M00_2 <- C.LwsAnae_f_M00_1
C.LwsAnae_f_M00_2[is.na(C.LwsAnae_f_M00_2)] <- 0
C.LwsAnae_f_M00_3 <- C.LwsAnae_f_M00_2
C.LwsAnae_f_M00_3$TimeDays <- as.numeric(C.LwsAnae_f_M00_3$TimeDays)
C.LwsAnae_f_M00 <- dplyr::filter(C.LwsAnae_f_M00_3, !TimeDays == 4)

### C - STEP 5a - plotting: barplots:
C.gpb_LwsAnae_F01 <- ggbarplot(data = C.LwsAnae_s_F01, title = C.ttl_LwsAnae,
                               palette = C.vec_hex_F01, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gpb_LwsAnae_M00 <- ggbarplot(data = C.LwsAnae_s_M00, title = C.ttl_LwsAnae,
                               palette = C.vec_hex_M00, color = greydient[1], 
                               font.family = fnt, width = C.sze_wid_bar, 
                               x = "ConsortiumAbrv", y = "prop", ylab = C.yttl, 
                               fill = "ConsortiumAbrv", add = "mean_se", 
                               error.plot = "upper_errorbar") +
  labs(x = C.bars_xttl, subtitle = C.bars_stl) +
  scale_x_discrete(breaks = C.bars_xbrk, labels = C.bars_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5b - plotting: lineplots
C.gpl_LwsAnae_F01 <- ggplot(data = C.LwsAnae_f_F01,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsAnae, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_F01) +
  scale_fill_manual(values = C.vec_fil_F01) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics
C.gpl_LwsAnae_M00 <- ggplot(data = C.LwsAnae_f_M00,
                            aes(x = TimeDays, y = mean_prop,
                                group = ConsortiumAbrv, 
                                color = ConsortiumAbrv)) +
  labs(title = C.ttl_LwsAnae, y = C.yttl, 
       x = C.line_xttl, subtitle = C.line_stl) +
  geom_segment(inherit.aes = F, data = C.line_time,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[7], size = C.sze_tme_lne) +
  geom_label(inherit.aes = F, data = C.line_time, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_fec, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[3],
             vjust = "middle", hjust = "middle", size = C.sze_tme_txt) +
  geom_line(aes(linetype = ConsortiumAbrv), size = C.sze_lne_fec) +
  geom_point(aes(x = TimeDays, y = mean_prop, color = ConsortiumAbrv,
                 fill = ConsortiumAbrv), shape = shp_fec, size = C.sze_pts_fec,
             stroke = C.sze_stk_fec) +
  scale_color_manual(values = C.vec_hex_M00) +
  scale_fill_manual(values = C.vec_fil_M00) +
  scale_linetype_manual(values = C.vec_lne_grp) +
  scale_x_continuous(limits = C.line_xlim, breaks = C.line_xbrk, 
                     labels = C.line_xlab) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  theme_pubr() +
  C.Dynamics

### C - STEP 5c - plotting: stripcharts with boxplots
C.gps_LwsAnae_F01 <- ggstripchart(data = C.LwsAnae_EPD_F01,
                                  title = C.ttl_LwsAnae, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_F01,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_F01, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_F01) +
  scale_fill_manual(values = C.shp_fil_F01) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_F01) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics
C.gps_LwsAnae_M00 <- ggstripchart(data = C.LwsAnae_EPD_M00,
                                  title = C.ttl_LwsAnae, font.family = fnt, 
                                  y = "prop", shape = "TypeLab",
                                  x = "ConsortiumType", 
                                  color = "ConsortiumType",
                                  fill = "ConsortiumType", 
                                  position = C.strip_jitr, order = C.strp_xbrk, 
                                  size = C.sze_pts_EPD, add = "boxplot", 
                                  add.params = list(width = C.sze_bxp_wid,
                                                    size = C.sze_bxp_lwd,
                                                    alpha = 0,
                                                    color = C.box_hex_M00,
                                                    fill = greydient[8])) +
  labs(x = C.strp_xttl_M00, subtitle = C.strp_stl) +
  geom_segment(inherit.aes = F, data = C.strip_line,
               aes(x = x, xend = x, y = y_lne, yend = 0), linetype = "dotted",
               color = greydient[6], size = C.sze_typ_lne) +
  geom_label(inherit.aes = F, data = C.strip_text, family = fnt,
             aes(x = x, y = y_lab, label = label), label.r = unit(0, "mm"),
             label.padding = unit(C.sze_pad_EPD, "mm"), label.size = NA,
             fill = greydient[8], color = greydient[2],
             vjust = "middle", hjust = "middle", size = C.sze_typ_txt, lineheight = C.sze_typ_hgt) +
  scale_shape_manual(values = C.shp_typ) +
  scale_color_manual(values = C.shp_hex_M00) +
  scale_fill_manual(values = C.shp_fil_M00) +
  scale_x_discrete(breaks = C.strp_xbrk, labels = C.strp_xlab_M00) +
  scale_y_continuous(limits = C.ylim, breaks = C.ybrk, labels = C.ylab) +
  C.Dynamics

### C - STEP 5d - plotting: arrange

# arrange bars, lines, and strips (matched by sex of human donor)
C.gga_LwsAnae_F01 <- ggarrange(C.gpb_LwsAnae_F01, 
                               (C.gpl_LwsAnae_F01 + C.LineSave), 
                               (C.gps_LwsAnae_F01 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")
C.gga_LwsAnae_M00 <- ggarrange(C.gpb_LwsAnae_M00, 
                               (C.gpl_LwsAnae_M00 + C.LineSave), 
                               (C.gps_LwsAnae_M00 + C.StrpSave),
                               labels = c("", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = C.sze_gga_wid, 
                               align = "h")

### C - WRITE OUTPUTS

# NOTE: each plot is 55mm in height

C.ofv_plot_LwsAnae_F01 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsAnae_F01.pdf"
C.ofv_plot_LwsAnae_M00 <- "TransFaunation/vault/taxa_in_CRC/plot_LwsAnae_M00.pdf"

# relative paths from wd for outputs headed for storage in the central "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsAnae_F01, plot = C.gga_LwsAnae_F01)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = C.ofv_plot_LwsAnae_M00, plot = C.gga_LwsAnae_M00)

### C - arrange plots ----

# arrange arranged taxa by sex of human donor
C.gga_taxa_F01_part01 <- ggarrange(C.gga_LwsAkke_F01, 
                                   C.gga_LwsAlis_F01, 
                                   C.gga_LwsBact_F01, 
                                   C.gga_LwsBlau_F01, 
                                   C.gga_LwsEsSh_F01,
                                   C.gga_LwsPara_F01,
                                   C.gga_LwsBfra_F01,
                                   C.gga_LwsRose_F01,
                                   labels = "AUTO", font.label = pan_fnt,
                                   ncol = 1, nrow = 8, align = "h")
C.gga_taxa_F01_part02 <- ggarrange(C.gga_LwsRcon_F01, 
                                   C.gga_LwsBuni_F01, 
                                   C.gga_LwsAnae_F01, 
                                   C.gga_LwsFaec_F01, 
                                   C.gga_LwsRumi_F01, 
                                   labels = c("I", "J", "K", "L", "M"), 
                                   font.label = pan_fnt,
                                   ncol = 1, nrow = 5, align = "h")

C.gga_taxa_M00_part01 <- ggarrange(C.gga_LwsAkke_M00, 
                                   C.gga_LwsAlis_M00, 
                                   C.gga_LwsBact_M00, 
                                   C.gga_LwsBlau_M00, 
                                   C.gga_LwsEsSh_M00,
                                   C.gga_LwsPara_M00,
                                   C.gga_LwsBfra_M00,
                                   C.gga_LwsRose_M00,
                                   ncol = 1, nrow = 8, align = "h")
C.gga_taxa_M00_part02 <- ggarrange(C.gga_LwsRcon_M00, 
                                   C.gga_LwsBuni_M00, 
                                   C.gga_LwsAnae_M00, 
                                   C.gga_LwsFaec_M00, 
                                   C.gga_LwsRumi_M00, 
                                   labels = NULL, 
                                   ncol = 1, nrow = 5, align = "h")

### C - WRITE OUTPUTS - FINAL ----

C.ofv_plot_F01_part01 <- "TransFaunation/vault/taxa_in_CRC/plot_taxa_F01_part01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 440,
       filename = C.ofv_plot_F01_part01, plot = C.gga_taxa_F01_part01)
C.ofv_plot_F01_part02 <- "TransFaunation/vault/taxa_in_CRC/plot_taxa_F01_part02.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 275,
       filename = C.ofv_plot_F01_part02, plot = C.gga_taxa_F01_part02)


C.ofv_plot_M00_part01 <- "TransFaunation/vault/taxa_in_CRC/plot_taxa_M00_part01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 440,
       filename = C.ofv_plot_M00_part01, plot = C.gga_taxa_M00_part01)
C.ofv_plot_M00_part02 <- "TransFaunation/vault/taxa_in_CRC/plot_taxa_M00_part02.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 275,
       filename = C.ofv_plot_M00_part02, plot = C.gga_taxa_M00_part02)

# save workspace
C.obj <- ls(pattern = "C.")
C.lst <- c(C.obj[grep(pattern = "C.", x = C.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, C.obj_from_A, COSMOS)
save(list = C.lst, file = C.ofv_wksp)

### ************************************
### BEGIN Section D ----
### ************************************

# note for kdp: taken from TransFaunation/section_D.R

#### preface ----

# R packages accessed via require:
require(ALDEx2, quietly = T)
require(BiocParallel, quietly = T)
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)
require(ggraph, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
load(file = A.ofv_wksp)
load(file = S.ofv_wksp)
load(file = ofv_COSMOS_wksp)

D.ofv_wksp <- paste(sep = "", path_vault, "/WS_", "SectionD.RData")

### ************************************
### ^^^D - STEP  1 - format for ALDEx2 ----
### ************************************

# provide provenance for information gathering at end of section:
D.prov_secstep_BS1 <- "Section D - STEP 1"
D.prov_heading_BS1 <- "format for ALDEx2"
# ^^^D.prov_output_obj_BS1 <- "" # this object is output to the vault'
# ^^^info and prov have not actually been added to this script

# ^^^NOTE: if the environment is empty; there are some requirements:

# create a new version of the objects needed from sections A; S; ...
# ... and create a vector naming those objects (used when saving the workspace)
D.seq_tax <- A.seq_tax
D.EBTKS_abs_pro <- A.EBTKS_abs_pro
D.asv_s_F01_0 <- S.abs_s_F01
D.asv_E_F01_0 <- S.abs_E_F01
D.asv_P_F01_0 <- S.abs_P_F01
D.asv_D_F01_0 <- S.abs_D_F01
D.asv_E_M00_0 <- S.abs_E_M00
D.asv_s_M00_0 <- S.abs_s_M00
D.asv_P_M00_0 <- S.abs_P_M00
D.asv_D_M00_0 <- S.abs_D_M00
D.obj_from_AS <- c("A.seq_tax", "A.EBTKS_abs_pro",
                   "S.abs_s_F01", "S.abs_E_F01", "S.abs_P_F01", "S.abs_D_F01", 
                   "S.abs_E_M00", "S.abs_s_M00", "S.abs_P_M00", "S.abs_D_M00")

# create data.frames to store SILVA lineage and relevant assignment levels info
# (1) retain relevant columns (renaming them in the process)
# (2) reduce data.frames to retain single assignments

D.slv_phy_0 <- dplyr::select(D.EBTKS_abs_pro, taxon = int.slv.Phylum)
D.slv_fam_0 <- dplyr::select(D.EBTKS_abs_pro, taxon = int.slv.Family)
D.slv_lws_0 <- dplyr::select(D.EBTKS_abs_pro, taxon = int.slv.lws.txn, 
                             level = int.slv.lws.lvl)

D.slv_phy <- dplyr::distinct(D.slv_phy_0, .keep_all = T)
D.slv_fam <- dplyr::distinct(D.slv_fam_0, .keep_all = T)
D.slv_lws <- dplyr::distinct(D.slv_lws_0, .keep_all = T)

# data will be formatted for testing from four perspectives:
# (asv) all ASVs present in the respective feature table
# (phy) collapsed at the Phylum level for the SILVA database
# (fam) collapsed at the Family level for the SILVA database
# (lws) collapsed at the lowest assigned level for the SILVA database

# for phy/fam/lws, collapse taxonomy and subset tables as follows:
# (1) retain relevant cols (renaming them in the process)
# (2) collapse (sum/combine) identical assignments at the respective level
# (3) subset by human donor (inoculum) and then by sample type
# (4) subset by sample type

D.phy_EBTKS_abs_0 <- dplyr::select(D.EBTKS_abs_pro, taxon = int.slv.Phylum,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))
D.fam_EBTKS_abs_0 <- dplyr::select(D.EBTKS_abs_pro, taxon = int.slv.Family,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))
D.lws_EBTKS_abs_0 <- dplyr::select(D.EBTKS_abs_pro, taxon = int.slv.lws.txn,
                                   dplyr::contains("F0", ignore.case = F),
                                   dplyr::contains("M0", ignore.case = F))

D.phy_EBTKS_abs_1 <- dplyr::summarise_all(
  dplyr::group_by(D.phy_EBTKS_abs_0, taxon), .fun = sum)
D.fam_EBTKS_abs_1 <- dplyr::summarise_all(
  dplyr::group_by(D.fam_EBTKS_abs_0, taxon), .fun = sum)
D.lws_EBTKS_abs_1 <- dplyr::summarise_all(
  dplyr::group_by(D.lws_EBTKS_abs_0, taxon), .fun = sum)

D.phy_F01 <- dplyr::select(D.phy_EBTKS_abs_1, taxon,
                           dplyr::contains("F0", ignore.case = F))
D.phy_M00 <- dplyr::select(D.phy_EBTKS_abs_1, taxon,
                           dplyr::contains("M0", ignore.case = F))
D.fam_F01 <- dplyr::select(D.fam_EBTKS_abs_1, taxon,
                           dplyr::contains("F0", ignore.case = F))
D.fam_M00 <- dplyr::select(D.fam_EBTKS_abs_1, taxon,
                           dplyr::contains("M0", ignore.case = F))
D.lws_F01 <- dplyr::select(D.lws_EBTKS_abs_1, taxon,
                           dplyr::contains("F0", ignore.case = F))
D.lws_M00 <- dplyr::select(D.lws_EBTKS_abs_1, taxon,
                           dplyr::contains("M0", ignore.case = F))

D.phy_s_F01_0 <- dplyr::select(D.phy_F01, taxon,
                               dplyr::contains("si", ignore.case = F))
D.phy_E_F01_0 <- dplyr::select(D.phy_F01, taxon,
                               dplyr::contains("CE", ignore.case = F))
D.phy_P_F01_0 <- dplyr::select(D.phy_F01, taxon,
                               dplyr::contains("CO_P", ignore.case = F))
D.phy_D_F01_0 <- dplyr::select(D.phy_F01, taxon,
                               dplyr::contains("CO_D", ignore.case = F))
D.phy_s_M00_0 <- dplyr::select(D.phy_M00, taxon,
                               dplyr::contains("si", ignore.case = F))
D.phy_E_M00_0 <- dplyr::select(D.phy_M00, taxon,
                               dplyr::contains("CE", ignore.case = F))
D.phy_P_M00_0 <- dplyr::select(D.phy_M00, taxon,
                               dplyr::contains("CO_P", ignore.case = F))
D.phy_D_M00_0 <- dplyr::select(D.phy_M00, taxon,
                               dplyr::contains("CO_D", ignore.case = F))
D.fam_s_F01_0 <- dplyr::select(D.fam_F01, taxon,
                               dplyr::contains("si", ignore.case = F))
D.fam_E_F01_0 <- dplyr::select(D.fam_F01, taxon,
                               dplyr::contains("CE", ignore.case = F))
D.fam_P_F01_0 <- dplyr::select(D.fam_F01, taxon,
                               dplyr::contains("CO_P", ignore.case = F))
D.fam_D_F01_0 <- dplyr::select(D.fam_F01, taxon,
                               dplyr::contains("CO_D", ignore.case = F))
D.fam_s_M00_0 <- dplyr::select(D.fam_M00, taxon,
                               dplyr::contains("si", ignore.case = F))
D.fam_E_M00_0 <- dplyr::select(D.fam_M00, taxon,
                               dplyr::contains("CE", ignore.case = F))
D.fam_P_M00_0 <- dplyr::select(D.fam_M00, taxon,
                               dplyr::contains("CO_P", ignore.case = F))
D.fam_D_M00_0 <- dplyr::select(D.fam_M00, taxon,
                               dplyr::contains("CO_D", ignore.case = F))
D.lws_s_F01_0 <- dplyr::select(D.lws_F01, taxon,
                               dplyr::contains("si", ignore.case = F))
D.lws_E_F01_0 <- dplyr::select(D.lws_F01, taxon,
                               dplyr::contains("CE", ignore.case = F))
D.lws_P_F01_0 <- dplyr::select(D.lws_F01, taxon,
                               dplyr::contains("CO_P", ignore.case = F))
D.lws_D_F01_0 <- dplyr::select(D.lws_F01, taxon,
                               dplyr::contains("CO_D", ignore.case = F))
D.lws_s_M00_0 <- dplyr::select(D.lws_M00, taxon,
                               dplyr::contains("si", ignore.case = F))
D.lws_E_M00_0 <- dplyr::select(D.lws_M00, taxon,
                               dplyr::contains("CE", ignore.case = F))
D.lws_P_M00_0 <- dplyr::select(D.lws_M00, taxon,
                               dplyr::contains("CO_P", ignore.case = F))
D.lws_D_M00_0 <- dplyr::select(D.lws_M00, taxon,
                               dplyr::contains("CO_D", ignore.case = F))

# format to obtain conditions used for ALDEx2, this process occurs as follows:
# (1) create copies to avoid overwriting the originals
# (2) convert column FeatureID or column taxon into row.names
# (3) remove unneeded columns
# (4) create vectors defining the pairwise conditions to be tested
# NOTE: vectors must be in exact order as the cols in the respective data.frames

D.asv_s_F01_1 <- D.asv_s_F01_0
D.asv_E_F01_1 <- D.asv_E_F01_0
D.asv_P_F01_1 <- D.asv_P_F01_0
D.asv_D_F01_1 <- D.asv_D_F01_0
D.asv_s_M00_1 <- D.asv_s_M00_0
D.asv_E_M00_1 <- D.asv_E_M00_0
D.asv_P_M00_1 <- D.asv_P_M00_0
D.asv_D_M00_1 <- D.asv_D_M00_0
D.phy_s_F01_1 <- as.data.frame(D.phy_s_F01_0)
D.phy_E_F01_1 <- as.data.frame(D.phy_E_F01_0)
D.phy_P_F01_1 <- as.data.frame(D.phy_P_F01_0)
D.phy_D_F01_1 <- as.data.frame(D.phy_D_F01_0)
D.phy_s_M00_1 <- as.data.frame(D.phy_s_M00_0)
D.phy_E_M00_1 <- as.data.frame(D.phy_E_M00_0)
D.phy_P_M00_1 <- as.data.frame(D.phy_P_M00_0)
D.phy_D_M00_1 <- as.data.frame(D.phy_D_M00_0)
D.fam_s_F01_1 <- as.data.frame(D.fam_s_F01_0)
D.fam_E_F01_1 <- as.data.frame(D.fam_E_F01_0)
D.fam_P_F01_1 <- as.data.frame(D.fam_P_F01_0)
D.fam_D_F01_1 <- as.data.frame(D.fam_D_F01_0)
D.fam_s_M00_1 <- as.data.frame(D.fam_s_M00_0)
D.fam_E_M00_1 <- as.data.frame(D.fam_E_M00_0)
D.fam_P_M00_1 <- as.data.frame(D.fam_P_M00_0)
D.fam_D_M00_1 <- as.data.frame(D.fam_D_M00_0)
D.lws_s_F01_1 <- as.data.frame(D.lws_s_F01_0)
D.lws_E_F01_1 <- as.data.frame(D.lws_E_F01_0)
D.lws_P_F01_1 <- as.data.frame(D.lws_P_F01_0)
D.lws_D_F01_1 <- as.data.frame(D.lws_D_F01_0)
D.lws_s_M00_1 <- as.data.frame(D.lws_s_M00_0)
D.lws_E_M00_1 <- as.data.frame(D.lws_E_M00_0)
D.lws_P_M00_1 <- as.data.frame(D.lws_P_M00_0)
D.lws_D_M00_1 <- as.data.frame(D.lws_D_M00_0)

row.names(D.asv_s_F01_1) <- D.asv_s_F01_1$FeatureID
row.names(D.asv_E_F01_1) <- D.asv_E_F01_1$FeatureID
row.names(D.asv_P_F01_1) <- D.asv_P_F01_1$FeatureID
row.names(D.asv_D_F01_1) <- D.asv_D_F01_1$FeatureID
row.names(D.asv_s_M00_1) <- D.asv_s_M00_1$FeatureID
row.names(D.asv_E_M00_1) <- D.asv_E_M00_1$FeatureID
row.names(D.asv_P_M00_1) <- D.asv_P_M00_1$FeatureID
row.names(D.asv_D_M00_1) <- D.asv_D_M00_1$FeatureID
row.names(D.phy_s_F01_1) <- D.phy_s_F01_1$taxon
row.names(D.phy_E_F01_1) <- D.phy_E_F01_1$taxon
row.names(D.phy_P_F01_1) <- D.phy_P_F01_1$taxon
row.names(D.phy_D_F01_1) <- D.phy_D_F01_1$taxon
row.names(D.phy_s_M00_1) <- D.phy_s_M00_1$taxon
row.names(D.phy_E_M00_1) <- D.phy_E_M00_1$taxon
row.names(D.phy_P_M00_1) <- D.phy_P_M00_1$taxon
row.names(D.phy_D_M00_1) <- D.phy_D_M00_1$taxon
row.names(D.fam_s_F01_1) <- D.fam_s_F01_1$taxon
row.names(D.fam_E_F01_1) <- D.fam_E_F01_1$taxon
row.names(D.fam_P_F01_1) <- D.fam_P_F01_1$taxon
row.names(D.fam_D_F01_1) <- D.fam_D_F01_1$taxon
row.names(D.fam_s_M00_1) <- D.fam_s_M00_1$taxon
row.names(D.fam_E_M00_1) <- D.fam_E_M00_1$taxon
row.names(D.fam_P_M00_1) <- D.fam_P_M00_1$taxon
row.names(D.fam_D_M00_1) <- D.fam_D_M00_1$taxon
row.names(D.lws_s_F01_1) <- D.lws_s_F01_1$taxon
row.names(D.lws_E_F01_1) <- D.lws_E_F01_1$taxon
row.names(D.lws_P_F01_1) <- D.lws_P_F01_1$taxon
row.names(D.lws_D_F01_1) <- D.lws_D_F01_1$taxon
row.names(D.lws_s_M00_1) <- D.lws_s_M00_1$taxon
row.names(D.lws_E_M00_1) <- D.lws_E_M00_1$taxon
row.names(D.lws_P_M00_1) <- D.lws_P_M00_1$taxon
row.names(D.lws_D_M00_1) <- D.lws_D_M00_1$taxon

D.asv_s_F01 <- dplyr::select(D.asv_s_F01_1, -dplyr::one_of(com_col))
D.asv_E_F01 <- dplyr::select(D.asv_E_F01_1, -dplyr::one_of(com_col))
D.asv_P_F01 <- dplyr::select(D.asv_P_F01_1, -dplyr::one_of(com_col))
D.asv_D_F01 <- dplyr::select(D.asv_D_F01_1, -dplyr::one_of(com_col))
D.asv_s_M00 <- dplyr::select(D.asv_s_M00_1, -dplyr::one_of(com_col))
D.asv_E_M00 <- dplyr::select(D.asv_E_M00_1, -dplyr::one_of(com_col))
D.asv_P_M00 <- dplyr::select(D.asv_P_M00_1, -dplyr::one_of(com_col))
D.asv_D_M00 <- dplyr::select(D.asv_D_M00_1, -dplyr::one_of(com_col))
D.phy_s_F01 <- dplyr::select(D.phy_s_F01_1, -taxon)
D.phy_E_F01 <- dplyr::select(D.phy_E_F01_1, -taxon)
D.phy_P_F01 <- dplyr::select(D.phy_P_F01_1, -taxon)
D.phy_D_F01 <- dplyr::select(D.phy_D_F01_1, -taxon)
D.phy_s_M00 <- dplyr::select(D.phy_s_M00_1, -taxon)
D.phy_E_M00 <- dplyr::select(D.phy_E_M00_1, -taxon)
D.phy_P_M00 <- dplyr::select(D.phy_P_M00_1, -taxon)
D.phy_D_M00 <- dplyr::select(D.phy_D_M00_1, -taxon)
D.fam_s_F01 <- dplyr::select(D.fam_s_F01_1, -taxon)
D.fam_E_F01 <- dplyr::select(D.fam_E_F01_1, -taxon)
D.fam_P_F01 <- dplyr::select(D.fam_P_F01_1, -taxon)
D.fam_D_F01 <- dplyr::select(D.fam_D_F01_1, -taxon)
D.fam_s_M00 <- dplyr::select(D.fam_s_M00_1, -taxon)
D.fam_E_M00 <- dplyr::select(D.fam_E_M00_1, -taxon)
D.fam_P_M00 <- dplyr::select(D.fam_P_M00_1, -taxon)
D.fam_D_M00 <- dplyr::select(D.fam_D_M00_1, -taxon)
D.lws_s_F01 <- dplyr::select(D.lws_s_F01_1, -taxon)
D.lws_E_F01 <- dplyr::select(D.lws_E_F01_1, -taxon)
D.lws_P_F01 <- dplyr::select(D.lws_P_F01_1, -taxon)
D.lws_D_F01 <- dplyr::select(D.lws_D_F01_1, -taxon)
D.lws_s_M00 <- dplyr::select(D.lws_s_M00_1, -taxon)
D.lws_E_M00 <- dplyr::select(D.lws_E_M00_1, -taxon)
D.lws_P_M00 <- dplyr::select(D.lws_P_M00_1, -taxon)
D.lws_D_M00 <- dplyr::select(D.lws_D_M00_1, -taxon)

D.asv_s_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.asv_s_F01)))),
                     rep("F01R", length(grep("F01R", names(D.asv_s_F01)))))
D.asv_E_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.asv_E_F01)))),
                     rep("F01R", length(grep("F01R", names(D.asv_E_F01)))))
D.asv_P_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.asv_P_F01)))),
                     rep("F01R", length(grep("F01R", names(D.asv_P_F01)))))
D.asv_D_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.asv_D_F01)))),
                     rep("F01R", length(grep("F01R", names(D.asv_D_F01)))))
D.asv_s_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.asv_s_M00)))),
                     rep("M02R", length(grep("M02R", names(D.asv_s_M00)))))
D.asv_E_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.asv_E_M00)))),
                     rep("M02R", length(grep("M02R", names(D.asv_E_M00)))))
D.asv_P_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.asv_P_M00)))),
                     rep("M02R", length(grep("M02R", names(D.asv_P_M00)))))
D.asv_D_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.asv_D_M00)))),
                     rep("M02R", length(grep("M02R", names(D.asv_D_M00)))))
D.phy_s_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.phy_s_F01)))),
                     rep("F01R", length(grep("F01R", names(D.phy_s_F01)))))
D.phy_E_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.phy_E_F01)))),
                     rep("F01R", length(grep("F01R", names(D.phy_E_F01)))))
D.phy_P_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.phy_P_F01)))),
                     rep("F01R", length(grep("F01R", names(D.phy_P_F01)))))
D.phy_D_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.phy_D_F01)))),
                     rep("F01R", length(grep("F01R", names(D.phy_D_F01)))))
D.phy_s_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.phy_s_M00)))),
                     rep("M02R", length(grep("M02R", names(D.phy_s_M00)))))
D.phy_E_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.phy_E_M00)))),
                     rep("M02R", length(grep("M02R", names(D.phy_E_M00)))))
D.phy_P_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.phy_P_M00)))),
                     rep("M02R", length(grep("M02R", names(D.phy_P_M00)))))
D.phy_D_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.phy_D_M00)))),
                     rep("M02R", length(grep("M02R", names(D.phy_D_M00)))))
D.fam_s_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.fam_s_F01)))),
                     rep("F01R", length(grep("F01R", names(D.fam_s_F01)))))
D.fam_E_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.fam_E_F01)))),
                     rep("F01R", length(grep("F01R", names(D.fam_E_F01)))))
D.fam_P_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.fam_P_F01)))),
                     rep("F01R", length(grep("F01R", names(D.fam_P_F01)))))
D.fam_D_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.fam_D_F01)))),
                     rep("F01R", length(grep("F01R", names(D.fam_D_F01)))))
D.fam_s_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.fam_s_M00)))),
                     rep("M02R", length(grep("M02R", names(D.fam_s_M00)))))
D.fam_E_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.fam_E_M00)))),
                     rep("M02R", length(grep("M02R", names(D.fam_E_M00)))))
D.fam_P_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.fam_P_M00)))),
                     rep("M02R", length(grep("M02R", names(D.fam_P_M00)))))
D.fam_D_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.fam_D_M00)))),
                     rep("M02R", length(grep("M02R", names(D.fam_D_M00)))))
D.lws_s_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.lws_s_F01)))),
                     rep("F01R", length(grep("F01R", names(D.lws_s_F01)))))
D.lws_E_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.lws_E_F01)))),
                     rep("F01R", length(grep("F01R", names(D.lws_E_F01)))))
D.lws_P_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.lws_P_F01)))),
                     rep("F01R", length(grep("F01R", names(D.lws_P_F01)))))
D.lws_D_F01_cnd <- c(rep("F01C", length(grep("F01C", names(D.lws_D_F01)))),
                     rep("F01R", length(grep("F01R", names(D.lws_D_F01)))))
D.lws_s_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.lws_s_M00)))),
                     rep("M02R", length(grep("M02R", names(D.lws_s_M00)))))
D.lws_E_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.lws_E_M00)))),
                     rep("M02R", length(grep("M02R", names(D.lws_E_M00)))))
D.lws_P_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.lws_P_M00)))),
                     rep("M02R", length(grep("M02R", names(D.lws_P_M00)))))
D.lws_D_M00_cnd <- c(rep("M01C", length(grep("M01C", names(D.lws_D_M00)))),
                     rep("M02R", length(grep("M02R", names(D.lws_D_M00)))))

### ************************************
### D - STEP  2 - perform testing ----
### ************************************

# this process occurs as follows:
# (1) clr transform n = 1000 monte carlo (mc) instances
# (2) perform statistical tests and adjust p-values for multiple comparisons
# (3) calculate effect sizes for each mc instance, report the expected values
# NOTE: for stool inoculum samples from human female, the tests must be paired
# NOTE: 1-2-3 took ~10 minutes minute on my machine

D.asv_clr_s_F01 <- aldex.clr(D.asv_s_F01, conds = D.asv_s_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.asv_clr_E_F01 <- aldex.clr(D.asv_E_F01, conds = D.asv_E_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.asv_clr_P_F01 <- aldex.clr(D.asv_P_F01, conds = D.asv_P_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.asv_clr_D_F01 <- aldex.clr(D.asv_D_F01, conds = D.asv_D_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.asv_clr_s_M00 <- aldex.clr(D.asv_s_M00, conds = D.asv_s_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.asv_clr_E_M00 <- aldex.clr(D.asv_E_M00, conds = D.asv_E_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.asv_clr_P_M00 <- aldex.clr(D.asv_P_M00, conds = D.asv_P_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.asv_clr_D_M00 <- aldex.clr(D.asv_D_M00, conds = D.asv_D_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_s_F01 <- aldex.clr(D.phy_s_F01, conds = D.phy_s_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_E_F01 <- aldex.clr(D.phy_E_F01, conds = D.phy_E_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_P_F01 <- aldex.clr(D.phy_P_F01, conds = D.phy_P_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_D_F01 <- aldex.clr(D.phy_D_F01, conds = D.phy_D_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_s_M00 <- aldex.clr(D.phy_s_M00, conds = D.phy_s_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_E_M00 <- aldex.clr(D.phy_E_M00, conds = D.phy_E_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_P_M00 <- aldex.clr(D.phy_P_M00, conds = D.phy_P_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.phy_clr_D_M00 <- aldex.clr(D.phy_D_M00, conds = D.phy_D_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_s_F01 <- aldex.clr(D.fam_s_F01, conds = D.fam_s_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_E_F01 <- aldex.clr(D.fam_E_F01, conds = D.fam_E_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_P_F01 <- aldex.clr(D.fam_P_F01, conds = D.fam_P_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_D_F01 <- aldex.clr(D.fam_D_F01, conds = D.fam_D_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_s_M00 <- aldex.clr(D.fam_s_M00, conds = D.fam_s_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_E_M00 <- aldex.clr(D.fam_E_M00, conds = D.fam_E_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_P_M00 <- aldex.clr(D.fam_P_M00, conds = D.fam_P_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.fam_clr_D_M00 <- aldex.clr(D.fam_D_M00, conds = D.fam_D_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_s_F01 <- aldex.clr(D.lws_s_F01, conds = D.lws_s_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_E_F01 <- aldex.clr(D.lws_E_F01, conds = D.lws_E_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_P_F01 <- aldex.clr(D.lws_P_F01, conds = D.lws_P_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_D_F01 <- aldex.clr(D.lws_D_F01, conds = D.lws_D_F01_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_s_M00 <- aldex.clr(D.lws_s_M00, conds = D.lws_s_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_E_M00 <- aldex.clr(D.lws_E_M00, conds = D.lws_E_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_P_M00 <- aldex.clr(D.lws_P_M00, conds = D.lws_P_M00_cnd, useMC = T,
                             mc.samples = 1000)
D.lws_clr_D_M00 <- aldex.clr(D.lws_D_M00, conds = D.lws_D_M00_cnd, useMC = T,
                             mc.samples = 1000)

D.asv_tst_s_F01 <- aldex.ttest(D.asv_clr_s_F01, paired.test = T)
D.asv_tst_E_F01 <- aldex.ttest(D.asv_clr_E_F01, paired.test = F)
D.asv_tst_P_F01 <- aldex.ttest(D.asv_clr_P_F01, paired.test = F)
D.asv_tst_D_F01 <- aldex.ttest(D.asv_clr_D_F01, paired.test = F)
D.asv_tst_s_M00 <- aldex.ttest(D.asv_clr_s_M00, paired.test = F)
D.asv_tst_E_M00 <- aldex.ttest(D.asv_clr_E_M00, paired.test = F)
D.asv_tst_P_M00 <- aldex.ttest(D.asv_clr_P_M00, paired.test = F)
D.asv_tst_D_M00 <- aldex.ttest(D.asv_clr_D_M00, paired.test = F)
D.phy_tst_s_F01 <- aldex.ttest(D.phy_clr_s_F01, paired.test = T)
D.phy_tst_E_F01 <- aldex.ttest(D.phy_clr_E_F01, paired.test = F)
D.phy_tst_P_F01 <- aldex.ttest(D.phy_clr_P_F01, paired.test = F)
D.phy_tst_D_F01 <- aldex.ttest(D.phy_clr_D_F01, paired.test = F)
D.phy_tst_s_M00 <- aldex.ttest(D.phy_clr_s_M00, paired.test = F)
D.phy_tst_E_M00 <- aldex.ttest(D.phy_clr_E_M00, paired.test = F)
D.phy_tst_P_M00 <- aldex.ttest(D.phy_clr_P_M00, paired.test = F)
D.phy_tst_D_M00 <- aldex.ttest(D.phy_clr_D_M00, paired.test = F)
D.fam_tst_s_F01 <- aldex.ttest(D.fam_clr_s_F01, paired.test = T)
D.fam_tst_E_F01 <- aldex.ttest(D.fam_clr_E_F01, paired.test = F)
D.fam_tst_P_F01 <- aldex.ttest(D.fam_clr_P_F01, paired.test = F)
D.fam_tst_D_F01 <- aldex.ttest(D.fam_clr_D_F01, paired.test = F)
D.fam_tst_s_M00 <- aldex.ttest(D.fam_clr_s_M00, paired.test = F)
D.fam_tst_E_M00 <- aldex.ttest(D.fam_clr_E_M00, paired.test = F)
D.fam_tst_P_M00 <- aldex.ttest(D.fam_clr_P_M00, paired.test = F)
D.fam_tst_D_M00 <- aldex.ttest(D.fam_clr_D_M00, paired.test = F)
D.lws_tst_s_F01 <- aldex.ttest(D.lws_clr_s_F01, paired.test = T)
D.lws_tst_E_F01 <- aldex.ttest(D.lws_clr_E_F01, paired.test = F)
D.lws_tst_P_F01 <- aldex.ttest(D.lws_clr_P_F01, paired.test = F)
D.lws_tst_D_F01 <- aldex.ttest(D.lws_clr_D_F01, paired.test = F)
D.lws_tst_s_M00 <- aldex.ttest(D.lws_clr_s_M00, paired.test = F)
D.lws_tst_E_M00 <- aldex.ttest(D.lws_clr_E_M00, paired.test = F)
D.lws_tst_P_M00 <- aldex.ttest(D.lws_clr_P_M00, paired.test = F)
D.lws_tst_D_M00 <- aldex.ttest(D.lws_clr_D_M00, paired.test = F)

D.asv_eff_s_F01 <- aldex.effect(D.asv_clr_s_F01, include.sample.summary = T, 
                                useMC = T)
D.asv_eff_E_F01 <- aldex.effect(D.asv_clr_E_F01, include.sample.summary = T, 
                                useMC = T)
D.asv_eff_P_F01 <- aldex.effect(D.asv_clr_P_F01, include.sample.summary = T, 
                                useMC = T)
D.asv_eff_D_F01 <- aldex.effect(D.asv_clr_D_F01, include.sample.summary = T, 
                                useMC = T)
D.asv_eff_s_M00 <- aldex.effect(D.asv_clr_s_M00, include.sample.summary = T, 
                                useMC = T)
D.asv_eff_E_M00 <- aldex.effect(D.asv_clr_E_M00, include.sample.summary = T, 
                                useMC = T)
D.asv_eff_P_M00 <- aldex.effect(D.asv_clr_P_M00, include.sample.summary = T, 
                                useMC = T)
D.asv_eff_D_M00 <- aldex.effect(D.asv_clr_D_M00, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_s_F01 <- aldex.effect(D.phy_clr_s_F01, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_E_F01 <- aldex.effect(D.phy_clr_E_F01, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_P_F01 <- aldex.effect(D.phy_clr_P_F01, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_D_F01 <- aldex.effect(D.phy_clr_D_F01, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_s_M00 <- aldex.effect(D.phy_clr_s_M00, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_E_M00 <- aldex.effect(D.phy_clr_E_M00, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_P_M00 <- aldex.effect(D.phy_clr_P_M00, include.sample.summary = T, 
                                useMC = T)
D.phy_eff_D_M00 <- aldex.effect(D.phy_clr_D_M00, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_s_F01 <- aldex.effect(D.fam_clr_s_F01, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_E_F01 <- aldex.effect(D.fam_clr_E_F01, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_P_F01 <- aldex.effect(D.fam_clr_P_F01, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_D_F01 <- aldex.effect(D.fam_clr_D_F01, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_s_M00 <- aldex.effect(D.fam_clr_s_M00, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_E_M00 <- aldex.effect(D.fam_clr_E_M00, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_P_M00 <- aldex.effect(D.fam_clr_P_M00, include.sample.summary = T, 
                                useMC = T)
D.fam_eff_D_M00 <- aldex.effect(D.fam_clr_D_M00, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_s_F01 <- aldex.effect(D.lws_clr_s_F01, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_E_F01 <- aldex.effect(D.lws_clr_E_F01, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_P_F01 <- aldex.effect(D.lws_clr_P_F01, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_D_F01 <- aldex.effect(D.lws_clr_D_F01, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_s_M00 <- aldex.effect(D.lws_clr_s_M00, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_E_M00 <- aldex.effect(D.lws_clr_E_M00, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_P_M00 <- aldex.effect(D.lws_clr_P_M00, include.sample.summary = T, 
                                useMC = T)
D.lws_eff_D_M00 <- aldex.effect(D.lws_clr_D_M00, include.sample.summary = T, 
                                useMC = T)

### ************************************
### ^^^D - STEP  3 - format results ----
### ************************************

# ^^^ = code is mostly not commented in this step

# this process occurs as follows:
# (1) combine tests and effects data.frames
# (2)^^^add back in taxanomic info
# (3)^^^
# (4)^^^
# (5)^^^
# NOTE: note needed for Phylum or Family, just a name change

D.asv_rst_s_F01_0 <- data.frame(D.asv_tst_s_F01, D.asv_eff_s_F01)
D.asv_rst_E_F01_0 <- data.frame(D.asv_tst_E_F01, D.asv_eff_E_F01)
D.asv_rst_P_F01_0 <- data.frame(D.asv_tst_P_F01, D.asv_eff_P_F01)
D.asv_rst_D_F01_0 <- data.frame(D.asv_tst_D_F01, D.asv_eff_D_F01)
D.phy_rst_s_F01_0 <- data.frame(D.phy_tst_s_F01, D.phy_eff_s_F01)
D.phy_rst_E_F01_0 <- data.frame(D.phy_tst_E_F01, D.phy_eff_E_F01)
D.phy_rst_P_F01_0 <- data.frame(D.phy_tst_P_F01, D.phy_eff_P_F01)
D.phy_rst_D_F01_0 <- data.frame(D.phy_tst_D_F01, D.phy_eff_D_F01)
D.fam_rst_s_F01_0 <- data.frame(D.fam_tst_s_F01, D.fam_eff_s_F01)
D.fam_rst_E_F01_0 <- data.frame(D.fam_tst_E_F01, D.fam_eff_E_F01)
D.fam_rst_P_F01_0 <- data.frame(D.fam_tst_P_F01, D.fam_eff_P_F01)
D.fam_rst_D_F01_0 <- data.frame(D.fam_tst_D_F01, D.fam_eff_D_F01)
D.lws_rst_s_F01_0 <- data.frame(D.lws_tst_s_F01, D.lws_eff_s_F01)
D.lws_rst_E_F01_0 <- data.frame(D.lws_tst_E_F01, D.lws_eff_E_F01)
D.lws_rst_P_F01_0 <- data.frame(D.lws_tst_P_F01, D.lws_eff_P_F01)
D.lws_rst_D_F01_0 <- data.frame(D.lws_tst_D_F01, D.lws_eff_D_F01)
D.asv_rst_s_M00_0 <- data.frame(D.asv_tst_s_M00, D.asv_eff_s_M00)
D.asv_rst_E_M00_0 <- data.frame(D.asv_tst_E_M00, D.asv_eff_E_M00)
D.asv_rst_P_M00_0 <- data.frame(D.asv_tst_P_M00, D.asv_eff_P_M00)
D.asv_rst_D_M00_0 <- data.frame(D.asv_tst_D_M00, D.asv_eff_D_M00)
D.phy_rst_s_M00_0 <- data.frame(D.phy_tst_s_M00, D.phy_eff_s_M00)
D.phy_rst_E_M00_0 <- data.frame(D.phy_tst_E_M00, D.phy_eff_E_M00)
D.phy_rst_P_M00_0 <- data.frame(D.phy_tst_P_M00, D.phy_eff_P_M00)
D.phy_rst_D_M00_0 <- data.frame(D.phy_tst_D_M00, D.phy_eff_D_M00)
D.fam_rst_s_M00_0 <- data.frame(D.fam_tst_s_M00, D.fam_eff_s_M00)
D.fam_rst_E_M00_0 <- data.frame(D.fam_tst_E_M00, D.fam_eff_E_M00)
D.fam_rst_P_M00_0 <- data.frame(D.fam_tst_P_M00, D.fam_eff_P_M00)
D.fam_rst_D_M00_0 <- data.frame(D.fam_tst_D_M00, D.fam_eff_D_M00)
D.lws_rst_s_M00_0 <- data.frame(D.lws_tst_s_M00, D.lws_eff_s_M00)
D.lws_rst_E_M00_0 <- data.frame(D.lws_tst_E_M00, D.lws_eff_E_M00)
D.lws_rst_P_M00_0 <- data.frame(D.lws_tst_P_M00, D.lws_eff_P_M00)
D.lws_rst_D_M00_0 <- data.frame(D.lws_tst_D_M00, D.lws_eff_D_M00)

D.asv_rst_s_F01_1 <- D.asv_rst_s_F01_0
D.asv_rst_E_F01_1 <- D.asv_rst_E_F01_0 
D.asv_rst_P_F01_1 <- D.asv_rst_P_F01_0
D.asv_rst_D_F01_1 <- D.asv_rst_D_F01_0
D.phy_rst_s_F01_1 <- D.phy_rst_s_F01_0
D.phy_rst_E_F01_1 <- D.phy_rst_E_F01_0 
D.phy_rst_P_F01_1 <- D.phy_rst_P_F01_0
D.phy_rst_D_F01_1 <- D.phy_rst_D_F01_0
D.fam_rst_s_F01_1 <- D.fam_rst_s_F01_0
D.fam_rst_E_F01_1 <- D.fam_rst_E_F01_0 
D.fam_rst_P_F01_1 <- D.fam_rst_P_F01_0
D.fam_rst_D_F01_1 <- D.fam_rst_D_F01_0
D.lws_rst_s_F01_1 <- D.lws_rst_s_F01_0
D.lws_rst_E_F01_1 <- D.lws_rst_E_F01_0 
D.lws_rst_P_F01_1 <- D.lws_rst_P_F01_0
D.lws_rst_D_F01_1 <- D.lws_rst_D_F01_0
D.asv_rst_s_M00_1 <- D.asv_rst_s_M00_0
D.asv_rst_E_M00_1 <- D.asv_rst_E_M00_0 
D.asv_rst_P_M00_1 <- D.asv_rst_P_M00_0
D.asv_rst_D_M00_1 <- D.asv_rst_D_M00_0
D.phy_rst_s_M00_1 <- D.phy_rst_s_M00_0
D.phy_rst_E_M00_1 <- D.phy_rst_E_M00_0 
D.phy_rst_P_M00_1 <- D.phy_rst_P_M00_0
D.phy_rst_D_M00_1 <- D.phy_rst_D_M00_0
D.fam_rst_s_M00_1 <- D.fam_rst_s_M00_0
D.fam_rst_E_M00_1 <- D.fam_rst_E_M00_0 
D.fam_rst_P_M00_1 <- D.fam_rst_P_M00_0
D.fam_rst_D_M00_1 <- D.fam_rst_D_M00_0
D.lws_rst_s_M00_1 <- D.lws_rst_s_M00_0
D.lws_rst_E_M00_1 <- D.lws_rst_E_M00_0 
D.lws_rst_P_M00_1 <- D.lws_rst_P_M00_0
D.lws_rst_D_M00_1 <- D.lws_rst_D_M00_0

D.asv_rst_s_F01_1$FeatureID <- row.names(D.asv_rst_s_F01_1)
D.asv_rst_E_F01_1$FeatureID <- row.names(D.asv_rst_E_F01_1) 
D.asv_rst_P_F01_1$FeatureID <- row.names(D.asv_rst_P_F01_1)
D.asv_rst_D_F01_1$FeatureID <- row.names(D.asv_rst_D_F01_1)
D.phy_rst_s_F01_1$taxon <- row.names(D.phy_rst_s_F01_1)
D.phy_rst_E_F01_1$taxon <- row.names(D.phy_rst_E_F01_1) 
D.phy_rst_P_F01_1$taxon <- row.names(D.phy_rst_P_F01_1)
D.phy_rst_D_F01_1$taxon <- row.names(D.phy_rst_D_F01_1)
D.fam_rst_s_F01_1$taxon <- row.names(D.fam_rst_s_F01_1)
D.fam_rst_E_F01_1$taxon <- row.names(D.fam_rst_E_F01_1) 
D.fam_rst_P_F01_1$taxon <- row.names(D.fam_rst_P_F01_1)
D.fam_rst_D_F01_1$taxon <- row.names(D.fam_rst_D_F01_1)
D.lws_rst_s_F01_1$taxon <- row.names(D.lws_rst_s_F01_1)
D.lws_rst_E_F01_1$taxon <- row.names(D.lws_rst_E_F01_1) 
D.lws_rst_P_F01_1$taxon <- row.names(D.lws_rst_P_F01_1)
D.lws_rst_D_F01_1$taxon <- row.names(D.lws_rst_D_F01_1)
D.asv_rst_s_M00_1$FeatureID <- row.names(D.asv_rst_s_M00_1)
D.asv_rst_E_M00_1$FeatureID <- row.names(D.asv_rst_E_M00_1)
D.asv_rst_P_M00_1$FeatureID <- row.names(D.asv_rst_P_M00_1)
D.asv_rst_D_M00_1$FeatureID <- row.names(D.asv_rst_D_M00_1)
D.phy_rst_s_M00_1$taxon <- row.names(D.phy_rst_s_M00_1)
D.phy_rst_E_M00_1$taxon <- row.names(D.phy_rst_E_M00_1)
D.phy_rst_P_M00_1$taxon <- row.names(D.phy_rst_P_M00_1)
D.phy_rst_D_M00_1$taxon <- row.names(D.phy_rst_D_M00_1)
D.fam_rst_s_M00_1$taxon <- row.names(D.fam_rst_s_M00_1)
D.fam_rst_E_M00_1$taxon <- row.names(D.fam_rst_E_M00_1)
D.fam_rst_P_M00_1$taxon <- row.names(D.fam_rst_P_M00_1)
D.fam_rst_D_M00_1$taxon <- row.names(D.fam_rst_D_M00_1)
D.lws_rst_s_M00_1$taxon <- row.names(D.lws_rst_s_M00_1)
D.lws_rst_E_M00_1$taxon <- row.names(D.lws_rst_E_M00_1)
D.lws_rst_P_M00_1$taxon <- row.names(D.lws_rst_P_M00_1)
D.lws_rst_D_M00_1$taxon <- row.names(D.lws_rst_D_M00_1)

D.asv_rst_s_F01_2 <- dplyr::select(D.asv_rst_s_F01_1, FeatureID,
                                   dplyr::everything())
D.asv_rst_E_F01_2 <- dplyr::select(D.asv_rst_E_F01_1, FeatureID,
                                   dplyr::everything()) 
D.asv_rst_P_F01_2 <- dplyr::select(D.asv_rst_P_F01_1, FeatureID,
                                   dplyr::everything())
D.asv_rst_D_F01_2 <- dplyr::select(D.asv_rst_D_F01_1, FeatureID,
                                   dplyr::everything())
D.phy_rst_s_F01_2 <- dplyr::select(D.phy_rst_s_F01_1, taxon,
                                   dplyr::everything())
D.phy_rst_E_F01_2 <- dplyr::select(D.phy_rst_E_F01_1, taxon,
                                   dplyr::everything()) 
D.phy_rst_P_F01_2 <- dplyr::select(D.phy_rst_P_F01_1, taxon,
                                   dplyr::everything())
D.phy_rst_D_F01_2 <- dplyr::select(D.phy_rst_D_F01_1, taxon,
                                   dplyr::everything())
D.fam_rst_s_F01_2 <- dplyr::select(D.fam_rst_s_F01_1, taxon,
                                   dplyr::everything())
D.fam_rst_E_F01_2 <- dplyr::select(D.fam_rst_E_F01_1, taxon,
                                   dplyr::everything()) 
D.fam_rst_P_F01_2 <- dplyr::select(D.fam_rst_P_F01_1, taxon,
                                   dplyr::everything())
D.fam_rst_D_F01_2 <- dplyr::select(D.fam_rst_D_F01_1, taxon,
                                   dplyr::everything())
D.lws_rst_s_F01_2 <- dplyr::select(D.lws_rst_s_F01_1, taxon,
                                   dplyr::everything())
D.lws_rst_E_F01_2 <- dplyr::select(D.lws_rst_E_F01_1, taxon,
                                   dplyr::everything()) 
D.lws_rst_P_F01_2 <- dplyr::select(D.lws_rst_P_F01_1, taxon,
                                   dplyr::everything())
D.lws_rst_D_F01_2 <- dplyr::select(D.lws_rst_D_F01_1, taxon,
                                   dplyr::everything())
D.asv_rst_s_M00_2 <- dplyr::select(D.asv_rst_s_M00_1, FeatureID,
                                   dplyr::everything())
D.asv_rst_E_M00_2 <- dplyr::select(D.asv_rst_E_M00_1, FeatureID,
                                   dplyr::everything())
D.asv_rst_P_M00_2 <- dplyr::select(D.asv_rst_P_M00_1, FeatureID,
                                   dplyr::everything())
D.asv_rst_D_M00_2 <- dplyr::select(D.asv_rst_D_M00_1, FeatureID,
                                   dplyr::everything())
D.phy_rst_s_M00_2 <- dplyr::select(D.phy_rst_s_M00_1, taxon,
                                   dplyr::everything())
D.phy_rst_E_M00_2 <- dplyr::select(D.phy_rst_E_M00_1, taxon,
                                   dplyr::everything())
D.phy_rst_P_M00_2 <- dplyr::select(D.phy_rst_P_M00_1, taxon,
                                   dplyr::everything())
D.phy_rst_D_M00_2 <- dplyr::select(D.phy_rst_D_M00_1, taxon,
                                   dplyr::everything())
D.fam_rst_s_M00_2 <- dplyr::select(D.fam_rst_s_M00_1, taxon,
                                   dplyr::everything())
D.fam_rst_E_M00_2 <- dplyr::select(D.fam_rst_E_M00_1, taxon,
                                   dplyr::everything())
D.fam_rst_P_M00_2 <- dplyr::select(D.fam_rst_P_M00_1, taxon,
                                   dplyr::everything())
D.fam_rst_D_M00_2 <- dplyr::select(D.fam_rst_D_M00_1, taxon,
                                   dplyr::everything())
D.lws_rst_s_M00_2 <- dplyr::select(D.lws_rst_s_M00_1, taxon, 
                                   dplyr::everything())
D.lws_rst_E_M00_2 <- dplyr::select(D.lws_rst_E_M00_1, taxon,
                                   dplyr::everything())
D.lws_rst_P_M00_2 <- dplyr::select(D.lws_rst_P_M00_1, taxon,
                                   dplyr::everything())
D.lws_rst_D_M00_2 <- dplyr::select(D.lws_rst_D_M00_1, taxon,
                                   dplyr::everything())

D.asv_rst_s_F01_3 <- merge(x = D.asv_rst_s_F01_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.asv_rst_E_F01_3 <- merge(x = D.asv_rst_E_F01_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.asv_rst_P_F01_3 <- merge(x = D.asv_rst_P_F01_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.asv_rst_D_F01_3 <- merge(x = D.asv_rst_D_F01_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.phy_rst_s_F01_3 <- D.phy_rst_s_F01_2
D.phy_rst_E_F01_3 <- D.phy_rst_E_F01_2
D.phy_rst_P_F01_3 <- D.phy_rst_P_F01_2
D.phy_rst_D_F01_3 <- D.phy_rst_D_F01_2
D.fam_rst_s_F01_3 <- D.fam_rst_s_F01_2 
D.fam_rst_E_F01_3 <- D.fam_rst_E_F01_2 
D.fam_rst_P_F01_3 <- D.fam_rst_P_F01_2 
D.fam_rst_D_F01_3 <- D.fam_rst_D_F01_2
D.lws_rst_s_F01_3 <- merge(x = D.lws_rst_s_F01_2, y = D.slv_lws, 
                           sort = F, by = "taxon")
D.lws_rst_E_F01_3 <- merge(x = D.lws_rst_E_F01_2, y = D.slv_lws,  
                           sort = F, by = "taxon")
D.lws_rst_P_F01_3 <- merge(x = D.lws_rst_P_F01_2, y = D.slv_lws, 
                           sort = F, by = "taxon")
D.lws_rst_D_F01_3 <- merge(x = D.lws_rst_D_F01_2, y = D.slv_lws, 
                           sort = F, by = "taxon")
D.asv_rst_s_M00_3 <- merge(x = D.asv_rst_s_M00_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.asv_rst_E_M00_3 <- merge(x = D.asv_rst_E_M00_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.asv_rst_P_M00_3 <- merge(x = D.asv_rst_P_M00_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.asv_rst_D_M00_3 <- merge(x = D.asv_rst_D_M00_2, y = D.seq_tax, all.y = F, 
                           sort = F, by = "FeatureID")
D.phy_rst_s_M00_3 <- D.phy_rst_s_M00_2
D.phy_rst_E_M00_3 <- D.phy_rst_E_M00_2
D.phy_rst_P_M00_3 <- D.phy_rst_P_M00_2
D.phy_rst_D_M00_3 <- D.phy_rst_D_M00_2
D.fam_rst_s_M00_3 <- D.fam_rst_s_M00_2 
D.fam_rst_E_M00_3 <- D.fam_rst_E_M00_2  
D.fam_rst_P_M00_3 <- D.fam_rst_P_M00_2
D.fam_rst_D_M00_3 <- D.fam_rst_D_M00_2
D.lws_rst_s_M00_3 <- merge(x = D.lws_rst_s_M00_2, y = D.slv_lws, 
                           sort = F, by = "taxon")
D.lws_rst_E_M00_3 <- merge(x = D.lws_rst_E_M00_2, y = D.slv_lws, 
                           sort = F, by = "taxon")
D.lws_rst_P_M00_3 <- merge(x = D.lws_rst_P_M00_2, y = D.slv_lws, 
                           sort = F, by = "taxon")
D.lws_rst_D_M00_3 <- merge(x = D.lws_rst_D_M00_2, y = D.slv_lws,  
                           sort = F, by = "taxon")

# filter to retain Features with wi.eBH values less than or equal to 0.1
D.asv_sig_s_F01 <- dplyr::filter(D.asv_rst_s_F01_3, wi.eBH <= 0.1)
D.asv_sig_E_F01 <- dplyr::filter(D.asv_rst_E_F01_3, wi.eBH <= 0.1)
D.asv_sig_P_F01 <- dplyr::filter(D.asv_rst_P_F01_3, wi.eBH <= 0.1)
D.asv_sig_D_F01 <- dplyr::filter(D.asv_rst_D_F01_3, wi.eBH <= 0.1)
D.phy_sig_s_F01 <- dplyr::filter(D.phy_rst_s_F01_3, wi.eBH <= 0.1)
D.phy_sig_E_F01 <- dplyr::filter(D.phy_rst_E_F01_3, wi.eBH <= 0.1)
D.phy_sig_P_F01 <- dplyr::filter(D.phy_rst_P_F01_3, wi.eBH <= 0.1)
D.phy_sig_D_F01 <- dplyr::filter(D.phy_rst_D_F01_3, wi.eBH <= 0.1)
D.fam_sig_s_F01 <- dplyr::filter(D.fam_rst_s_F01_3, wi.eBH <= 0.1)
D.fam_sig_E_F01 <- dplyr::filter(D.fam_rst_E_F01_3, wi.eBH <= 0.1)
D.fam_sig_P_F01 <- dplyr::filter(D.fam_rst_P_F01_3, wi.eBH <= 0.1)
D.fam_sig_D_F01 <- dplyr::filter(D.fam_rst_D_F01_3, wi.eBH <= 0.1)
D.lws_sig_s_F01 <- dplyr::filter(D.lws_rst_s_F01_3, wi.eBH <= 0.1)
D.lws_sig_E_F01 <- dplyr::filter(D.lws_rst_E_F01_3, wi.eBH <= 0.1)
D.lws_sig_P_F01 <- dplyr::filter(D.lws_rst_P_F01_3, wi.eBH <= 0.1)
D.lws_sig_D_F01 <- dplyr::filter(D.lws_rst_D_F01_3, wi.eBH <= 0.1)
D.asv_sig_s_M00 <- dplyr::filter(D.asv_rst_s_M00_3, wi.eBH <= 0.1)
D.asv_sig_E_M00 <- dplyr::filter(D.asv_rst_E_M00_3, wi.eBH <= 0.1)
D.asv_sig_P_M00 <- dplyr::filter(D.asv_rst_P_M00_3, wi.eBH <= 0.1)
D.asv_sig_D_M00 <- dplyr::filter(D.asv_rst_D_M00_3, wi.eBH <= 0.1)
D.phy_sig_s_M00 <- dplyr::filter(D.phy_rst_s_M00_3, wi.eBH <= 0.1)
D.phy_sig_E_M00 <- dplyr::filter(D.phy_rst_E_M00_3, wi.eBH <= 0.1)
D.phy_sig_P_M00 <- dplyr::filter(D.phy_rst_P_M00_3, wi.eBH <= 0.1)
D.phy_sig_D_M00 <- dplyr::filter(D.phy_rst_D_M00_3, wi.eBH <= 0.1)
D.fam_sig_s_M00 <- dplyr::filter(D.fam_rst_s_M00_3, wi.eBH <= 0.1)
D.fam_sig_E_M00 <- dplyr::filter(D.fam_rst_E_M00_3, wi.eBH <= 0.1)
D.fam_sig_P_M00 <- dplyr::filter(D.fam_rst_P_M00_3, wi.eBH <= 0.1)
D.fam_sig_D_M00 <- dplyr::filter(D.fam_rst_D_M00_3, wi.eBH <= 0.1)
D.lws_sig_s_M00 <- dplyr::filter(D.lws_rst_s_M00_3, wi.eBH <= 0.1)
D.lws_sig_E_M00 <- dplyr::filter(D.lws_rst_E_M00_3, wi.eBH <= 0.1)
D.lws_sig_P_M00 <- dplyr::filter(D.lws_rst_P_M00_3, wi.eBH <= 0.1)
D.lws_sig_D_M00 <- dplyr::filter(D.lws_rst_D_M00_3, wi.eBH <= 0.1)

# count number of differentially abundant features at the specified threshold
# NOTE: numbers may change between runs (see vignette for ALDEx2)
print(nrow(D.asv_sig_s_F01)) # 0
print(nrow(D.asv_sig_E_F01)) # 4 OR 5
print(nrow(D.asv_sig_P_F01)) # 14
print(nrow(D.asv_sig_D_F01)) # 0
print(nrow(D.asv_sig_s_M00)) # 53 or 54
print(nrow(D.asv_sig_E_M00)) # 8
print(nrow(D.asv_sig_P_M00)) # 0
print(nrow(D.asv_sig_D_M00)) # 0
print(nrow(D.phy_sig_s_F01)) # 0
print(nrow(D.phy_sig_E_F01)) # 0
print(nrow(D.phy_sig_P_F01)) # 0
print(nrow(D.phy_sig_D_F01)) # 0
print(nrow(D.phy_sig_s_M00)) # 1
print(nrow(D.phy_sig_E_M00)) # 0
print(nrow(D.phy_sig_P_M00)) # 0
print(nrow(D.phy_sig_D_M00)) # 0
print(nrow(D.fam_sig_s_F01)) # 0
print(nrow(D.fam_sig_E_F01)) # 0
print(nrow(D.fam_sig_P_F01)) # 4
print(nrow(D.fam_sig_D_F01)) # 0
print(nrow(D.fam_sig_s_M00)) # 7 OR 8
print(nrow(D.fam_sig_E_M00)) # 2
print(nrow(D.fam_sig_P_M00)) # 0
print(nrow(D.fam_sig_D_M00)) # 0
print(nrow(D.lws_sig_s_F01)) # 0
print(nrow(D.lws_sig_E_F01)) # 4
print(nrow(D.lws_sig_P_F01)) # 12
print(nrow(D.lws_sig_D_F01)) # 0
print(nrow(D.lws_sig_s_M00)) # 39 or 40
print(nrow(D.lws_sig_E_M00)) # 8 or 9
print(nrow(D.lws_sig_P_M00)) # 0
print(nrow(D.lws_sig_D_M00)) # 2 to 7

# ^^^both the sigs data.frames and the full data.frames will be formattedsss
# ^^^sigs become Table X while the full data.frames will be filtered for plots

# for sigs, retain columns of interest
D.asv_frm_sig_E_F01_0 <- dplyr::select(D.asv_sig_E_F01, 
                                       rab.win.F01C, rab.win.F01R, diff.btw, 
                                       wi.ep, wi.eBH, int.ggs.lws.txn, 
                                       int.slv.lws.txn, RepSeq, int.ggs.tax, 
                                       int.slv.tax, FeatureID)
D.asv_frm_sig_P_F01_0 <- dplyr::select(D.asv_sig_P_F01, 
                                       rab.win.F01C, rab.win.F01R, diff.btw, 
                                       wi.ep, wi.eBH, int.ggs.lws.txn, 
                                       int.slv.lws.txn, RepSeq, int.ggs.tax, 
                                       int.slv.tax, FeatureID)
D.asv_frm_sig_s_M00_0 <- dplyr::select(D.asv_sig_s_M00, 
                                       rab.win.M01C, rab.win.M02R, diff.btw, 
                                       wi.ep, wi.eBH, int.ggs.lws.txn, 
                                       int.slv.lws.txn, RepSeq, int.ggs.tax, 
                                       int.slv.tax, FeatureID)
D.asv_frm_sig_E_M00_0 <- dplyr::select(D.asv_sig_E_M00, 
                                       rab.win.M01C, rab.win.M02R, diff.btw, 
                                       wi.ep, wi.eBH, int.ggs.lws.txn, 
                                       int.slv.lws.txn, RepSeq, int.ggs.tax, 
                                       int.slv.tax, FeatureID)
D.phy_frm_sig_s_M00_0 <- dplyr::select(D.phy_sig_s_M00,
                                       rab.win.M01C, rab.win.M02R, diff.btw,
                                       wi.ep, wi.eBH, taxon)
D.fam_frm_sig_P_F01_0 <- dplyr::select(D.fam_sig_P_F01,
                                       rab.win.F01C, rab.win.F01R, diff.btw,
                                       wi.ep, wi.eBH, taxon)
D.fam_frm_sig_s_M00_0 <- dplyr::select(D.fam_sig_s_M00,
                                       rab.win.M01C, rab.win.M02R, diff.btw,
                                       wi.ep, wi.eBH, taxon)
D.fam_frm_sig_E_M00_0 <- dplyr::select(D.fam_sig_E_M00,
                                       rab.win.M01C, rab.win.M02R, diff.btw,
                                       wi.ep, wi.eBH, taxon)
D.lws_frm_sig_E_F01_0 <- dplyr::select(D.lws_sig_E_F01,
                                       rab.win.F01C, rab.win.F01R, diff.btw,
                                       wi.ep, wi.eBH, taxon, level)
D.lws_frm_sig_P_F01_0 <- dplyr::select(D.lws_sig_P_F01,
                                       rab.win.F01C, rab.win.F01R, diff.btw,
                                       wi.ep, wi.eBH, taxon, level)
D.lws_frm_sig_s_M00_0 <- dplyr::select(D.lws_sig_s_M00,
                                       rab.win.M01C, rab.win.M02R, diff.btw,
                                       wi.ep, wi.eBH, taxon, level)
D.lws_frm_sig_E_M00_0 <- dplyr::select(D.lws_sig_E_M00,
                                       rab.win.M01C, rab.win.M02R, diff.btw,
                                       wi.ep, wi.eBH, taxon, level)
D.lws_frm_sig_D_M00_0 <- dplyr::select(D.lws_sig_D_M00,
                                       rab.win.M01C, rab.win.M02R, diff.btw,
                                       wi.ep, wi.eBH, taxon, level)

# determine direction of increase (i.e. which group was the ASV increased in?)
# this determination will use values in column 'diff.btw' as follows:
# comparison C vs R: + value = increased in R; - value = increased in C
# NOTE: this can be double checked/confirmed by looking at rows in 'rab.win.'

D.positive <- "RMC"
D.negative <- "CMC"

D.asv_frm_sig_E_F01_1 <- D.asv_frm_sig_E_F01_0
D.asv_frm_sig_P_F01_1 <- D.asv_frm_sig_P_F01_0
D.fam_frm_sig_P_F01_1 <- D.fam_frm_sig_P_F01_0
D.lws_frm_sig_E_F01_1 <- D.lws_frm_sig_E_F01_0
D.lws_frm_sig_P_F01_1 <- D.lws_frm_sig_P_F01_0
D.asv_frm_sig_s_M00_1 <- D.asv_frm_sig_s_M00_0
D.asv_frm_sig_E_M00_1 <- D.asv_frm_sig_E_M00_0
D.phy_frm_sig_s_M00_1 <- D.phy_frm_sig_s_M00_0
D.fam_frm_sig_s_M00_1 <- D.fam_frm_sig_s_M00_0
D.fam_frm_sig_E_M00_1 <- D.fam_frm_sig_E_M00_0
D.lws_frm_sig_s_M00_1 <- D.lws_frm_sig_s_M00_0
D.lws_frm_sig_E_M00_1 <- D.lws_frm_sig_E_M00_0
D.lws_frm_sig_D_M00_1 <- D.lws_frm_sig_D_M00_0
D.asv_rst_s_F01_4 <- D.asv_rst_s_F01_3
D.asv_rst_E_F01_4 <- D.asv_rst_E_F01_3
D.asv_rst_P_F01_4 <- D.asv_rst_P_F01_3
D.asv_rst_D_F01_4 <- D.asv_rst_D_F01_3
D.phy_rst_s_F01_4 <- D.phy_rst_s_F01_3
D.phy_rst_E_F01_4 <- D.phy_rst_E_F01_3
D.phy_rst_P_F01_4 <- D.phy_rst_P_F01_3
D.phy_rst_D_F01_4 <- D.phy_rst_D_F01_3
D.fam_rst_s_F01_4 <- D.fam_rst_s_F01_3 
D.fam_rst_E_F01_4 <- D.fam_rst_E_F01_3 
D.fam_rst_P_F01_4 <- D.fam_rst_P_F01_3 
D.fam_rst_D_F01_4 <- D.fam_rst_D_F01_3
D.lws_rst_s_F01_4 <- D.lws_rst_s_F01_3
D.lws_rst_E_F01_4 <- D.lws_rst_E_F01_3
D.lws_rst_P_F01_4 <- D.lws_rst_P_F01_3
D.lws_rst_D_F01_4 <- D.lws_rst_D_F01_3
D.asv_rst_s_M00_4 <- D.asv_rst_s_M00_3
D.asv_rst_E_M00_4 <- D.asv_rst_E_M00_3
D.asv_rst_P_M00_4 <- D.asv_rst_P_M00_3
D.asv_rst_D_M00_4 <- D.asv_rst_D_M00_3
D.phy_rst_s_M00_4 <- D.phy_rst_s_M00_3
D.phy_rst_E_M00_4 <- D.phy_rst_E_M00_3
D.phy_rst_P_M00_4 <- D.phy_rst_P_M00_3
D.phy_rst_D_M00_4 <- D.phy_rst_D_M00_3
D.fam_rst_s_M00_4 <- D.fam_rst_s_M00_3 
D.fam_rst_E_M00_4 <- D.fam_rst_E_M00_3 
D.fam_rst_P_M00_4 <- D.fam_rst_P_M00_3 
D.fam_rst_D_M00_4 <- D.fam_rst_D_M00_3
D.lws_rst_s_M00_4 <- D.lws_rst_s_M00_3
D.lws_rst_E_M00_4 <- D.lws_rst_E_M00_3
D.lws_rst_P_M00_4 <- D.lws_rst_P_M00_3
D.lws_rst_D_M00_4 <- D.lws_rst_D_M00_3

D.asv_frm_sig_E_F01_1$enrichment <- ifelse(
  D.asv_frm_sig_E_F01_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_frm_sig_P_F01_1$enrichment <- ifelse(
  D.asv_frm_sig_P_F01_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_frm_sig_P_F01_1$enrichment <- ifelse(
  D.fam_frm_sig_P_F01_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_frm_sig_E_F01_1$enrichment <- ifelse(
  D.lws_frm_sig_E_F01_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_frm_sig_P_F01_1$enrichment <- ifelse(
  D.lws_frm_sig_P_F01_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_frm_sig_s_M00_1$enrichment <- ifelse(
  D.asv_frm_sig_s_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_frm_sig_E_M00_1$enrichment <- ifelse(
  D.asv_frm_sig_E_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_frm_sig_s_M00_1$enrichment <- ifelse(
  D.phy_frm_sig_s_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_frm_sig_s_M00_1$enrichment <- ifelse(
  D.fam_frm_sig_s_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_frm_sig_E_M00_1$enrichment <- ifelse(
  D.fam_frm_sig_E_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_frm_sig_s_M00_1$enrichment <- ifelse(
  D.lws_frm_sig_s_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_frm_sig_E_M00_1$enrichment <- ifelse(
  D.lws_frm_sig_E_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_frm_sig_D_M00_1$enrichment <- ifelse(
  D.lws_frm_sig_D_M00_1[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_s_F01_4$enrichment <- ifelse(
  D.asv_rst_s_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_E_F01_4$enrichment <- ifelse(
  D.asv_rst_E_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_P_F01_4$enrichment <- ifelse(
  D.asv_rst_P_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_D_F01_4$enrichment <- ifelse(
  D.asv_rst_D_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_s_F01_4$enrichment <- ifelse(
  D.phy_rst_s_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_E_F01_4$enrichment <- ifelse(
  D.phy_rst_E_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_P_F01_4$enrichment <- ifelse(
  D.phy_rst_P_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_D_F01_4$enrichment <- ifelse(
  D.phy_rst_D_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_rst_s_F01_4$enrichment <- ifelse(
  D.fam_rst_s_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_rst_E_F01_4$enrichment <- ifelse(
  D.fam_rst_E_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative) 
D.fam_rst_P_F01_4$enrichment <- ifelse(
  D.fam_rst_P_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_rst_D_F01_4$enrichment <- ifelse(
  D.fam_rst_D_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_s_F01_4$enrichment <- ifelse(
  D.lws_rst_s_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_E_F01_4$enrichment <- ifelse(
  D.lws_rst_E_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_P_F01_4$enrichment <- ifelse(
  D.lws_rst_P_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_D_F01_4$enrichment <- ifelse(
  D.lws_rst_D_F01_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_s_M00_4$enrichment <- ifelse(
  D.asv_rst_s_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_E_M00_4$enrichment <- ifelse(
  D.asv_rst_E_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_P_M00_4$enrichment <- ifelse(
  D.asv_rst_P_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.asv_rst_D_M00_4$enrichment <- ifelse(
  D.asv_rst_D_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_s_M00_4$enrichment <- ifelse(
  D.phy_rst_s_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_E_M00_4$enrichment <- ifelse(
  D.phy_rst_E_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_P_M00_4$enrichment <- ifelse(
  D.phy_rst_P_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.phy_rst_D_M00_4$enrichment <- ifelse(
  D.phy_rst_D_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_rst_s_M00_4$enrichment <- ifelse(
  D.fam_rst_s_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_rst_E_M00_4$enrichment <- ifelse(
  D.fam_rst_E_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative) 
D.fam_rst_P_M00_4$enrichment <- ifelse(
  D.fam_rst_P_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.fam_rst_D_M00_4$enrichment <- ifelse(
  D.fam_rst_D_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_s_M00_4$enrichment <- ifelse(
  D.lws_rst_s_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_E_M00_4$enrichment <- ifelse(
  D.lws_rst_E_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_P_M00_4$enrichment <- ifelse(
  D.lws_rst_P_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)
D.lws_rst_D_M00_4$enrichment <- ifelse(
  D.lws_rst_D_M00_4[, "diff.btw"] > 0, yes = D.positive, no = D.negative)

# create new cols w/ info that will allow for smooth combination of results
D.asv_frm_sig_E_F01_2 <- D.asv_frm_sig_E_F01_1
D.asv_frm_sig_P_F01_2 <- D.asv_frm_sig_P_F01_1
D.fam_frm_sig_P_F01_2 <- D.fam_frm_sig_P_F01_1
D.lws_frm_sig_E_F01_2 <- D.lws_frm_sig_E_F01_1
D.lws_frm_sig_P_F01_2 <- D.lws_frm_sig_P_F01_1
D.asv_frm_sig_s_M00_2 <- D.asv_frm_sig_s_M00_1
D.asv_frm_sig_E_M00_2 <- D.asv_frm_sig_E_M00_1
D.phy_frm_sig_s_M00_2 <- D.phy_frm_sig_s_M00_1
D.fam_frm_sig_s_M00_2 <- D.fam_frm_sig_s_M00_1
D.fam_frm_sig_E_M00_2 <- D.fam_frm_sig_E_M00_1
D.lws_frm_sig_s_M00_2 <- D.lws_frm_sig_s_M00_1
D.lws_frm_sig_E_M00_2 <- D.lws_frm_sig_E_M00_1
D.lws_frm_sig_D_M00_2 <- D.lws_frm_sig_D_M00_1
D.asv_rst_s_F01_5 <- D.asv_rst_s_F01_4
D.asv_rst_E_F01_5 <- D.asv_rst_E_F01_4
D.asv_rst_P_F01_5 <- D.asv_rst_P_F01_4
D.asv_rst_D_F01_5 <- D.asv_rst_D_F01_4
D.phy_rst_s_F01_5 <- D.phy_rst_s_F01_4
D.phy_rst_E_F01_5 <- D.phy_rst_E_F01_4
D.phy_rst_P_F01_5 <- D.phy_rst_P_F01_4
D.phy_rst_D_F01_5 <- D.phy_rst_D_F01_4
D.fam_rst_s_F01_5 <- D.fam_rst_s_F01_4 
D.fam_rst_E_F01_5 <- D.fam_rst_E_F01_4 
D.fam_rst_P_F01_5 <- D.fam_rst_P_F01_4 
D.fam_rst_D_F01_5 <- D.fam_rst_D_F01_4
D.lws_rst_s_F01_5 <- D.lws_rst_s_F01_4
D.lws_rst_E_F01_5 <- D.lws_rst_E_F01_4
D.lws_rst_P_F01_5 <- D.lws_rst_P_F01_4
D.lws_rst_D_F01_5 <- D.lws_rst_D_F01_4
D.asv_rst_s_M00_5 <- D.asv_rst_s_M00_4
D.asv_rst_E_M00_5 <- D.asv_rst_E_M00_4
D.asv_rst_P_M00_5 <- D.asv_rst_P_M00_4
D.asv_rst_D_M00_5 <- D.asv_rst_D_M00_4
D.phy_rst_s_M00_5 <- D.phy_rst_s_M00_4
D.phy_rst_E_M00_5 <- D.phy_rst_E_M00_4
D.phy_rst_P_M00_5 <- D.phy_rst_P_M00_4
D.phy_rst_D_M00_5 <- D.phy_rst_D_M00_4
D.fam_rst_s_M00_5 <- D.fam_rst_s_M00_4 
D.fam_rst_E_M00_5 <- D.fam_rst_E_M00_4 
D.fam_rst_P_M00_5 <- D.fam_rst_P_M00_4 
D.fam_rst_D_M00_5 <- D.fam_rst_D_M00_4
D.lws_rst_s_M00_5 <- D.lws_rst_s_M00_4
D.lws_rst_E_M00_5 <- D.lws_rst_E_M00_4
D.lws_rst_P_M00_5 <- D.lws_rst_P_M00_4
D.lws_rst_D_M00_5 <- D.lws_rst_D_M00_4

D.comp_btw <- "RMC vs CMC"
D.asv_frm_sig_E_F01_2$comparison <- D.comp_btw
D.asv_frm_sig_P_F01_2$comparison <- D.comp_btw
D.fam_frm_sig_P_F01_2$comparison <- D.comp_btw
D.lws_frm_sig_E_F01_2$comparison <- D.comp_btw
D.lws_frm_sig_P_F01_2$comparison <- D.comp_btw
D.asv_frm_sig_s_M00_2$comparison <- D.comp_btw
D.asv_frm_sig_E_M00_2$comparison <- D.comp_btw
D.phy_frm_sig_s_M00_2$comparison <- D.comp_btw
D.fam_frm_sig_s_M00_2$comparison <- D.comp_btw
D.fam_frm_sig_E_M00_2$comparison <- D.comp_btw
D.lws_frm_sig_s_M00_2$comparison <- D.comp_btw
D.lws_frm_sig_E_M00_2$comparison <- D.comp_btw
D.lws_frm_sig_D_M00_2$comparison <- D.comp_btw
D.asv_rst_s_F01_5$comparison <- D.comp_btw
D.asv_rst_E_F01_5$comparison <- D.comp_btw
D.asv_rst_P_F01_5$comparison <- D.comp_btw
D.asv_rst_D_F01_5$comparison <- D.comp_btw
D.phy_rst_s_F01_5$comparison <- D.comp_btw
D.phy_rst_E_F01_5$comparison <- D.comp_btw
D.phy_rst_P_F01_5$comparison <- D.comp_btw
D.phy_rst_D_F01_5$comparison <- D.comp_btw
D.fam_rst_s_F01_5$comparison <- D.comp_btw 
D.fam_rst_E_F01_5$comparison <- D.comp_btw 
D.fam_rst_P_F01_5$comparison <- D.comp_btw 
D.fam_rst_D_F01_5$comparison <- D.comp_btw
D.lws_rst_s_F01_5$comparison <- D.comp_btw
D.lws_rst_E_F01_5$comparison <- D.comp_btw
D.lws_rst_P_F01_5$comparison <- D.comp_btw
D.lws_rst_D_F01_5$comparison <- D.comp_btw
D.asv_rst_s_M00_5$comparison <- D.comp_btw
D.asv_rst_E_M00_5$comparison <- D.comp_btw
D.asv_rst_P_M00_5$comparison <- D.comp_btw
D.asv_rst_D_M00_5$comparison <- D.comp_btw
D.phy_rst_s_M00_5$comparison <- D.comp_btw
D.phy_rst_E_M00_5$comparison <- D.comp_btw
D.phy_rst_P_M00_5$comparison <- D.comp_btw
D.phy_rst_D_M00_5$comparison <- D.comp_btw
D.fam_rst_s_M00_5$comparison <- D.comp_btw 
D.fam_rst_E_M00_5$comparison <- D.comp_btw 
D.fam_rst_P_M00_5$comparison <- D.comp_btw 
D.fam_rst_D_M00_5$comparison <- D.comp_btw
D.lws_rst_s_M00_5$comparison <- D.comp_btw
D.lws_rst_E_M00_5$comparison <- D.comp_btw
D.lws_rst_P_M00_5$comparison <- D.comp_btw
D.lws_rst_D_M00_5$comparison <- D.comp_btw

D.typ_ino <- "human.stool.inoculum"
D.typ_cec <- "cecum"
D.typ_prx <- "colon-prox"
D.typ_dtl <- "colon-dist"
D.asv_frm_sig_E_F01_2$type <- D.typ_cec
D.asv_frm_sig_P_F01_2$type <- D.typ_prx
D.fam_frm_sig_P_F01_2$type <- D.typ_prx
D.lws_frm_sig_E_F01_2$type <- D.typ_cec
D.lws_frm_sig_P_F01_2$type <- D.typ_prx
D.asv_frm_sig_s_M00_2$type <- D.typ_ino
D.asv_frm_sig_E_M00_2$type <- D.typ_cec
D.phy_frm_sig_s_M00_2$type <- D.typ_ino
D.fam_frm_sig_s_M00_2$type <- D.typ_ino
D.fam_frm_sig_E_M00_2$type <- D.typ_cec
D.lws_frm_sig_s_M00_2$type <- D.typ_ino
D.lws_frm_sig_E_M00_2$type <- D.typ_cec
D.lws_frm_sig_D_M00_2$type <- D.typ_dtl
D.asv_rst_s_F01_5$type <- D.typ_ino
D.asv_rst_E_F01_5$type <- D.typ_cec
D.asv_rst_P_F01_5$type <- D.typ_prx
D.asv_rst_D_F01_5$type <- D.typ_dtl
D.phy_rst_s_F01_5$type <- D.typ_ino
D.phy_rst_E_F01_5$type <- D.typ_cec
D.phy_rst_P_F01_5$type <- D.typ_prx
D.phy_rst_D_F01_5$type <- D.typ_dtl
D.fam_rst_s_F01_5$type <- D.typ_ino
D.fam_rst_E_F01_5$type <- D.typ_cec
D.fam_rst_P_F01_5$type <- D.typ_prx 
D.fam_rst_D_F01_5$type <- D.typ_dtl
D.lws_rst_s_F01_5$type <- D.typ_ino
D.lws_rst_E_F01_5$type <- D.typ_cec
D.lws_rst_P_F01_5$type <- D.typ_prx
D.lws_rst_D_F01_5$type <- D.typ_dtl
D.asv_rst_s_M00_5$type <- D.typ_ino
D.asv_rst_E_M00_5$type <- D.typ_cec
D.asv_rst_P_M00_5$type <- D.typ_prx
D.asv_rst_D_M00_5$type <- D.typ_dtl
D.phy_rst_s_M00_5$type <- D.typ_ino
D.phy_rst_E_M00_5$type <- D.typ_cec
D.phy_rst_P_M00_5$type <- D.typ_prx
D.phy_rst_D_M00_5$type <- D.typ_dtl
D.fam_rst_s_M00_5$type <- D.typ_ino
D.fam_rst_E_M00_5$type <- D.typ_cec
D.fam_rst_P_M00_5$type <- D.typ_prx 
D.fam_rst_D_M00_5$type <- D.typ_dtl
D.lws_rst_s_M00_5$type <- D.typ_ino
D.lws_rst_E_M00_5$type <- D.typ_cec
D.lws_rst_P_M00_5$type <- D.typ_prx
D.lws_rst_D_M00_5$type <- D.typ_dtl

D.consortium_F01 <- "derived from human female"
D.consortium_M00 <- "derived from human males"
D.asv_frm_sig_E_F01_2$group <- D.consortium_F01
D.asv_frm_sig_P_F01_2$group <- D.consortium_F01
D.fam_frm_sig_P_F01_2$group <- D.consortium_F01
D.lws_frm_sig_E_F01_2$group <- D.consortium_F01
D.lws_frm_sig_P_F01_2$group <- D.consortium_F01
D.asv_frm_sig_s_M00_2$group <- D.consortium_M00
D.asv_frm_sig_E_M00_2$group <- D.consortium_M00
D.phy_frm_sig_s_M00_2$group <- D.consortium_M00
D.fam_frm_sig_s_M00_2$group <- D.consortium_M00
D.fam_frm_sig_E_M00_2$group <- D.consortium_M00
D.lws_frm_sig_s_M00_2$group <- D.consortium_M00
D.lws_frm_sig_E_M00_2$group <- D.consortium_M00
D.lws_frm_sig_D_M00_2$group <- D.consortium_M00
D.asv_rst_s_F01_5$group <- D.consortium_F01
D.asv_rst_E_F01_5$group <- D.consortium_F01
D.asv_rst_P_F01_5$group <- D.consortium_F01
D.asv_rst_D_F01_5$group <- D.consortium_F01
D.phy_rst_s_F01_5$group <- D.consortium_F01
D.phy_rst_E_F01_5$group <- D.consortium_F01
D.phy_rst_P_F01_5$group <- D.consortium_F01
D.phy_rst_D_F01_5$group <- D.consortium_F01
D.fam_rst_s_F01_5$group <- D.consortium_F01
D.fam_rst_E_F01_5$group <- D.consortium_F01
D.fam_rst_P_F01_5$group <- D.consortium_F01 
D.fam_rst_D_F01_5$group <- D.consortium_F01
D.lws_rst_s_F01_5$group <- D.consortium_F01
D.lws_rst_E_F01_5$group <- D.consortium_F01
D.lws_rst_P_F01_5$group <- D.consortium_F01
D.lws_rst_D_F01_5$group <- D.consortium_F01
D.asv_rst_s_M00_5$group <- D.consortium_M00
D.asv_rst_E_M00_5$group <- D.consortium_M00
D.asv_rst_P_M00_5$group <- D.consortium_M00
D.asv_rst_D_M00_5$group <- D.consortium_M00
D.phy_rst_s_M00_5$group <- D.consortium_M00
D.phy_rst_E_M00_5$group <- D.consortium_M00
D.phy_rst_P_M00_5$group <- D.consortium_M00
D.phy_rst_D_M00_5$group <- D.consortium_M00
D.fam_rst_s_M00_5$group <- D.consortium_M00
D.fam_rst_E_M00_5$group <- D.consortium_M00
D.fam_rst_P_M00_5$group <- D.consortium_M00 
D.fam_rst_D_M00_5$group <- D.consortium_M00
D.lws_rst_s_M00_5$group <- D.consortium_M00
D.lws_rst_E_M00_5$group <- D.consortium_M00
D.lws_rst_P_M00_5$group <- D.consortium_M00
D.lws_rst_D_M00_5$group <- D.consortium_M00

# for rst (but not for sigs)
D.HumanDonorSex_F01 <- "female"
D.HumanDonorSex_M00 <- "male"
D.asv_rst_s_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.asv_rst_E_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.asv_rst_P_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.asv_rst_D_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.phy_rst_s_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.phy_rst_E_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.phy_rst_P_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.phy_rst_D_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.fam_rst_s_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.fam_rst_E_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.fam_rst_P_F01_5$HumanDonorSex <- D.HumanDonorSex_F01 
D.fam_rst_D_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.lws_rst_s_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.lws_rst_E_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.lws_rst_P_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.lws_rst_D_F01_5$HumanDonorSex <- D.HumanDonorSex_F01
D.asv_rst_s_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.asv_rst_E_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.asv_rst_P_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.asv_rst_D_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.phy_rst_s_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.phy_rst_E_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.phy_rst_P_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.phy_rst_D_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.fam_rst_s_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.fam_rst_E_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.fam_rst_P_M00_5$HumanDonorSex <- D.HumanDonorSex_M00 
D.fam_rst_D_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.lws_rst_s_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.lws_rst_E_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.lws_rst_P_M00_5$HumanDonorSex <- D.HumanDonorSex_M00
D.lws_rst_D_M00_5$HumanDonorSex <- D.HumanDonorSex_M00

# phy/fam/lws specific:
D.clp_phy <- "Phylum Level"
D.clp_fam <- "Family Level"
D.clp_lws <- "Lowest Assignment"
D.fam_frm_sig_P_F01_2$collapse <- D.clp_fam
D.lws_frm_sig_E_F01_2$collapse <- D.clp_lws
D.lws_frm_sig_P_F01_2$collapse <- D.clp_lws
D.phy_frm_sig_s_M00_2$collapse <- D.clp_phy
D.fam_frm_sig_s_M00_2$collapse <- D.clp_fam
D.fam_frm_sig_E_M00_2$collapse <- D.clp_fam
D.lws_frm_sig_s_M00_2$collapse <- D.clp_lws
D.lws_frm_sig_E_M00_2$collapse <- D.clp_lws
D.lws_frm_sig_D_M00_2$collapse <- D.clp_lws
D.phy_rst_s_F01_5$collapse <- D.clp_phy
D.phy_rst_E_F01_5$collapse <- D.clp_phy
D.phy_rst_P_F01_5$collapse <- D.clp_phy
D.phy_rst_D_F01_5$collapse <- D.clp_phy
D.fam_rst_s_F01_5$collapse <- D.clp_fam
D.fam_rst_E_F01_5$collapse <- D.clp_fam
D.fam_rst_P_F01_5$collapse <- D.clp_fam
D.fam_rst_D_F01_5$collapse <- D.clp_fam
D.lws_rst_s_F01_5$collapse <- D.clp_lws
D.lws_rst_E_F01_5$collapse <- D.clp_lws
D.lws_rst_P_F01_5$collapse <- D.clp_lws
D.lws_rst_D_F01_5$collapse <- D.clp_lws
D.phy_rst_s_M00_5$collapse <- D.clp_phy
D.phy_rst_E_M00_5$collapse <- D.clp_phy
D.phy_rst_P_M00_5$collapse <- D.clp_phy
D.phy_rst_D_M00_5$collapse <- D.clp_phy
D.fam_rst_s_M00_5$collapse <- D.clp_fam
D.fam_rst_E_M00_5$collapse <- D.clp_fam
D.fam_rst_P_M00_5$collapse <- D.clp_fam
D.fam_rst_D_M00_5$collapse <- D.clp_fam
D.lws_rst_s_M00_5$collapse <- D.clp_lws
D.lws_rst_E_M00_5$collapse <- D.clp_lws
D.lws_rst_P_M00_5$collapse <- D.clp_lws
D.lws_rst_D_M00_5$collapse <- D.clp_lws

# reorder, rename, and drop unneeded columns
D.asv_frm_sig_E_F01_3 <- dplyr::select(D.asv_frm_sig_E_F01_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = int.slv.lws.txn, 
                                       clr_abund_CMC = rab.win.F01C,
                                       clr_abund_RMC = rab.win.F01R,
                                       ASV_sequence = RepSeq,
                                       taxon_Greengenes = int.ggs.lws.txn, 
                                       lineage_Greengenes = int.ggs.tax, 
                                       lineage_SILVA = int.slv.tax,
                                       FeatureID)
D.asv_frm_sig_P_F01_3 <- dplyr::select(D.asv_frm_sig_P_F01_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = int.slv.lws.txn, 
                                       clr_abund_CMC = rab.win.F01C,
                                       clr_abund_RMC = rab.win.F01R,
                                       ASV_sequence = RepSeq,
                                       taxon_Greengenes = int.ggs.lws.txn, 
                                       lineage_Greengenes = int.ggs.tax, 
                                       lineage_SILVA = int.slv.tax,
                                       FeatureID)
D.fam_frm_sig_P_F01_3 <- dplyr::select(D.fam_frm_sig_P_F01_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.F01C,
                                       clr_abund_RMC = rab.win.F01R)
D.lws_frm_sig_E_F01_3 <- dplyr::select(D.lws_frm_sig_E_F01_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.F01C,
                                       clr_abund_RMC = rab.win.F01R)
D.lws_frm_sig_P_F01_3 <- dplyr::select(D.lws_frm_sig_P_F01_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.F01C,
                                       clr_abund_RMC = rab.win.F01R)
D.asv_frm_sig_s_M00_3 <- dplyr::select(D.asv_frm_sig_s_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = int.slv.lws.txn, 
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R,
                                       ASV_sequence = RepSeq,
                                       taxon_Greengenes = int.ggs.lws.txn, 
                                       lineage_Greengenes = int.ggs.tax, 
                                       lineage_SILVA = int.slv.tax,
                                       FeatureID) 
D.asv_frm_sig_E_M00_3 <- dplyr::select(D.asv_frm_sig_E_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = int.slv.lws.txn, 
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R,
                                       ASV_sequence = RepSeq,
                                       taxon_Greengenes = int.ggs.lws.txn, 
                                       lineage_Greengenes = int.ggs.tax, 
                                       lineage_SILVA = int.slv.tax,
                                       FeatureID) 
D.phy_frm_sig_s_M00_3 <- dplyr::select(D.phy_frm_sig_s_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R)
D.fam_frm_sig_s_M00_3 <- dplyr::select(D.fam_frm_sig_s_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R)
D.fam_frm_sig_E_M00_3 <- dplyr::select(D.fam_frm_sig_E_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R)
D.lws_frm_sig_s_M00_3 <- dplyr::select(D.lws_frm_sig_s_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R)
D.lws_frm_sig_E_M00_3 <- dplyr::select(D.lws_frm_sig_E_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R)
D.lws_frm_sig_D_M00_3 <- dplyr::select(D.lws_frm_sig_D_M00_2, group, type, 
                                       comparison, enrichment, 
                                       P = wi.ep, BH_P = wi.eBH, 
                                       log2_fold_diff = diff.btw,
                                       taxon_SILVA = taxon, 
                                       collapse,
                                       clr_abund_CMC = rab.win.M01C,
                                       clr_abund_RMC = rab.win.M02R)
D.asv_rst_s_F01 <- dplyr::select(D.asv_rst_s_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.asv_rst_E_F01 <- dplyr::select(D.asv_rst_E_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.asv_rst_P_F01 <- dplyr::select(D.asv_rst_P_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.asv_rst_D_F01 <- dplyr::select(D.asv_rst_D_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.phy_rst_s_F01 <- dplyr::select(D.phy_rst_s_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.phy_rst_E_F01 <- dplyr::select(D.phy_rst_E_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.phy_rst_P_F01 <- dplyr::select(D.phy_rst_P_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.phy_rst_D_F01 <- dplyr::select(D.phy_rst_D_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.fam_rst_s_F01 <- dplyr::select(D.fam_rst_s_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.fam_rst_E_F01 <- dplyr::select(D.fam_rst_E_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.fam_rst_P_F01 <- dplyr::select(D.fam_rst_P_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.fam_rst_D_F01 <- dplyr::select(D.fam_rst_D_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.lws_rst_s_F01 <- dplyr::select(D.lws_rst_s_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.lws_rst_E_F01 <- dplyr::select(D.lws_rst_E_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.lws_rst_P_F01 <- dplyr::select(D.lws_rst_P_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.lws_rst_D_F01 <- dplyr::select(D.lws_rst_D_F01_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.F01C,
                                 clr_abund_RMC = rab.win.F01R,
                                 comparison)
D.asv_rst_s_M00 <- dplyr::select(D.asv_rst_s_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.asv_rst_E_M00 <- dplyr::select(D.asv_rst_E_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.asv_rst_P_M00 <- dplyr::select(D.asv_rst_P_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.asv_rst_D_M00 <- dplyr::select(D.asv_rst_D_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = int.slv.lws.txn, 
                                 ASV_sequence = RepSeq,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 taxon_Greengenes = int.ggs.lws.txn, 
                                 lineage_Greengenes = int.ggs.tax, 
                                 lineage_SILVA = int.slv.tax,
                                 FeatureID, comparison)
D.phy_rst_s_M00 <- dplyr::select(D.phy_rst_s_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.phy_rst_E_M00 <- dplyr::select(D.phy_rst_E_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.phy_rst_P_M00 <- dplyr::select(D.phy_rst_P_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.phy_rst_D_M00 <- dplyr::select(D.phy_rst_D_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.fam_rst_s_M00 <- dplyr::select(D.fam_rst_s_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.fam_rst_E_M00 <- dplyr::select(D.fam_rst_E_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.fam_rst_P_M00 <- dplyr::select(D.fam_rst_P_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.fam_rst_D_M00 <- dplyr::select(D.fam_rst_D_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.lws_rst_s_M00 <- dplyr::select(D.lws_rst_s_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.lws_rst_E_M00 <- dplyr::select(D.lws_rst_E_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.lws_rst_P_M00 <- dplyr::select(D.lws_rst_P_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)
D.lws_rst_D_M00 <- dplyr::select(D.lws_rst_D_M00_5, 
                                 group, HumanDonorSex, type, enrichment, 
                                 P = wi.ep, 
                                 BH_P = wi.eBH, 
                                 log2_fold_diff = diff.btw,
                                 taxon_SILVA = taxon, 
                                 collapse,
                                 clr_abund_CMC = rab.win.M01C,
                                 clr_abund_RMC = rab.win.M02R,
                                 comparison)

# for sigs data.frames:
# (1) sort by direction of increase (enrichment)
# (2) combine above formatted data.frames
# (3) round all numeric values to the three decimal places

D.asv_frm_sig_E_F01 <- D.asv_frm_sig_E_F01_3[
  order(D.asv_frm_sig_E_F01_3$enrichment, decreasing = T), ]
D.asv_frm_sig_P_F01 <- D.asv_frm_sig_P_F01_3[
  order(D.asv_frm_sig_P_F01_3$enrichment, decreasing = T), ]
D.fam_frm_sig_P_F01 <- D.fam_frm_sig_P_F01_3[
  order(D.fam_frm_sig_P_F01_3$enrichment, decreasing = T), ]
D.lws_frm_sig_E_F01 <- D.lws_frm_sig_E_F01_3[
  order(D.lws_frm_sig_E_F01_3$enrichment, decreasing = T), ]
D.lws_frm_sig_P_F01 <- D.lws_frm_sig_P_F01_3[
  order(D.lws_frm_sig_P_F01_3$enrichment, decreasing = T), ]
D.asv_frm_sig_s_M00 <- D.asv_frm_sig_s_M00_3[
  order(D.asv_frm_sig_s_M00_3$enrichment, decreasing = T), ]
D.asv_frm_sig_E_M00 <- D.asv_frm_sig_E_M00_3[
  order(D.asv_frm_sig_E_M00_3$enrichment, decreasing = T), ]
D.phy_frm_sig_s_M00 <- D.phy_frm_sig_s_M00_3[
  order(D.phy_frm_sig_s_M00_3$enrichment, decreasing = T), ]
D.fam_frm_sig_s_M00 <- D.fam_frm_sig_s_M00_3[
  order(D.fam_frm_sig_s_M00_3$enrichment, decreasing = T), ]
D.fam_frm_sig_E_M00 <- D.fam_frm_sig_E_M00_3[
  order(D.fam_frm_sig_E_M00_3$enrichment, decreasing = T), ]
D.lws_frm_sig_s_M00 <- D.lws_frm_sig_s_M00_3[
  order(D.lws_frm_sig_s_M00_3$enrichment, decreasing = T), ]
D.lws_frm_sig_E_M00 <- D.lws_frm_sig_E_M00_3[
  order(D.lws_frm_sig_E_M00_3$enrichment, decreasing = T), ]
D.lws_frm_sig_D_M00 <- D.lws_frm_sig_D_M00_3[
  order(D.lws_frm_sig_D_M00_3$enrichment, decreasing = T), ]

D.alx_sig_asv_0 <- rbind(D.asv_frm_sig_E_F01, D.asv_frm_sig_P_F01, 
                         D.asv_frm_sig_s_M00, D.asv_frm_sig_E_M00)
D.alx_sig_tax_0 <- rbind(D.fam_frm_sig_P_F01, 
                         D.lws_frm_sig_E_F01, 
                         D.lws_frm_sig_P_F01, 
                         D.phy_frm_sig_s_M00, 
                         D.fam_frm_sig_s_M00, 
                         D.fam_frm_sig_E_M00, 
                         D.lws_frm_sig_s_M00, 
                         D.lws_frm_sig_E_M00, 
                         D.lws_frm_sig_D_M00)

D.alx_sig_asv <- dplyr::mutate_if(D.alx_sig_asv_0, is.numeric, round, 
                                  digits = 3)
D.alx_sig_tax <- dplyr::mutate_if(D.alx_sig_tax_0, is.numeric, round, 
                                  digits = 3)

# for full data.frames:
# (1) combine phy/fam/lws by sample type
# (2) combine asv and tax by sex of human donor
# (3) combine all together

D.alx_all_tax_s_F01 <- rbind(D.phy_rst_s_F01, D.fam_rst_s_F01, D.lws_rst_s_F01)
D.alx_all_tax_E_F01 <- rbind(D.phy_rst_E_F01, D.fam_rst_E_F01, D.lws_rst_E_F01)
D.alx_all_tax_P_F01 <- rbind(D.phy_rst_P_F01, D.fam_rst_P_F01, D.lws_rst_P_F01)
D.alx_all_tax_D_F01 <- rbind(D.phy_rst_D_F01, D.fam_rst_D_F01, D.lws_rst_D_F01)
D.alx_all_tax_s_M00 <- rbind(D.phy_rst_s_M00, D.fam_rst_s_M00, D.lws_rst_s_M00)
D.alx_all_tax_E_M00 <- rbind(D.phy_rst_E_M00, D.fam_rst_E_M00, D.lws_rst_E_M00)
D.alx_all_tax_P_M00 <- rbind(D.phy_rst_P_M00, D.fam_rst_P_M00, D.lws_rst_P_M00)
D.alx_all_tax_D_M00 <- rbind(D.phy_rst_D_M00, D.fam_rst_D_M00, D.lws_rst_D_M00)

D.alx_all_asv_F01 <- rbind(D.asv_rst_s_F01, D.asv_rst_E_F01, D.asv_rst_P_F01, 
                           D.asv_rst_D_F01)
D.alx_all_asv_M00 <- rbind(D.asv_rst_s_M00, D.asv_rst_E_M00, D.asv_rst_P_M00, 
                           D.asv_rst_D_M00)
D.alx_all_tax_F01 <- rbind(D.alx_all_tax_s_F01, D.alx_all_tax_E_F01, 
                           D.alx_all_tax_P_F01, D.alx_all_tax_D_F01)
D.alx_all_tax_M00 <- rbind(D.alx_all_tax_s_M00, D.alx_all_tax_E_M00, 
                           D.alx_all_tax_P_M00, D.alx_all_tax_D_M00)

D.alx_all_asv <- rbind(D.alx_all_asv_F01, D.alx_all_asv_M00)
D.alx_all_tax <- rbind(D.alx_all_tax_F01, D.alx_all_tax_M00)

### ************************************
### D - WRITE OUTPUTS ----
### ************************************

D.ofv_alx_sig_asv <- "TransFaunation/vault/table_alx_asv_sig.txt"
D.ofv_alx_sig_tax <- "TransFaunation/vault/table_alx_tax_sig.txt"
write.table(sep = "\t", row.names = F, x = D.alx_sig_asv,
            file = D.ofv_alx_sig_asv)
write.table(sep = "\t", row.names = F, x = D.alx_sig_tax,
            file = D.ofv_alx_sig_tax)

D.ofv_alx_all_asv <- "TransFaunation/vault/table_alx_asv_all.txt"
D.ofv_alx_all_tax <- "TransFaunation/vault/table_alx_tax_all.txt"
write.table(sep = "\t", row.names = F, x = D.alx_all_asv,
            file = D.ofv_alx_all_asv)
write.table(sep = "\t", row.names = F, x = D.alx_all_tax,
            file = D.ofv_alx_all_tax)

# save workspace
D.obj <- ls(pattern = "D.")
D.lst <- c(D.obj[grep(pattern = "D.", x = D.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, D.obj_from_AS, COSMOS)
save(list = D.lst, file = D.ofv_wksp)

### ************************************
### BEGIN Section F ----
### ************************************

# note for kdp: taken from TransFaunation/section_F.R

# NOTE: code is essentially completely uncommented
# NOTE: for final version of this code; look at inputs and format obj from secD

### preface ----

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
load(file = ofv_COSMOS_wksp)

F.ofv_wksp <- paste(sep = "", path_vault, "/WS_", "SectionF.RData")

# relative paths from wd for inputs
F.ifp_dat_F01 <- "TransFaunation/asv_log2_plots/alx_log2_F01.txt"
F.ifp_dat_M00 <- "TransFaunation/asv_log2_plots/alx_log2_M00.txt"
F.ifp_cons <- "TransFaunation/asv_log2_plots/alx_log2_cons.txt"

F.dat_F01_0 <- read.table(F.ifp_dat_F01, sep = "\t", header = T, as.is = T, 
                          stringsAsFactors = F)
F.dat_M00_0 <- read.table(F.ifp_dat_M00, sep = "\t", header = T, as.is = T, 
                          stringsAsFactors = F)

F.cons_0 <- read.table(F.ifp_cons, sep = "\t", header = T, as.is = T, 
                       stringsAsFactors = F)

### STEP 1 ---- 

# convert col x to character and define breaks and labels for x axis
F.dat_F01_1 <- F.dat_F01_0
F.dat_F01_1$x <- as.character(F.dat_F01_1$x)

F.xbrk_F01 <- F.dat_F01_1$x
F.xlab_F01 <- F.dat_F01_1$taxon

F.dat_M00_1 <- F.dat_M00_0
F.dat_M00_1$x <- as.character(F.dat_M00_1$x)
F.xbrk_M00 <- F.dat_M00_1$x
F.xlab_M00 <- F.dat_M00_1$taxon

F.cons_1 <- F.cons_0
F.cons_1$x <- as.character(F.cons_1$x)
F.xbrk_cons <- F.cons_1$x
F.xlab_cons <- F.cons_1$type

### parameters ---- 

# y-axis limits, breaks, labels
F.ylim <- c(-12, 12)
F.ybrk <- c(-10, -5, 0, 5, 10)
F.ylab <- c("-10", " -5", "0", " 5", "10")
F.yexp <- c(0, 0)

F.yttl <- "log2 fold difference"

# x-axis expand
F.xexp <- c(0, 0.6)

# bar sizing
F.wid <- 0.66

# x axis ticks
F.xtck_F01 <- c(rep("solid", times = 11), "blank",
                rep("solid", times = 02), "blank")

F.xtck_M00 <- c(rep("solid", times = 03), "blank",
                rep("solid", times = 10), "blank")


F.xtck_cons <- c(rep("solid", times = 4), "blank",
                 rep("solid", times = 4), "blank")

# vector for medium colors by consortium group
F.hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2])
F.hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2])
F.hex_cons <- c(hex_F01_cmc[2], hex_M00_cmc[2], hex_F01_rmc[2], hex_M00_rmc[2])

# legend breaks and labels
F.lbrk <- c("CMC", "RMC")
F.llab_F01 <- c("CMC-f", "RMC-f")
F.llab_M00 <- c("CMC-m", "RMC-m")

F.lbrk_cons <- c("CMC-f", "RMC-f", "CMC-m", "RMC-m")
F.llab_cons <- c("CMC-f", "RMC-f", "CMC-m", "RMC-m")

# sizing
F.sze_sig <- 2.12 # size for significance text
F.sze_atl <- 06 # size for axis titles (x and y)
F.sze_atx <- 06
F.sze_ltx <- 05 # size for legend text
F.sze_vlne <- 0.5
F.vlne_sclr <- 0.2

# custom theme parameters
F.log2 <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]),
  axis.text.x = element_text(size = F.sze_atx),
  axis.text.y = element_text(face = "italic", size = F.sze_atx),
  axis.title.x = element_text(size = F.sze_atl),
  axis.title.y = element_blank(),
  legend.background = element_rect(fill = NA),
  legend.key = element_rect(fill = NA, color = NA),
  legend.box.margin = margin(-3, 0, -6, 0, "mm"),
  legend.key.size = unit(2, "mm"),
  legend.title = element_blank(),
  legend.text = element_text(
    margin = margin(0.5, 0, 0.5, -2, "mm"), size = F.sze_ltx, 
    color = greydient[1], hjust = 0),
  panel.ontop = T,
  panel.background = element_rect(fill = NA))


F.log2_cons <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]),
  axis.text.x = element_text(size = F.sze_atx),
  axis.text.y = element_text(size = F.sze_atx),
  axis.title.x = element_text(size = F.sze_atl),
  axis.title.y = element_blank(),
  legend.background = element_rect(fill = NA),
  legend.key = element_rect(fill = NA, color = NA),
  legend.box.margin = margin(-3, 0, -6, 0, "mm"),
  legend.key.size = unit(2, "mm"),
  legend.title = element_blank(),
  legend.text = element_text(
    margin = margin(0.5, 0, 0.5, -2, "mm"), size = F.sze_ltx, 
    color = greydient[1], hjust = 0),
  panel.ontop = T,
  panel.background = element_rect(fill = NA))

# plotting  (non conserved) ----

F.gpr_log2_F01 <- ggbarplot(data = F.dat_F01_1, width = F.wid, 
                            font.family = fnt,
                            x = "x", y = "log2fd", fill = "enrichment",
                            color = greydient[8], orientation = "horiz") +
  labs(y = F.yttl) +
  geom_vline(xintercept = 12 + F.vlne_sclr, size = F.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  scale_x_discrete(breaks = F.xbrk_F01, labels = F.xlab_F01) +
  scale_y_continuous(limits = F.ylim, breaks = F.ybrk, labels = F.ylab,
                     expand = F.yexp) +
  scale_fill_manual(values = F.hex_F01, breaks = F.lbrk, labels = F.llab_F01) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = F.xtck_F01)) +
  F.log2

F.gpr_log2_M00 <- ggbarplot(data = F.dat_M00_1, width = F.wid,
                            font.family = fnt,
                            x = "x", y = "log2fd", fill = "enrichment",
                            color = greydient[8], orientation = "horiz") +
  geom_vline(xintercept = 04 + F.vlne_sclr, size = F.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  labs(y = F.yttl) +
  scale_x_discrete(breaks = F.xbrk_M00, labels = F.xlab_M00) +
  scale_y_continuous(limits = F.ylim, breaks = F.ybrk, labels = F.ylab,
                     expand = F.yexp) +
  scale_fill_manual(values = F.hex_M00, breaks = F.lbrk, labels = F.llab_M00) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = F.xtck_M00)) +
  F.log2


# plotting  (conserved) ----

F.ysclr <- 0.5

F.gpr_log2_cons <- ggbarplot(data = F.cons_1, width = F.wid,
                             font.family = fnt,
                             x = "x", y = "log2fd", fill = "enrichment.group",
                             color = greydient[8], orientation = "horiz") +
  labs(y = F.yttl) +
  geom_vline(xintercept = 5 + F.vlne_sclr, size = F.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  geom_text(aes(x = x, label = sig_BH_P,
                y = ifelse(log2fd > 0, yes = log2fd + F.ysclr,
                           no = log2fd - F.ysclr)),
            size = F.sze_sig, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = F.xbrk_cons, labels = F.xlab_cons) +
  scale_y_continuous(limits = F.ylim, breaks = F.ybrk, labels = F.ylab,
                     expand = F.yexp) +
  scale_fill_manual(values = F.hex_cons, breaks = F.lbrk_cons, 
                    labels = F.llab_cons) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = F.xtck_cons)) +
  F.log2_cons

### plotting arrange ----

F.gga_log2 <- ggarrange(F.gpr_log2_F01, F.gpr_log2_M00, F.gpr_log2_cons,
                        labels = c("A", "B", "C"), font.label = pan_fnt, 
                        ncol = 1, nrow = 3, align = "hv",
                        heights = c(1.5, 1.5, 1))

### save ----

F.ofv_plot_F01 <- "TransFaunation/vault/plot_alx_log2_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 55,
       filename = F.ofv_plot_F01, plot = F.gpr_log2_F01)
F.ofv_plot_M00 <- "TransFaunation/vault/plot_alx_log2_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 55,
       filename = F.ofv_plot_M00, plot = F.gpr_log2_M00)

F.ofv_gpr_log2_cons <- "TransFaunation/vault/plot_alx_log2_cons.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 37,
       filename = F.ofv_gpr_log2_cons, plot = F.gpr_log2_cons)

F.ofv_plot <- "TransFaunation/vault/plot_alx_log2.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 145,
       filename = F.ofv_plot, plot = F.gga_log2)

# save workspace
F.obj <- ls(pattern = "F.")
F.lst <- c(F.obj[grep(pattern = "F.", x = F.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS)
save(list = F.lst, file = F.ofv_wksp)

### ************************************
### BEGIN Section J ----
### ************************************

# note for kdp: taken from TransFaunation/section_B.R

#### preface ----

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
load(file = ofv_COSMOS_wksp)

J.ofv_wksp <- paste(sep = "", path_vault, "/WS_", "SectionJ.RData")

# relative paths from wd for inputs
J.ifp_mtb_prox <- "TransFaunation/mtbs/cons_mtb_prox.txt"
J.ifp_mtb_dist <- "TransFaunation/mtbs/cons_mtb_dist.txt"
J.ifp_mtb_plas <- "TransFaunation/mtbs/cons_mtb_plas.txt"
J.ifp_mtb_cros <- "TransFaunation/mtbs/cons_mtb_cros.txt"

# read in the data
J.mtb_prox_0 <- read.table(J.ifp_mtb_prox, sep = "\t", header = T, as.is = T, 
                           stringsAsFactors = F, check.names = F, quote = "",
                           comment.char = "")
J.mtb_dist_0 <- read.table(J.ifp_mtb_dist, sep = "\t", header = T, as.is = T, 
                           stringsAsFactors = F, check.names = F, quote = "",
                           comment.char = "")
J.mtb_plas_0 <- read.table(J.ifp_mtb_plas, sep = "\t", header = T, as.is = T, 
                           stringsAsFactors = F, check.names = F, quote = "",
                           comment.char = "")
J.mtb_cros_0 <- read.table(J.ifp_mtb_cros, sep = "\t", header = T, as.is = T, 
                           stringsAsFactors = F, check.names = F, quote = "",
                           comment.char = "")

# common columns for mtb data
J.col_mtb <- c("SuperPathway", "SubPathway", "BiochemicalName", "ChemicalID", 
               "CAS", "RI", "Mass", "PathwaySortOrder", "Platform", "CompID", 
               "KEGG", "HMDB", "PUBCHEM")

### STEP 1: ----

# convert fold differences to log2 fold differences
# and repace NA's with 0
J.mtb_prox_1 <- J.mtb_prox_0
J.mtb_dist_1 <- J.mtb_dist_0
J.mtb_plas_1 <- J.mtb_plas_0
J.mtb_prox_1$log2fd <- log2(J.mtb_prox_1$FoldDiff)
J.mtb_dist_1$log2fd <- log2(J.mtb_dist_1$FoldDiff)
J.mtb_plas_1$log2fd <- log2(J.mtb_plas_1$FoldDiff)
J.mtb_prox_2 <- J.mtb_prox_1
J.mtb_dist_2 <- J.mtb_dist_1
J.mtb_plas_2 <- J.mtb_plas_1
J.mtb_prox_2[is.na(J.mtb_prox_2)] <- 0
J.mtb_dist_2[is.na(J.mtb_dist_2)] <- 0
J.mtb_plas_2[is.na(J.mtb_plas_2)] <- 0

# subset by sex of human donor
J.mtb_prox_F01_0 <- dplyr::filter(J.mtb_prox_2, group == "human female")
J.mtb_prox_M00_0 <- dplyr::filter(J.mtb_prox_2, group == "human male")
J.mtb_dist_F01_0 <- dplyr::filter(J.mtb_dist_2, group == "human female")
J.mtb_dist_M00_0 <- dplyr::filter(J.mtb_dist_2, group == "human male")
J.mtb_plas_F01_0 <- dplyr::filter(J.mtb_plas_2, group == "human female")
J.mtb_plas_M00_0 <- dplyr::filter(J.mtb_plas_2, group == "human male")

# convert col x to character and define breaks and labels for x axis
J.mtb_prox_F01_1 <- J.mtb_prox_F01_0
J.mtb_prox_M00_1 <- J.mtb_prox_M00_0
J.mtb_dist_F01_1 <- J.mtb_dist_F01_0
J.mtb_dist_M00_1 <- J.mtb_dist_M00_0
J.mtb_plas_F01_1 <- J.mtb_plas_F01_0
J.mtb_plas_M00_1 <- J.mtb_plas_M00_0
J.mtb_prox_F01_1$x <- as.character(J.mtb_prox_F01_1$x)
J.mtb_prox_M00_1$x <- as.character(J.mtb_prox_M00_1$x)
J.mtb_dist_F01_1$x <- as.character(J.mtb_dist_F01_1$x)
J.mtb_dist_M00_1$x <- as.character(J.mtb_dist_M00_1$x)
J.mtb_plas_F01_1$x <- as.character(J.mtb_plas_F01_1$x)
J.mtb_plas_M00_1$x <- as.character(J.mtb_plas_M00_1$x)

J.xbrk_mtb_prox_F01_1 <- J.mtb_prox_F01_1$x
J.xbrk_mtb_prox_M00_1 <- J.mtb_prox_M00_1$x
J.xbrk_mtb_dist_F01_1 <- J.mtb_dist_F01_1$x
J.xbrk_mtb_dist_M00_1 <- J.mtb_dist_M00_1$x
J.xbrk_mtb_plas_F01_1 <- J.mtb_plas_F01_1$x
J.xbrk_mtb_plas_M00_1 <- J.mtb_plas_M00_1$x

J.xlab_mtb_prox_F01_1 <- J.mtb_prox_F01_1$BiochemicalName
J.xlab_mtb_prox_M00_1 <- J.mtb_prox_M00_1$BiochemicalName
J.xlab_mtb_dist_F01_1 <- J.mtb_dist_F01_1$BiochemicalName
J.xlab_mtb_dist_M00_1 <- J.mtb_dist_M00_1$BiochemicalName
J.xlab_mtb_plas_F01_1 <- J.mtb_plas_F01_1$BiochemicalName
J.xlab_mtb_plas_M00_1 <- J.mtb_plas_M00_1$BiochemicalName

# the above for cros
J.mtb_cros_1 <- J.mtb_cros_0
J.mtb_cros_1$log2fd <- log2(J.mtb_cros_1$FoldDiff)
J.mtb_cros_2 <- J.mtb_cros_1
J.mtb_cros_2[is.na(J.mtb_cros_2)] <- 0
J.mtb_cros_F01_0 <- dplyr::filter(J.mtb_cros_2, group == "human female")
J.mtb_cros_M00_0 <- dplyr::filter(J.mtb_cros_2, group == "human male")
J.mtb_cros_F01_1 <- J.mtb_cros_F01_0
J.mtb_cros_M00_1 <- J.mtb_cros_M00_0
J.mtb_cros_F01_1$x <- as.character(J.mtb_cros_F01_1$x)
J.mtb_cros_M00_1$x <- as.character(J.mtb_cros_M00_1$x)
J.xbrk_mtb_cros_F01_1 <- J.mtb_cros_F01_1$x
J.xbrk_mtb_cros_M00_1 <- J.mtb_cros_M00_1$x
J.xlab_mtb_cros_F01_1 <- J.mtb_cros_F01_1$TypeLabel
J.xlab_mtb_cros_M00_1 <- J.mtb_cros_M00_1$TypeLabel

### parameters ---- 

# y-axis limits, breaks, labels
J.ylim <- c(-7, 7)
J.ybrk <- c(-6, -3, 0, 3, 6)
J.ylab <- c("-6", " -3", "0", " 3", "6")
J.yexp <- c(0, 0)

J.yttl <- "log2 fold difference"

# x-axis expand
J.xexp <- c(0, 0.6)

# bar sizing
J.wid <- 0.66

# scalers for significance asterisks
J.sclr_mtb <- 0.25
J.sze_sig <- 1.88

# x axis ticks
J.xtck_mtb_prox <- c(rep("solid", times = 11), "blank",
                     "solid", "blank",
                     "solid", "blank",
                     rep("solid", times = 07), "blank")
J.xtck_mtb_dist <- c(rep("solid", times = 02), 
                     "blank", "solid", "blank",
                     rep("solid", times = 05), "blank",
                     rep("solid", times = 02), "blank",
                     rep("solid", times = 04), "blank",
                     rep("solid", times = 05), "blank")
J.xtck_mtb_plas <- c(rep("solid", times = 04), "blank",
                     "solid", "blank",
                     "solid", "blank",
                     rep("solid", times = 09), "blank",
                     "solid", "blank",
                     rep("solid", times = 03), "blank")
J.xtck_mtb_cros <- c(rep("solid", times = 03), "blank", "blank",
                     rep("solid", times = 02), "blank", "blank",
                     rep("solid", times = 02), "blank", "blank",
                     rep("solid", times = 02), "blank", "blank")

J.xhex_mtb_prox <- c(rep(greydient[1], times = 23), greydient[8])
J.xhex_mtb_dist <- c(rep(greydient[1], times = 24), greydient[8])
J.xhex_mtb_plas <- c(rep(greydient[1], times = 24), greydient[8])
J.xhex_mtb_cros <- c(rep(greydient[1], times = 16), greydient[8])

# vector for medium colors by consortium group
J.hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2])
J.hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2])

# sizing
J.sze_atl <- 06 # size for axis titles (x and y)
J.sze_atx <- 07
J.sze_ltx <- 07 # size for legend text
J.sze_vlne <- 0.5

# horizontal line location and scaler
# NOTE: lines are vertical before changing orientation
J.vlne_mtb_prox <- c(16, 14, 12)
J.vlne_mtb_dist <- c(19, 14, 11, 05, 03)
J.vlne_mtb_plas <- c(21, 19, 09, 07, 05)
J.vlne_mtb_cros <- c(13, 09, 05)
J.vlne_sclr <- 0.2

# metabolite labels for cross conserved plot
# J.mtb_cros_lab_0 <- dplyr::distinct(J.mtb_cros_F01_1, BiochemicalName,
#                                     .keep_all = T)
# J.mtb_cros_lab_1 <- dplyr::select(J.mtb_cros_lab_0, BiochemicalName, x)
# J.mtb_cros_lab_2 <- dplyr::filter(J.mtb_cros_lab_1, 
#                                   !BiochemicalName == "blank space")
# J.mtb_cros_lab_3 <- J.mtb_cros_lab_2
# # J.mtb_cros_lab_3$x <- as.numeric(J.mtb_cros_lab_3$x + 0.5)

J.mtb_cros_lab_0 <- dplyr::distinct(J.mtb_cros_F01_1, BiochemicalLabel,
                                    .keep_all = T)
J.mtb_cros_lab_1 <- dplyr::select(J.mtb_cros_lab_0, BiochemicalLabel, x)
J.mtb_cros_lab_2 <- dplyr::filter(
  J.mtb_cros_lab_1, !BiochemicalLabel == "4-hydroxyphenylpyruvate")
J.mtb_cros_lab <- dplyr::filter(J.mtb_cros_lab_2, !BiochemicalLabel == "blank")

J.mtb_cros_lab_hyd <- dplyr::filter(
  J.mtb_cros_lab_1, BiochemicalLabel == "4-hydroxyphenylpyruvate")

# custom theme parameters
J.log2 <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]),
  axis.text = element_text(size = J.sze_atx),
  axis.title.x = element_text(size = J.sze_atl),
  axis.title.y = element_blank(),
  legend.background = element_rect(fill = NA),
  legend.key = element_rect(fill = NA, color = NA),
  legend.box.margin = margin(-3, 0, -6, 0, "mm"),
  legend.key.size = unit(2, "mm"),
  legend.title = element_text(size = J.sze_ltx),
  legend.text = element_text(
    margin = margin(0.5, 0, 0.5, 0, "mm"), size = J.sze_ltx, 
    color = greydient[1], hjust = 0),
  panel.ontop = T,
  panel.background = element_rect(fill = NA))

# plotting: mtb colon-prox ----

# disregard the Warning message:
# "Removed X rows containing missing values (position_stack)."
# these are the blank spaces inserted to separate metabolites by chemical class

J.gpr_mtb_prox_F01_1 <- ggbarplot(data = J.mtb_prox_F01_1, width = J.wid, 
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_prox + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  # geom_text(aes(x = x, label = Sig,
  #               y = ifelse(log2fd > 0, yes = log2fd + J.sclr_mtb,
  #                          no = log2fd - J.sclr_mtb)), 
  #           size = J.sze_sig, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_prox_F01_1, 
                   labels = J.xlab_mtb_prox_F01_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_F01, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_prox),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_prox)) +
  J.log2
# print(J.gpr_mtb_prox_F01_1)

J.gpr_mtb_prox_M00_1 <- ggbarplot(data = J.mtb_prox_M00_1, width = J.wid, 
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_prox + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  # geom_text(aes(x = x, label = Sig,
  #               y = ifelse(log2fd > 0, yes = log2fd + J.sclr_mtb,
  #                          no = log2fd - J.sclr_mtb)), 
  #           size = J.sze_sig, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_prox_M00_1, 
                   labels = J.xlab_mtb_prox_M00_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_M00, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_prox),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_prox)) +
  J.log2
# print(J.gpr_mtb_prox_M00_1)

# plotting: mtb colon-dist ----

# disregard the Warning message:
# "Removed X rows containing missing values (position_stack)."
# these are the blank spaces inserted to separate metabolites by chemical class

J.gpr_mtb_dist_F01_1 <- ggbarplot(data = J.mtb_dist_F01_1, width = J.wid, 
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_dist + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  # geom_text(aes(x = x, label = Sig,
  #               y = ifelse(log2fd > 0, yes = log2fd + J.sclr_mtb,
  #                          no = log2fd - J.sclr_mtb)), 
  #           size = J.sze_sig, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_dist_F01_1, 
                   labels = J.xlab_mtb_dist_F01_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_F01, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_dist),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_dist)) +
  J.log2
# print(J.gpr_mtb_dist_F01_1)

J.gpr_mtb_dist_M00_1 <- ggbarplot(data = J.mtb_dist_M00_1, width = J.wid, 
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_dist + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  # geom_text(aes(x = x, label = Sig,
  #               y = ifelse(log2fd > 0, yes = log2fd + J.sclr_mtb,
  #                          no = log2fd - J.sclr_mtb)), 
  #           size = J.sze_sig, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_dist_M00_1, 
                   labels = J.xlab_mtb_dist_M00_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_M00, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_dist),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_dist)) +
  J.log2
# print(J.gpr_mtb_dist_M00_1)

# plotting: mtb plasma ----

# disregard the Warning message:
# "Removed X rows containing missing values (position_stack)."
# these are the blank spaces inserted to separate metabolites by chemical class

J.gpr_mtb_plas_F01_1 <- ggbarplot(data = J.mtb_plas_F01_1, width = J.wid,
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_plas + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  # geom_text(aes(x = x, label = Sig,
  #               y = ifelse(log2fd > 0, yes = log2fd + J.sclr_mtb,
  #                          no = log2fd - J.sclr_mtb)),
  #           size = J.sze_sig, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_plas_F01_1,
                   labels = J.xlab_mtb_plas_F01_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_F01, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_plas),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_plas)) +
  J.log2
# print(J.gpr_mtb_plas_F01_1)

J.gpr_mtb_plas_M00_1 <- ggbarplot(data = J.mtb_plas_M00_1, width = J.wid,
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_plas + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  # geom_text(aes(x = x, label = Sig,
  #               y = ifelse(log2fd > 0, yes = log2fd + J.sclr_mtb,
  #                          no = log2fd - J.sclr_mtb)),
  #           size = J.sze_sig, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_plas_M00_1,
                   labels = J.xlab_mtb_plas_M00_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_M00, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_plas),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_plas)) +
  J.log2
# print(J.gpr_mtb_plas_M00_1)

# plotting: mtb cross ----

J.sze_mtb_lab <- 1.95
J.sze_mtb_lab_hyd <- 1.6
J.scl_mtb_lab <- 1

J.gpr_mtb_cros_F01_1 <- ggbarplot(data = J.mtb_cros_F01_1, width = J.wid,
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_cros + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  geom_text(data = J.mtb_cros_lab, aes(x = x, label = BiochemicalLabel),
            y = 0, size = J.sze_mtb_lab, family = fnt, color = greydient[1]) +
  geom_text(data = J.mtb_cros_lab_hyd, aes(x = x, label = BiochemicalLabel),
            y = 0, size = J.sze_mtb_lab_hyd, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_cros_F01_1,
                   labels = J.xlab_mtb_cros_F01_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_F01, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_cros),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_cros)) +
  J.log2
# print(J.gpr_mtb_cros_F01_1)

J.gpr_mtb_cros_M00_1 <- ggbarplot(data = J.mtb_cros_M00_1, width = J.wid,
                                  x = "x", y = "log2fd", fill = "Increase",
                                  color = greydient[8], font.family = fnt,
                                  orientation = "horiz") +
  labs(y = J.yttl) +
  geom_vline(xintercept = J.vlne_mtb_cros + J.vlne_sclr, size = J.sze_vlne,
             color = greydient[5], linetype = "dotted") +
  geom_text(data = J.mtb_cros_lab, aes(x = x, label = BiochemicalLabel),
            y = 0, size = J.sze_mtb_lab, family = fnt, color = greydient[1]) +
  geom_text(data = J.mtb_cros_lab_hyd, aes(x = x, label = BiochemicalLabel),
            y = 0, size = J.sze_mtb_lab_hyd, family = fnt, color = greydient[1]) +
  scale_x_discrete(breaks = J.xbrk_mtb_cros_M00_1,
                   labels = J.xlab_mtb_cros_M00_1) +
  scale_y_continuous(limits = J.ylim, breaks = J.ybrk, labels = J.ylab,
                     expand = J.yexp) +
  scale_fill_manual(values = J.hex_M00, name = NULL) +
  border(color = greydient[1]) +
  theme(axis.ticks.y = element_line(linetype = J.xtck_mtb_cros),
        axis.text.y = element_text(size = J.sze_atx, color = J.xhex_mtb_cros)) +
  J.log2
# print(J.gpr_mtb_cros_M00_1)

# plotting: arrange ----

J.gga_mtb_prox <- ggarrange(J.gpr_mtb_prox_F01_1, J.gpr_mtb_prox_M00_1,
                            labels = c("A", ""), font.label = pan_fnt, 
                            ncol = 2, nrow = 1, align = "hv")
J.gga_mtb_dist <- ggarrange(J.gpr_mtb_dist_F01_1, J.gpr_mtb_dist_M00_1,
                            labels = c("B", ""), font.label = pan_fnt, 
                            ncol = 2, nrow = 1, align = "hv")
J.gga_mtb_plas <- ggarrange(J.gpr_mtb_plas_F01_1, J.gpr_mtb_plas_M00_1,
                            labels = c("C", ""), font.label = pan_fnt,
                            ncol = 2, nrow = 1, align = "hv")

### ************************************
### J - WRITE OUTPUTS ----
### ************************************

# J.gga_mtb_prox # 24
# J.gga_mtb_dist # 25
# J.gga_mtb_plas # 25
# J.gga_mtb_cros # 17

# 24 * 3.4375 # 82.5
# 25 * 3.4375 # 85.9375
# 17 * 3.4375 # 58.4375

# output file paths
J.ofp_plot_mtb_prox_F01 <- "TransFaunation/vault/plot_mtb_prox_F01.pdf"
J.ofp_plot_mtb_dist_F01 <- "TransFaunation/vault/plot_mtb_dist_F01.pdf"
J.ofp_plot_mtb_plas_F01 <- "TransFaunation/vault/plot_mtb_plas_F01.pdf"
J.ofp_plot_mtb_cros_F01 <- "TransFaunation/vault/plot_mtb_cros_F01.pdf"

J.ofp_plot_mtb_prox_M00 <- "TransFaunation/vault/plot_mtb_prox_M00.pdf"
J.ofp_plot_mtb_dist_M00 <- "TransFaunation/vault/plot_mtb_dist_M00.pdf"
J.ofp_plot_mtb_plas_M00 <- "TransFaunation/vault/plot_mtb_plas_M00.pdf"
J.ofp_plot_mtb_cros_M00 <- "TransFaunation/vault/plot_mtb_cros_M00.pdf"


ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 82.5,
       filename = J.ofp_plot_mtb_prox_F01, plot = J.gpr_mtb_prox_F01_1)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 85.9375,
       filename = J.ofp_plot_mtb_dist_F01, plot = J.gpr_mtb_dist_F01_1)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 82.5,
       filename = J.ofp_plot_mtb_plas_F01, plot = J.gpr_mtb_plas_F01_1)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 58.4375,
       filename = J.ofp_plot_mtb_cros_F01, plot = J.gpr_mtb_cros_F01_1)

ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 82.5,
       filename = J.ofp_plot_mtb_prox_M00, plot = J.gpr_mtb_prox_M00_1)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 85.9375,
       filename = J.ofp_plot_mtb_dist_M00, plot = J.gpr_mtb_dist_M00_1)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 82.5,
       filename = J.ofp_plot_mtb_plas_M00, plot = J.gpr_mtb_plas_M00_1)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_sgl, height = 58.4375,
       filename = J.ofp_plot_mtb_cros_M00, plot = J.gpr_mtb_cros_M00_1)

# write.table(sep = "\t", row.names = F, x = J.mtb_prox_F01_1, 
#             file = "J_mtb_prox_F01.txt")
# write.table(sep = "\t", row.names = F, x = J.mtb_dist_F01_1, 
#             file = "J_mtb_dist_F01.txt")
# write.table(sep = "\t", row.names = F, x = J.mtb_plas_F01_1, 
#             file = "J_mtb_plas_F01.txt")
# 
# write.table(sep = "\t", row.names = F, x = J.mtb_prox_M00_1, 
#             file = "J_mtb_prox_M00.txt")
# write.table(sep = "\t", row.names = F, x = J.mtb_dist_M00_1, 
#             file = "J_mtb_dist_M00.txt")
# write.table(sep = "\t", row.names = F, x = J.mtb_plas_M00_1, 
#             file = "J_mtb_plas_M00.txt")

# save workspace
J.obj <- ls(pattern = "J.")
J.lst <- c(J.obj[grep(pattern = "J.", x = J.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS)
save(list = J.lst, file = J.ofv_wksp)

### ************************************
### BEGIN Section L ----
### ************************************

# note for kdp: taken from TransFaunation/section_B.R

# NEED to add murine sex-based stats and plots <- is this still true?
# check verbiage and ensure that tumor/lesion is always referenced as lesion

### preface ----

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
load(file = ofv_COSMOS_wksp)

L.ofv_wksp <- paste(sep = "", path_vault, "/WS_", "SectionL.RData")

### ************************************
### L - FUNCTION ----
### ************************************

# create a vector naming all Section A functions (used when saving workspace)
L.function <- "format_wi_kw"

# ** note for KDP: version 0.1 ** #
# function takexs an input of one of the following classes:
# "htest" from kruskal.test()
# "pairwise.htest" from pairwise.wilcox.test() 
# "PMCMR" from PMCMR::posthoc.kruskal.dunn.test()
# and returns a formatted results data.frame

format_wi_kw <- function(data, digits, cutpoints, symbols, names_P = c("","")) {
  
  # internal checks for correct input types
  if (!inherits(data, "htest") && !inherits(data, "pairwise.htest") 
      && !inherits(data, "PMCMR")) {
    stop("data must be of class 'htest' or 'pairwise.htest' or 'PMCMR'")
  }
  if (!inherits(digits, "numeric")) {
    stop("input for digits must be class numeric")
  }
  if (!inherits(cutpoints, "numeric")) {
    stop("input for cutpoints must be class numeric")
  }
  if (!inherits(symbols, "character")) {
    stop("input for symbol must be class character")
  }
  if (!inherits(names_P, "character")) {
    stop("input for names_P must be class character")
  }
  if (!length(names_P) == 2) {
    stop("input for names_P must be vector of length 2")
  }
  
  # input data are from kruskal.test()
  # the process occurs as follows:
  # (_0) create a data.frame of relevant information and p.values
  # (_1) create col 'sig_P' with symbols for significance of P values
  # (_2) coerce col with symbols for significance of P values to character
  # (_3) round P values to the specified number of digits
  # (_4) reorder and retain relevant columns
  # (_5) rename columns 'P' and 'sig_P' using the input values for names_P
  if (inherits(data, "htest")) {
    kw_0 <- data.frame("Test" = "Kruskal-Wallis", "Comparison" = "global", 
                       "P" = data$p.value, stringsAsFactors = F)
    kw_1 <- dplyr::mutate(kw_0, sig_P = symnum(P, cutpoints = cutpoints, 
                                               symbols = symbols, corr = F))
    kw_2 <- dplyr::mutate(kw_1, sig_P = as.character(sig_P))
    kw_3 <- dplyr::mutate(kw_2, P = round(P, digits = digits))
    kw_4 <- dplyr::select(kw_3, Test, Comparison, P, sig_P)
    kw_5 <- dplyr::rename_at(kw_4, dplyr::vars(P, sig_P), ~names_P)
    return(kw_5)
  }
  # NOTE: the process is identical for both 'pairwise.htest' and 'PMCMR' ... 
  # and occurs as follows:
  # (_0) create a data.frame of the p.values
  # (_1) create a data.frame of additional relevant information and p.values
  # (_2) reshape the data.frame into a workable format
  # (_3) remove any rows with NA vals (these represent same group comparisons)
  # (_4) rename column 'variable' to 'Pair2'
  # (_5) create col 'Comparison' by combining rows in cols 'Pair1' & 'Pair2'
  # (_6) create col 'sig_P' with symbols for significance of P values
  # (_7) coerce col with symbols for significance of P values to character
  # (_8) round P values to the specified number of digits
  # (_9) reorder and retain relevant columns
  # (_10) rename columns 'P' and 'sig_P' using the input values for names_P
  
  # input data are from pairwise.wilcox.test() 
  if (inherits(data, "pairwise.htest")) {
    wi_0 <- data.frame(data$p.value, stringsAsFactors = F)
    wi_1 <- data.frame("Test" = "Wilcoxon", "Pair1" = row.names(wi_0), wi_0, 
                       stringsAsFactors = F)
    wi_2 <- reshape2::melt(wi_1, id.vars = c("Test", "Pair1"), value.name = "P")
    wi_3 <- dplyr::filter(wi_2, !is.na(P))
    wi_4 <- dplyr::rename(wi_3, Pair2 = variable)
    wi_5 <- dplyr::mutate(wi_4, 
                          Comparison = paste(Pair1, "vs", Pair2, sep = " "))
    wi_6 <- dplyr::mutate(wi_5, sig_P = symnum(P, cutpoints = cutpoints, 
                                               symbols = symbols, corr = F))
    wi_7 <- dplyr::mutate(wi_6, sig_P = as.character(sig_P))
    wi_8 <- dplyr::mutate(wi_7, P = round(P, digits = digits))
    wi_9 <- dplyr::select(wi_8, Test, Comparison, P, sig_P)
    wi_10 <- dplyr::rename_at(wi_9, dplyr::vars(P, sig_P), ~names_P)
    return(wi_10)
  }
  # input data are from PMCMR::posthoc.kruskal.dunn.test()
  if (inherits(data, "PMCMR")) {
    du_0 <- data.frame(data$p.value, stringsAsFactors = F)
    du_1 <- data.frame("Test" = "Dunn's Test", "Pair1" = row.names(du_0), du_0, 
                       stringsAsFactors = F)
    du_2 <- reshape2::melt(du_1, id.vars = c("Test", "Pair1"), value.name = "P")
    du_3 <- dplyr::filter(du_2, !is.na(P))
    du_4 <- dplyr::rename(du_3, Pair2 = variable)
    du_5 <- dplyr::mutate(du_4, 
                          Comparison = paste(Pair1, "vs", Pair2, sep = " "))
    du_6 <- dplyr::mutate(du_5, sig_P = symnum(P, cutpoints = cutpoints, 
                                               symbols = symbols, corr = F))
    du_7 <- dplyr::mutate(du_6, sig_P = as.character(sig_P))
    du_8 <- dplyr::mutate(du_7, P = round(P, digits = digits))
    du_9 <- dplyr::select(du_8, Test, Comparison, P, sig_P)
    du_10 <- dplyr::rename_at(du_9, dplyr::vars(P, sig_P), ~names_P)
    return(du_10)
  }
}
#
# example usage:
# new.df <- format_wi_kw(data = results, digits = 3, 
#                        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
#                        symbols = c("****", "***", "**", "*", " ")
#                        names_P = c("BH_P", "sig_BH_P"))

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

### ************************************
### L - STEP  1 - format data for statistics and plotting ----
### ************************************

# NOTE: lesion data are located in columns of the sample data

# format sample data for comparisons across consortium groups
# this process occurs as follows:
# (1;2) retain relevant cols; filter to isolate lesion data
# (3;4) copy data.frame & convert columns with lesion numbers to numeric

L.les_0 <- dplyr::select(smp_dat, LesionSampleID, ConsortiumAbrv,
                         ConsortiumAbrv, DonorGroupConsortium, 
                         ConsortiumAnimalSex, AnimalSex, HumanDonorSex, 
                         LesionFilter, LesionTotal, 
                         LesionTotPrx, LesionTotMid, LesionTotDtl, 
                         LesionSizeOnePrx, LesionSizeOneTwoPrx, 
                         LesionSizeTwoThreePrx, LesionSizeOneMid, 
                         LesionSizeOneTwoMid, LesionSizeTwoThreeMid, 
                         LesionSizeOneDtl, LesionSizeOneTwoDtl, 
                         LesionSizeTwoThreeDtl)
L.les_1 <- dplyr::filter(L.les_0, LesionFilter == "yay")

L.les_2 <- L.les_1
L.les_2$LesionTotal <- as.numeric(L.les_2$LesionTotal)
L.les_2$LesionTotPrx <- as.numeric(L.les_2$LesionTotPrx)
L.les_2$LesionTotMid <- as.numeric(L.les_2$LesionTotMid)
L.les_2$LesionTotDtl <- as.numeric(L.les_2$LesionTotDtl)
L.les_2$LesionSizeOnePrx <- as.numeric(L.les_2$LesionSizeOnePrx)
L.les_2$LesionSizeOneTwoPrx <- as.numeric(L.les_2$LesionSizeOneTwoPrx)
L.les_2$LesionSizeTwoThreePrx <- as.numeric(L.les_2$LesionSizeTwoThreePrx)
L.les_2$LesionSizeOneMid <- as.numeric(L.les_2$LesionSizeOneMid)
L.les_2$LesionSizeOneTwoMid <- as.numeric(L.les_2$LesionSizeOneTwoMid)
L.les_2$LesionSizeTwoThreeMid <- as.numeric(L.les_2$LesionSizeTwoThreeMid)
L.les_2$LesionSizeOneDtl <- as.numeric(L.les_2$LesionSizeOneDtl)
L.les_2$LesionSizeOneTwoDtl <- as.numeric(L.les_2$LesionSizeOneTwoDtl)
L.les_2$LesionSizeTwoThreeDtl <- as.numeric(L.les_2$LesionSizeTwoThreeDtl)

# perform normality test
L.nor_tot <- shapiro.test(L.les_2$LesionTotal)
L.nor_prx <- shapiro.test(L.les_2$LesionTotPrx)
L.nor_mid <- shapiro.test(L.les_2$LesionTotMid)
L.nor_dtl <- shapiro.test(L.les_2$LesionTotDtl)

if (!isTRUE(L.nor_tot$p.value > 0.05)) {
  message("L.nor_tot: normality test failed")
}
if (!isTRUE(L.nor_prx$p.value > 0.05)) {
  message("L.nor_prx: normality test failed")
}
if (!isTRUE(L.nor_mid$p.value > 0.05)) {
  message("L.nor_mid: normality test failed")
}
if (!isTRUE(L.nor_dtl$p.value > 0.05)) {
  message("L.nor_dtl: normality test failed")
}

# (1;2) subset to obtain data for desired comparisons
L.tot <- dplyr::select(L.les_2, LesionSampleID, LesionTotal,
                       HumanDonorSex, ConsortiumAbrv, 
                       ConsortiumAnimalSex, AnimalSex)
L.prx <- dplyr::select(L.les_2, LesionSampleID, LesionTotPrx,
                       HumanDonorSex, ConsortiumAbrv, 
                       ConsortiumAnimalSex, AnimalSex)
L.mid <- dplyr::select(L.les_2, LesionSampleID, LesionTotMid,
                       HumanDonorSex, ConsortiumAbrv,
                       ConsortiumAnimalSex, AnimalSex)
L.dtl <- dplyr::select(L.les_2, LesionSampleID, LesionTotDtl,
                       HumanDonorSex, ConsortiumAbrv, 
                       ConsortiumAnimalSex, AnimalSex)
L.tot_F01_0 <- dplyr::filter(L.tot, HumanDonorSex == "Female" |
                               HumanDonorSex == "not.app")
L.prx_F01_0 <- dplyr::filter(L.prx, HumanDonorSex == "Female" |
                               HumanDonorSex == "not.app")
L.mid_F01_0 <- dplyr::filter(L.mid, HumanDonorSex == "Female" |
                               HumanDonorSex == "not.app")
L.dtl_F01_0 <- dplyr::filter(L.dtl, HumanDonorSex == "Female" |
                               HumanDonorSex == "not.app")
L.tot_M00_0 <- dplyr::filter(L.tot, HumanDonorSex == "Male" |
                               HumanDonorSex == "not.app")
L.prx_M00_0 <- dplyr::filter(L.prx, HumanDonorSex == "Male" |
                               HumanDonorSex == "not.app")
L.mid_M00_0 <- dplyr::filter(L.mid, HumanDonorSex == "Male" |
                               HumanDonorSex == "not.app")
L.dtl_M00_0 <- dplyr::filter(L.dtl, HumanDonorSex == "Male" |
                               HumanDonorSex == "not.app")

# (3) convert column used as grouping into factors
# consortium group comparisons
L.tot_F01 <- dplyr::mutate(
  L.tot_F01_0, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))
L.tot_M00 <- dplyr::mutate(
  L.tot_M00_0, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))

# murine-sex based comparisons
L.tot_sex_F01 <- dplyr::mutate(
  L.tot_F01_0, 
  ConsortiumAnimalSex = factor(ConsortiumAnimalSex, 
                               levels = unique(ConsortiumAnimalSex)))
L.tot_sex_M00 <- dplyr::mutate(
  L.tot_M00_0, 
  ConsortiumAnimalSex = factor(ConsortiumAnimalSex, 
                               levels = unique(ConsortiumAnimalSex)))

# for plots, add relevant x-axis information to location specific data.frames
L.prx_F01_1 <- L.prx_F01_0
L.mid_F01_1 <- L.mid_F01_0
L.dtl_F01_1 <- L.dtl_F01_0
L.prx_M00_1 <- L.prx_M00_0
L.mid_M00_1 <- L.mid_M00_0
L.dtl_M00_1 <- L.dtl_M00_0

L.prx_F01_1$ConsortiumAbrvAlt[L.prx_F01_1$ConsortiumAbrv == "N"] <- "N.1"
L.prx_F01_1$ConsortiumAbrvAlt[L.prx_F01_1$ConsortiumAbrv == "CMC"] <- "CMC.1"
L.prx_F01_1$ConsortiumAbrvAlt[L.prx_F01_1$ConsortiumAbrv == "RMC"] <- "RMC.1"
L.mid_F01_1$ConsortiumAbrvAlt[L.mid_F01_1$ConsortiumAbrv == "N"] <- "N.2"
L.mid_F01_1$ConsortiumAbrvAlt[L.mid_F01_1$ConsortiumAbrv == "CMC"] <- "CMC.2"
L.mid_F01_1$ConsortiumAbrvAlt[L.mid_F01_1$ConsortiumAbrv == "RMC"] <- "RMC.2"
L.dtl_F01_1$ConsortiumAbrvAlt[L.dtl_F01_1$ConsortiumAbrv == "N"] <- "N.3"
L.dtl_F01_1$ConsortiumAbrvAlt[L.dtl_F01_1$ConsortiumAbrv == "CMC"] <- "CMC.3"
L.dtl_F01_1$ConsortiumAbrvAlt[L.dtl_F01_1$ConsortiumAbrv == "RMC"] <- "RMC.3"
L.prx_M00_1$ConsortiumAbrvAlt[L.prx_M00_1$ConsortiumAbrv == "N"] <- "N.1"
L.prx_M00_1$ConsortiumAbrvAlt[L.prx_M00_1$ConsortiumAbrv == "CMC"] <- "CMC.1"
L.prx_M00_1$ConsortiumAbrvAlt[L.prx_M00_1$ConsortiumAbrv == "RMC"] <- "RMC.1"
L.mid_M00_1$ConsortiumAbrvAlt[L.mid_M00_1$ConsortiumAbrv == "N"] <- "N.2"
L.mid_M00_1$ConsortiumAbrvAlt[L.mid_M00_1$ConsortiumAbrv == "CMC"] <- "CMC.2"
L.mid_M00_1$ConsortiumAbrvAlt[L.mid_M00_1$ConsortiumAbrv == "RMC"] <- "RMC.2"
L.dtl_M00_1$ConsortiumAbrvAlt[L.dtl_M00_1$ConsortiumAbrv == "N"] <- "N.3"
L.dtl_M00_1$ConsortiumAbrvAlt[L.dtl_M00_1$ConsortiumAbrv == "CMC"] <- "CMC.3"
L.dtl_M00_1$ConsortiumAbrvAlt[L.dtl_M00_1$ConsortiumAbrv == "RMC"] <- "RMC.3"

L.prx_F01_2 <- dplyr::rename(L.prx_F01_1, LesionTot = LesionTotPrx)
L.mid_F01_2 <- dplyr::rename(L.mid_F01_1, LesionTot = LesionTotMid)
L.dtl_F01_2 <- dplyr::rename(L.dtl_F01_1, LesionTot = LesionTotDtl)
L.prx_M00_2 <- dplyr::rename(L.prx_M00_1, LesionTot = LesionTotPrx)
L.mid_M00_2 <- dplyr::rename(L.mid_M00_1, LesionTot = LesionTotMid)
L.dtl_M00_2 <- dplyr::rename(L.dtl_M00_1, LesionTot = LesionTotDtl)

L.prx_F01 <- dplyr::mutate(
  L.prx_F01_2, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))
L.mid_F01 <- dplyr::mutate(
  L.mid_F01_2, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))
L.dtl_F01 <- dplyr::mutate(
  L.dtl_F01_2, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))
L.prx_M00 <- dplyr::mutate(
  L.prx_M00_2, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))
L.mid_M00 <- dplyr::mutate(
  L.mid_M00_2, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))
L.dtl_M00 <- dplyr::mutate(
  L.dtl_M00_2, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))

L.loc_F01_0 <- rbind(L.prx_F01, L.mid_F01, L.dtl_F01)
L.loc_M00_0 <- rbind(L.prx_M00, L.mid_M00, L.dtl_M00)

L.loc_F01 <- dplyr::mutate(
  L.loc_F01_0, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))
L.loc_M00 <- dplyr::mutate(
  L.loc_M00_0, 
  ConsortiumAbrv = factor(ConsortiumAbrv, levels = unique(ConsortiumAbrv)))

### ************************************
### L - STEP  2 - perform global tests and format results ----
### ************************************

# perform 'global' test (i.e. All vs. All); this process occurs as follows:
# (1) perform test (NOTE: kruskal-wallis used as there are three groups)
# (2) pass test outputs to the format_wi_kw() function to format results
# (3;4) add new cols 'BH_P' and 'sig_BH_P' for merging w/ pwse results later on

L.tot_F01_glbl_0 <- kruskal.test(LesionTotal ~ ConsortiumAbrv, data = L.tot_F01)
L.prx_F01_glbl_0 <- kruskal.test(LesionTot ~ ConsortiumAbrv, data = L.prx_F01)
L.mid_F01_glbl_0 <- kruskal.test(LesionTot ~ ConsortiumAbrv, data = L.mid_F01)
L.dtl_F01_glbl_0 <- kruskal.test(LesionTot ~ ConsortiumAbrv, data = L.dtl_F01)
L.tot_M00_glbl_0 <- kruskal.test(LesionTotal ~ ConsortiumAbrv, data = L.tot_M00)
L.prx_M00_glbl_0 <- kruskal.test(LesionTot ~ ConsortiumAbrv, data = L.prx_M00)
L.mid_M00_glbl_0 <- kruskal.test(LesionTot ~ ConsortiumAbrv, data = L.mid_M00)
L.dtl_M00_glbl_0 <- kruskal.test(LesionTot ~ ConsortiumAbrv, data = L.dtl_M00)

L.tot_F01_glbl_1 <- format_wi_kw(data = L.tot_F01_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))
L.prx_F01_glbl_1 <- format_wi_kw(data = L.prx_F01_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))
L.mid_F01_glbl_1 <- format_wi_kw(data = L.mid_F01_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))
L.dtl_F01_glbl_1 <- format_wi_kw(data = L.dtl_F01_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))
L.tot_M00_glbl_1 <- format_wi_kw(data = L.tot_M00_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))
L.prx_M00_glbl_1 <- format_wi_kw(data = L.prx_M00_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))
L.mid_M00_glbl_1 <- format_wi_kw(data = L.mid_M00_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))
L.dtl_M00_glbl_1 <- format_wi_kw(data = L.dtl_M00_glbl_0, 
                                 digits = 3, cutpoints = cutpts, 
                                 symbols = symbls, names_P = c("P", "sig_P"))

L.tot_F01_glbl <- L.tot_F01_glbl_1
L.prx_F01_glbl <- L.prx_F01_glbl_1
L.mid_F01_glbl <- L.mid_F01_glbl_1
L.dtl_F01_glbl <- L.dtl_F01_glbl_1
L.tot_M00_glbl <- L.tot_M00_glbl_1
L.prx_M00_glbl <- L.prx_M00_glbl_1
L.mid_M00_glbl <- L.mid_M00_glbl_1
L.dtl_M00_glbl <- L.dtl_M00_glbl_1

L.tot_F01_glbl$BH_P <- "-"
L.prx_F01_glbl$BH_P <- "-"
L.mid_F01_glbl$BH_P <- "-"
L.dtl_F01_glbl$BH_P <- "-"
L.tot_M00_glbl$BH_P <- "-"
L.prx_M00_glbl$BH_P <- "-"
L.mid_M00_glbl$BH_P <- "-"
L.dtl_M00_glbl$BH_P <- "-"
L.tot_F01_glbl$sig_BH_P <- "-"
L.prx_F01_glbl$sig_BH_P <- "-"
L.mid_F01_glbl$sig_BH_P <- "-"
L.dtl_F01_glbl$sig_BH_P <- "-"
L.tot_M00_glbl$sig_BH_P <- "-"
L.prx_M00_glbl$sig_BH_P <- "-"
L.mid_M00_glbl$sig_BH_P <- "-"
L.dtl_M00_glbl$sig_BH_P <- "-"

### ************************************
### L - STEP  3 - perform pairwise tests and format results ----
### ************************************

# perform pairwise comparisons; this process occurs as follows:
# (1) perform pairwise tests without multiple comparisons adjustment
# (2) perform pairwise tests with multiple comparisons adjustment (fdr/BH)
# (3) pass test outputs to the format_wi_kw() function to format results
# (4) merge appropriate _raw and _fdr data.frames by col 'Comparison'
# (5) combine appropriate _pwwi and _pwkw data.frames

L.tot_F01_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.tot_F01$LesionTotal, 
  g = L.tot_F01$ConsortiumAbrv, 
  p.adjust.method = "none")
L.prx_F01_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.prx_F01$LesionTot,
  g = L.prx_F01$ConsortiumAbrv,
  p.adjust.method = "none")
L.mid_F01_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.mid_F01$LesionTot,
  g = L.mid_F01$ConsortiumAbrv,
  p.adjust.method = "none")
L.dtl_F01_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.dtl_F01$LesionTot,
  g = L.dtl_F01$ConsortiumAbrv,
  p.adjust.method = "none")
L.tot_M00_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.tot_M00$LesionTotal,
  g = L.tot_M00$ConsortiumAbrv,
  p.adjust.method = "none")
L.prx_M00_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.prx_M00$LesionTot,
  g = L.prx_M00$ConsortiumAbrv,
  p.adjust.method = "none")
L.mid_M00_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.mid_M00$LesionTot,
  g = L.mid_M00$ConsortiumAbrv,
  p.adjust.method = "none")
L.dtl_M00_pwkw_raw_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.dtl_M00$LesionTot,
  g = L.dtl_M00$ConsortiumAbrv,
  p.adjust.method = "none")

L.tot_F01_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.tot_F01$LesionTotal, 
  g = L.tot_F01$ConsortiumAbrv, 
  p.adjust.method = "fdr")
L.prx_F01_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.prx_F01$LesionTot,
  g = L.prx_F01$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.mid_F01_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.mid_F01$LesionTot,
  g = L.mid_F01$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.dtl_F01_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.dtl_F01$LesionTot,
  g = L.dtl_F01$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.tot_M00_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.tot_M00$LesionTotal,
  g = L.tot_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.prx_M00_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.prx_M00$LesionTot,
  g = L.prx_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.mid_M00_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.mid_M00$LesionTot,
  g = L.mid_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.dtl_M00_pwkw_fdr_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = L.dtl_M00$LesionTot,
  g = L.dtl_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")

L.tot_F01_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.tot_F01$LesionTotal, 
  g = L.tot_F01$ConsortiumAbrv, 
  p.adjust.method = "none")
L.prx_F01_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.prx_F01$LesionTot,
  g = L.prx_F01$ConsortiumAbrv,
  p.adjust.method = "none")
L.mid_F01_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.mid_F01$LesionTot,
  g = L.mid_F01$ConsortiumAbrv,
  p.adjust.method = "none")
L.dtl_F01_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.dtl_F01$LesionTot,
  g = L.dtl_F01$ConsortiumAbrv,
  p.adjust.method = "none")
L.tot_M00_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.tot_M00$LesionTotal,
  g = L.tot_M00$ConsortiumAbrv,
  p.adjust.method = "none")
L.prx_M00_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.prx_M00$LesionTot,
  g = L.prx_M00$ConsortiumAbrv,
  p.adjust.method = "none")
L.mid_M00_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.mid_M00$LesionTot,
  g = L.mid_M00$ConsortiumAbrv,
  p.adjust.method = "none")
L.dtl_M00_pwwi_raw_0 <- pairwise.wilcox.test(
  x = L.dtl_M00$LesionTot,
  g = L.dtl_M00$ConsortiumAbrv,
  p.adjust.method = "none")

L.tot_F01_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.tot_F01$LesionTotal, 
  g = L.tot_F01$ConsortiumAbrv, 
  p.adjust.method = "fdr")
L.prx_F01_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.prx_F01$LesionTot,
  g = L.prx_F01$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.mid_F01_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.mid_F01$LesionTot,
  g = L.mid_F01$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.dtl_F01_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.dtl_F01$LesionTot,
  g = L.dtl_F01$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.tot_M00_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.tot_M00$LesionTotal,
  g = L.tot_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.prx_M00_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.prx_M00$LesionTot,
  g = L.prx_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.mid_M00_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.mid_M00$LesionTot,
  g = L.mid_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")
L.dtl_M00_pwwi_fdr_0 <- pairwise.wilcox.test(
  x = L.dtl_M00$LesionTot,
  g = L.dtl_M00$ConsortiumAbrv,
  p.adjust.method = "fdr")

L.tot_F01_pwkw_raw_1 <- format_wi_kw(data = L.tot_F01_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.prx_F01_pwkw_raw_1 <- format_wi_kw(data = L.prx_F01_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.mid_F01_pwkw_raw_1 <- format_wi_kw(data = L.mid_F01_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.dtl_F01_pwkw_raw_1 <- format_wi_kw(data = L.dtl_F01_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.tot_M00_pwkw_raw_1 <- format_wi_kw(data = L.tot_M00_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.prx_M00_pwkw_raw_1 <- format_wi_kw(data = L.prx_M00_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.mid_M00_pwkw_raw_1 <- format_wi_kw(data = L.mid_M00_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.dtl_M00_pwkw_raw_1 <- format_wi_kw(data = L.dtl_M00_pwkw_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))

L.tot_F01_pwkw_fdr_1 <- format_wi_kw(data = L.tot_F01_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.prx_F01_pwkw_fdr_1 <- format_wi_kw(data = L.prx_F01_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.mid_F01_pwkw_fdr_1 <- format_wi_kw(data = L.mid_F01_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.dtl_F01_pwkw_fdr_1 <- format_wi_kw(data = L.dtl_F01_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.tot_M00_pwkw_fdr_1 <- format_wi_kw(data = L.tot_M00_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.prx_M00_pwkw_fdr_1 <- format_wi_kw(data = L.prx_M00_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.mid_M00_pwkw_fdr_1 <- format_wi_kw(data = L.mid_M00_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.dtl_M00_pwkw_fdr_1 <- format_wi_kw(data = L.dtl_M00_pwkw_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))

L.tot_F01_pwwi_raw_1 <- format_wi_kw(data = L.tot_F01_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.prx_F01_pwwi_raw_1 <- format_wi_kw(data = L.prx_F01_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.mid_F01_pwwi_raw_1 <- format_wi_kw(data = L.mid_F01_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.dtl_F01_pwwi_raw_1 <- format_wi_kw(data = L.dtl_F01_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.tot_M00_pwwi_raw_1 <- format_wi_kw(data = L.tot_M00_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.prx_M00_pwwi_raw_1 <- format_wi_kw(data = L.prx_M00_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.mid_M00_pwwi_raw_1 <- format_wi_kw(data = L.mid_M00_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))
L.dtl_M00_pwwi_raw_1 <- format_wi_kw(data = L.dtl_M00_pwwi_raw_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("P", "sig_P"))

L.tot_F01_pwwi_fdr_1 <- format_wi_kw(data = L.tot_F01_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.prx_F01_pwwi_fdr_1 <- format_wi_kw(data = L.prx_F01_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.mid_F01_pwwi_fdr_1 <- format_wi_kw(data = L.mid_F01_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.dtl_F01_pwwi_fdr_1 <- format_wi_kw(data = L.dtl_F01_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.tot_M00_pwwi_fdr_1 <- format_wi_kw(data = L.tot_M00_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.prx_M00_pwwi_fdr_1 <- format_wi_kw(data = L.prx_M00_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.mid_M00_pwwi_fdr_1 <- format_wi_kw(data = L.mid_M00_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))
L.dtl_M00_pwwi_fdr_1 <- format_wi_kw(data = L.dtl_M00_pwwi_fdr_0, 
                                     digits = 3, cutpoints = cutpts, 
                                     symbols = symbls, 
                                     names_P = c("BH_P", "sig_BH_P"))

L.tot_F01_pwkw <- merge(x = L.tot_F01_pwkw_raw_1, y = L.tot_F01_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.prx_F01_pwkw <- merge(x = L.prx_F01_pwkw_raw_1, y = L.prx_F01_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.mid_F01_pwkw <- merge(x = L.mid_F01_pwkw_raw_1, y = L.mid_F01_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.dtl_F01_pwkw <- merge(x = L.dtl_F01_pwkw_raw_1, y = L.dtl_F01_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.tot_M00_pwkw <- merge(x = L.tot_M00_pwkw_raw_1, y = L.tot_M00_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.prx_M00_pwkw <- merge(x = L.prx_M00_pwkw_raw_1, y = L.prx_M00_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.mid_M00_pwkw <- merge(x = L.mid_M00_pwkw_raw_1, y = L.mid_M00_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.dtl_M00_pwkw <- merge(x = L.dtl_M00_pwkw_raw_1, y = L.dtl_M00_pwkw_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.tot_F01_pwwi <- merge(x = L.tot_F01_pwwi_raw_1, y = L.tot_F01_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.prx_F01_pwwi <- merge(x = L.prx_F01_pwwi_raw_1, y = L.prx_F01_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.mid_F01_pwwi <- merge(x = L.mid_F01_pwwi_raw_1, y = L.mid_F01_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.dtl_F01_pwwi <- merge(x = L.dtl_F01_pwwi_raw_1, y = L.dtl_F01_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.tot_M00_pwwi <- merge(x = L.tot_M00_pwwi_raw_1, y = L.tot_M00_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.prx_M00_pwwi <- merge(x = L.prx_M00_pwwi_raw_1, y = L.prx_M00_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.mid_M00_pwwi <- merge(x = L.mid_M00_pwwi_raw_1, y = L.mid_M00_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))
L.dtl_M00_pwwi <- merge(x = L.dtl_M00_pwwi_raw_1, y = L.dtl_M00_pwwi_fdr_1, 
                        sort = F, by = c("Test", "Comparison"))

L.tot_F01_pwse <- rbind(L.tot_F01_pwkw, L.tot_F01_pwwi)
L.prx_F01_pwse <- rbind(L.prx_F01_pwkw, L.prx_F01_pwwi)
L.mid_F01_pwse <- rbind(L.mid_F01_pwkw, L.mid_F01_pwwi)
L.dtl_F01_pwse <- rbind(L.dtl_F01_pwkw, L.dtl_F01_pwwi)
L.tot_M00_pwse <- rbind(L.tot_M00_pwkw, L.tot_M00_pwwi)
L.prx_M00_pwse <- rbind(L.prx_M00_pwkw, L.prx_M00_pwwi)
L.mid_M00_pwse <- rbind(L.mid_M00_pwkw, L.mid_M00_pwwi)
L.dtl_M00_pwse <- rbind(L.dtl_M00_pwkw, L.dtl_M00_pwwi)

### ************************************
### L - STEP  4 - combine global and pairwise results ----
### ************************************

# combine results data.frames as follows:
# (1) combine global and pwise results by Data type
# (2;3) add new column specifying Data type
# (4) combine results together

L.tot_F01_rslt_0 <- rbind(L.tot_F01_glbl, L.tot_F01_pwse)
L.prx_F01_rslt_0 <- rbind(L.prx_F01_glbl, L.prx_F01_pwse)
L.mid_F01_rslt_0 <- rbind(L.mid_F01_glbl, L.mid_F01_pwse)
L.dtl_F01_rslt_0 <- rbind(L.dtl_F01_glbl, L.dtl_F01_pwse)
L.tot_M00_rslt_0 <- rbind(L.tot_M00_glbl, L.tot_M00_pwse)
L.prx_M00_rslt_0 <- rbind(L.prx_M00_glbl, L.prx_M00_pwse)
L.mid_M00_rslt_0 <- rbind(L.mid_M00_glbl, L.mid_M00_pwse)
L.dtl_M00_rslt_0 <- rbind(L.dtl_M00_glbl, L.dtl_M00_pwse)

L.tot_F01_rslt <- L.tot_F01_rslt_0
L.prx_F01_rslt <- L.prx_F01_rslt_0
L.mid_F01_rslt <- L.mid_F01_rslt_0
L.dtl_F01_rslt <- L.dtl_F01_rslt_0
L.tot_M00_rslt <- L.tot_M00_rslt_0
L.prx_M00_rslt <- L.prx_M00_rslt_0
L.mid_M00_rslt <- L.mid_M00_rslt_0
L.dtl_M00_rslt <- L.dtl_M00_rslt_0
L.tot_F01_rslt$Data <- "total lesions"
L.prx_F01_rslt$Data <- "proximal lesions"
L.mid_F01_rslt$Data <- "middle lesions"
L.dtl_F01_rslt$Data <- "distal lesions"
L.tot_M00_rslt$Data <- "total lesions"
L.prx_M00_rslt$Data <- "proximal lesions"
L.mid_M00_rslt$Data <- "middle lesions"
L.dtl_M00_rslt$Data <- "distal lesions"

L.rslt_F01_0 <- rbind(L.tot_F01_rslt, L.prx_F01_rslt, L.mid_F01_rslt, 
                      L.dtl_F01_rslt)
L.rslt_M00_0 <- rbind(L.tot_M00_rslt, L.prx_M00_rslt, L.mid_M00_rslt, 
                      L.dtl_M00_rslt)

L.rslt_F01 <- L.rslt_F01_0
L.rslt_M00 <- L.rslt_M00_0
L.rslt_F01$Cohort <- "human female"
L.rslt_M00$Cohort <- "humane male"

L.rslt <- rbind(L.rslt_F01, L.rslt_M00)

### ************************************
### L - ^^^^STEP  5 - define plot parameters/customize plot aesthetics ----
### ************************************

# ^^^^ look at comments for when consortium group is discussed vs murine-sex

### consortium-based plots
## universal to both plot types:
# define comparisons for statistical tests used in plots
L.comps_grp <- list(c("N", "CMC"), c("CMC", "RMC"), c("N", "RMC"))

L.comps_alt <- list(c("N.1", "CMC.1"), c("CMC.1", "RMC.1"), c("N.1", "RMC.1"),
                    c("N.2", "CMC.2"), c("CMC.2", "RMC.2"), c("N.2", "RMC.2"),
                    c("N.3", "CMC.3"), c("CMC.3", "RMC.3"), c("N.3", "RMC.3"))

# x axis order
L.xord <- c("N", "CMC", "RMC")
L.xord_alt <- c("N.1", "CMC.1", "RMC.1", "N.2", "CMC.2", "RMC.2", 
                "N.3", "CMC.3", "RMC.3")

# y axis limits, breaks, labels
L.ylim <- c(0, 9.5)
L.ybrk <- c(0, 1, 2, 3, 4, 5, 6, 7)
L.ylab <- c("0", "", "2", "", "4", "", "6", "")

# plot labeling
L.ttl_prx <- "proximal"
L.ttl_mid <- "middle"
L.ttl_dtl <- "distal"

# ^^^create a data.frame to plot dotted lines and text labels for location
L.lab_tot <- data.frame(x = 2, y = L.ylim[2], l = "total", stringsAsFactors = F)

L.lne_loc <- data.frame(x = c(3.5, 6.5), ylne = L.ylim[2], stringsAsFactors = F)

L.lab_loc <- data.frame(x = c(2, 5, 8), y = L.ylim[2], 
                        l = c("proximal", "middle", "distal"), 
                        stringsAsFactors = F)

## specific to barplots:
# x axis labels
# L.xlab_bar <- c("No\ninoculum\n(N)\n \n", 
#                 "Control\nmicrobial\nconsortium\n(CMC)\n ", 
#                 "Rice bran\nmodified\nmicrobial\nconsortium\n(RMC)")
L.xlab_bar_F01 <- c("N", "CMC-f", "RMC-f")
L.xlab_bar_M00 <- c("N", "CMC-m", "RMC-m")

# vector for bar fill colors by consortium group
L.bar_fil_F01 <- c(hex_lyt_nmc, hex_F01_cmc[2], hex_F01_rmc[2])
L.bar_fil_M00 <- c(hex_lyt_nmc, hex_M00_cmc[1], hex_M00_rmc[2])

## specific to dotplots:
# x axis labels
# L.xlab_dot <. c("N\n \n \n \n", "CMC\n \n \n \n", "RMC\n \n \n \n")
L.xlab_dot_F01 <- c("N", "CMC-f", "RMC-f", 
                    "N", "CMC-f", "RMC-f", 
                    "N", "CMC-f", "RMC-f")
L.xlab_dot_M00 <- c("N", "CMC-m", "RMC-m", 
                    "N", "CMC-m", "RMC-m", 
                    "N", "CMC-m", "RMC-m")

# vector for shape outline colors by consortium group
L.shp_hex_F01_loc <- c(hex_drk_nmc, hex_F01_cmc[2], hex_F01_rmc[2],
                       hex_drk_nmc, hex_F01_cmc[1], hex_F01_rmc[1],
                       hex_drk_nmc, hex_F01_cmc[1], hex_F01_rmc[1])
L.shp_hex_M00_loc <- c(hex_drk_nmc, hex_M00_cmc[2], hex_M00_rmc[2],
                       hex_drk_nmc, hex_M00_cmc[1], hex_M00_rmc[1],
                       hex_drk_nmc, hex_M00_cmc[1], hex_M00_rmc[1])

# vector for shape fill colors by consortium group
L.shp_fil_F01_loc <- c(hex_lyt_nmc, hex_F01_cmc[5], hex_F01_rmc[5],
                       hex_lyt_nmc, hex_F01_cmc[3], hex_F01_rmc[3],
                       hex_lyt_nmc, hex_F01_cmc[2], hex_F01_rmc[2])

L.shp_fil_M00_loc <- c(hex_lyt_nmc, hex_M00_cmc[4], hex_M00_rmc[5],
                       hex_lyt_nmc, hex_M00_cmc[3], hex_M00_rmc[3],
                       hex_lyt_nmc, hex_M00_cmc[2], hex_M00_rmc[2])

### ^^^^murine sex-based plots
L.xord_sex_F01 <- c("N.Male", "CMC.Female", "RMC.Female", "RMC.Male")
L.xlab_sex_bar_F01 <- c("N\nMale\nmice", 
                        "CMC-f\nFemale\nmice", 
                        "RMC-f\nMale\nmice", "RMC-f\nFemale\nmice")
L.xord_sex_M00 <- c("N.Male", "CMC.Female", "CMC.Male", "RMC.Female", 
                    "RMC.Male")
L.xlab_sex_bar_M00 <- c("N\nMale\nmice", 
                        "CMC-m\nFemale\nmice", "CMC-m\nMale\nmice", 
                        "RMC-m\nMale\nmice", "RMC-m\nFemale\nmice")
L.bar_fil_sex_F01 <- c(hex_F01_cmc[2], hex_lyt_nmc, hex_F01_rmc[2])
L.bar_fil_sex_M00 <- c(hex_M00_cmc[2], hex_lyt_nmc, hex_M00_rmc[2])
L.comps_sex_F01 <- list(c("RMC.Female", "RMC.Male"))
L.comps_sex_M00 <- list(c("CMC.Female", "CMC.Male"),
                        c("RMC.Female", "RMC.Male"))
L.non_sig_sex_F01 <- data.frame(x = 3.5, y = 3.33, label = "ns", 
                                stringsAsFactors = F)
L.val_sig_sex_M00 <- data.frame(x = c(2.5, 4.5), y = c(7.5, 3.33),
                                label = c("*", "ns"), stringsAsFactors = F)
L.laby_sex_F01 <- 3.33
L.laby_sex_M00 <- c(7.5, 3.33)

# sizing parameters
L.sze_sts_lab <- 2.22 # (all) size for p value significance labels
L.sze_tip_len <- 0.00 # (all) size for bracket tip length
L.sze_bkt_lne <- 0.50 # (all) size for bracket lines
L.sze_loc_txt <- 2.88 # (all) - size for sample location text
L.sze_axs_txt <- 08 # (all) - size for axis text
L.sze_plt_ttl <- 10 # (all) - size for plot title
L.sze_bar_wid <- 0.51 # (bar) size for bar width
L.sze_bin_wid <- 0.27 # (dot) size for bin width
L.sze_pts_dot <- 0.81 # (dot) size for points
L.sze_loc_lne <- 0.21 # (dot) - size for sample location separator lines

L.sclr_sts <- 0.25 # label scaler for asterisks

L.sze_sex_sts_lab <- 1.51 # size for p value significance labels
L.sze_sex_tip_len <- 0.00 # size for bracket tip length
L.sze_sex_bkt_lne <- 0.42 # size for bracket lines
L.sze_sex_axs_txt <- 4 # size for axis text
L.sze_sex_bar_wid <- 0.42 # (bar) size for bar width

### ************************************
### L - STEP #a - plotting: ----
### ************************************

# NOTE: the stat_compare_means() function fails to change the font family ...
# ... when the comparisons argument is specified, so we will 'turn off' the ...
# ... labels for stat_compare_means() use the results generated in STEP 3 above
# NOTE: proceeding this way also allows us to plot the BH/fdr adjusted sigs ...
# ... rather than the undajusted values used by stat_compare_means()

# ^^^define x and y locations for each comparison
L.x_CvN <- 1.5
L.x_RvN <- 2.0
L.x_RvC <- 2.5

L.y_CvN <- 6.75
L.y_RvN <- 8.25
L.y_RvC <- 7.50

# ^^^create copies, then add col specifying blerh and combine prx, mid, dtl
L.sig_tot_F01_0 <- L.tot_F01_pwwi
L.sig_prx_F01_pwwi <- L.prx_F01_pwwi
L.sig_mid_F01_pwwi <- L.mid_F01_pwwi
L.sig_dtl_F01_pwwi <- L.dtl_F01_pwwi
L.sig_tot_M00_0 <- L.tot_M00_pwwi
L.sig_prx_M00_pwwi <- L.prx_M00_pwwi
L.sig_mid_M00_pwwi <- L.mid_M00_pwwi
L.sig_dtl_M00_pwwi <- L.dtl_M00_pwwi

L.sig_tot_F01_0$loc <- "tot"
L.sig_prx_F01_pwwi$loc <- "prx"
L.sig_mid_F01_pwwi$loc <- "mid"
L.sig_dtl_F01_pwwi$loc <- "dtl"
L.sig_tot_M00_0$loc <- "tot"
L.sig_prx_M00_pwwi$loc <- "prx"
L.sig_mid_M00_pwwi$loc <- "mid"
L.sig_dtl_M00_pwwi$loc <- "dtl"

# combine prx, mid, dtl
L.sig_loc_F01_0 <- rbind(L.sig_prx_F01_pwwi, L.sig_mid_F01_pwwi, 
                         L.sig_dtl_F01_pwwi)
L.sig_loc_M00_0 <- rbind(L.sig_prx_M00_pwwi, L.sig_mid_M00_pwwi, 
                         L.sig_dtl_M00_pwwi)

# ^^^create copies then add x and y locations and values for ordering the dfs
L.sig_tot_F01_1 <- L.sig_tot_F01_0
L.sig_loc_F01_1 <- L.sig_loc_F01_0
L.sig_tot_M00_1 <- L.sig_tot_M00_0
L.sig_loc_M00_1 <- L.sig_loc_M00_0

L.sig_tot_F01_1$x[L.sig_tot_F01_1$Comparison == "CMC vs N"] <- L.x_CvN
L.sig_tot_F01_1$x[L.sig_tot_F01_1$Comparison == "RMC vs N"] <- L.x_RvN
L.sig_tot_F01_1$x[L.sig_tot_F01_1$Comparison == "RMC vs CMC"] <- L.x_RvC
L.sig_tot_F01_1$y[L.sig_tot_F01_1$Comparison == "CMC vs N"] <- L.y_CvN
L.sig_tot_F01_1$y[L.sig_tot_F01_1$Comparison == "RMC vs N"] <- L.y_RvN
L.sig_tot_F01_1$y[L.sig_tot_F01_1$Comparison == "RMC vs CMC"] <- L.y_RvC
L.sig_tot_F01_1$ordr[L.sig_tot_F01_1$Comparison == "CMC vs N"] <- 1
L.sig_tot_F01_1$ordr[L.sig_tot_F01_1$Comparison == "RMC vs N"] <- 3
L.sig_tot_F01_1$ordr[L.sig_tot_F01_1$Comparison == "RMC vs CMC"] <- 2
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "CMC vs N" & 
                    L.sig_loc_F01_1$loc == "prx"] <- L.x_CvN
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "RMC vs N" & 
                    L.sig_loc_F01_1$loc == "prx"] <- L.x_RvN
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "RMC vs CMC" & 
                    L.sig_loc_F01_1$loc == "prx"] <- L.x_RvC
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "CMC vs N" & 
                    L.sig_loc_F01_1$loc == "mid"] <- L.x_CvN + 3
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "RMC vs N" & 
                    L.sig_loc_F01_1$loc == "mid"] <- L.x_RvN + 3
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "RMC vs CMC" & 
                    L.sig_loc_F01_1$loc == "mid"] <- L.x_RvC + 3
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "CMC vs N" & 
                    L.sig_loc_F01_1$loc == "dtl"] <- L.x_CvN + 6
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "RMC vs N" & 
                    L.sig_loc_F01_1$loc == "dtl"] <- L.x_RvN + 6
L.sig_loc_F01_1$x[L.sig_loc_F01_1$Comparison == "RMC vs CMC" & 
                    L.sig_loc_F01_1$loc == "dtl"] <- L.x_RvC + 6
L.sig_loc_F01_1$y[L.sig_loc_F01_1$Comparison == "CMC vs N"] <- L.y_CvN
L.sig_loc_F01_1$y[L.sig_loc_F01_1$Comparison == "RMC vs N"] <- L.y_RvN
L.sig_loc_F01_1$y[L.sig_loc_F01_1$Comparison == "RMC vs CMC"] <- L.y_RvC
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "CMC vs N" & 
                       L.sig_loc_F01_1$loc == "prx"] <- 1
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "RMC vs N" & 
                       L.sig_loc_F01_1$loc == "prx"] <- 3
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "RMC vs CMC" & 
                       L.sig_loc_F01_1$loc == "prx"] <- 2
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "CMC vs N" & 
                       L.sig_loc_F01_1$loc == "mid"] <- 4
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "RMC vs N" & 
                       L.sig_loc_F01_1$loc == "mid"] <- 6
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "RMC vs CMC" & 
                       L.sig_loc_F01_1$loc == "mid"] <- 5
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "CMC vs N" & 
                       L.sig_loc_F01_1$loc == "dtl"] <- 7
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "RMC vs N" & 
                       L.sig_loc_F01_1$loc == "dtl"] <- 9
L.sig_loc_F01_1$ordr[L.sig_loc_F01_1$Comparison == "RMC vs CMC" & 
                       L.sig_loc_F01_1$loc == "dtl"] <- 8

L.sig_tot_M00_1$x[L.sig_tot_M00_1$Comparison == "CMC vs N"] <- L.x_CvN
L.sig_tot_M00_1$x[L.sig_tot_M00_1$Comparison == "RMC vs N"] <- L.x_RvN
L.sig_tot_M00_1$x[L.sig_tot_M00_1$Comparison == "RMC vs CMC"] <- L.x_RvC
L.sig_tot_M00_1$y[L.sig_tot_M00_1$Comparison == "CMC vs N"] <- L.y_CvN
L.sig_tot_M00_1$y[L.sig_tot_M00_1$Comparison == "RMC vs N"] <- L.y_RvN
L.sig_tot_M00_1$y[L.sig_tot_M00_1$Comparison == "RMC vs CMC"] <- L.y_RvC
L.sig_tot_M00_1$ordr[L.sig_tot_M00_1$Comparison == "CMC vs N"] <- 1
L.sig_tot_M00_1$ordr[L.sig_tot_M00_1$Comparison == "RMC vs N"] <- 3
L.sig_tot_M00_1$ordr[L.sig_tot_M00_1$Comparison == "RMC vs CMC"] <- 2
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "CMC vs N" & 
                    L.sig_loc_M00_1$loc == "prx"] <- L.x_CvN
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "RMC vs N" & 
                    L.sig_loc_M00_1$loc == "prx"] <- L.x_RvN
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "RMC vs CMC" & 
                    L.sig_loc_M00_1$loc == "prx"] <- L.x_RvC
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "CMC vs N" & 
                    L.sig_loc_M00_1$loc == "mid"] <- L.x_CvN + 3
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "RMC vs N" & 
                    L.sig_loc_M00_1$loc == "mid"] <- L.x_RvN + 3
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "RMC vs CMC" & 
                    L.sig_loc_M00_1$loc == "mid"] <- L.x_RvC + 3
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "CMC vs N" & 
                    L.sig_loc_M00_1$loc == "dtl"] <- L.x_CvN + 6
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "RMC vs N" & 
                    L.sig_loc_M00_1$loc == "dtl"] <- L.x_RvN + 6
L.sig_loc_M00_1$x[L.sig_loc_M00_1$Comparison == "RMC vs CMC" & 
                    L.sig_loc_M00_1$loc == "dtl"] <- L.x_RvC + 6
L.sig_loc_M00_1$y[L.sig_loc_M00_1$Comparison == "CMC vs N"] <- L.y_CvN
L.sig_loc_M00_1$y[L.sig_loc_M00_1$Comparison == "RMC vs N"] <- L.y_RvN
L.sig_loc_M00_1$y[L.sig_loc_M00_1$Comparison == "RMC vs CMC"] <- L.y_RvC
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "CMC vs N" & 
                       L.sig_loc_M00_1$loc == "prx"] <- 1
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "RMC vs N" & 
                       L.sig_loc_M00_1$loc == "prx"] <- 3
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "RMC vs CMC" & 
                       L.sig_loc_M00_1$loc == "prx"] <- 2
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "CMC vs N" & 
                       L.sig_loc_M00_1$loc == "mid"] <- 4
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "RMC vs N" & 
                       L.sig_loc_M00_1$loc == "mid"] <- 6
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "RMC vs CMC" & 
                       L.sig_loc_M00_1$loc == "mid"] <- 5
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "CMC vs N" & 
                       L.sig_loc_M00_1$loc == "dtl"] <- 7
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "RMC vs N" & 
                       L.sig_loc_M00_1$loc == "dtl"] <- 9
L.sig_loc_M00_1$ordr[L.sig_loc_M00_1$Comparison == "RMC vs CMC" & 
                       L.sig_loc_M00_1$loc == "dtl"] <- 8

# order dfs col ordr defined above
L.sig_tot_F01_2 <- L.sig_tot_F01_1[order(L.sig_tot_F01_1$ordr), ]
L.sig_loc_F01_2 <- L.sig_loc_F01_1[order(L.sig_loc_F01_1$ordr), ]
L.sig_tot_M00_2 <- L.sig_tot_M00_1[order(L.sig_tot_M00_1$ordr), ]
L.sig_loc_M00_2 <- L.sig_loc_M00_1[order(L.sig_loc_M00_1$ordr), ]

# replace empty values in sigs col with ns for non-significant
L.sig_tot_F01 <- L.sig_tot_F01_2
L.sig_loc_F01 <- L.sig_loc_F01_2
L.sig_tot_M00 <- L.sig_tot_M00_2
L.sig_loc_M00 <- L.sig_loc_M00_2

L.sig_tot_F01$sig_BH_P[L.sig_tot_F01$sig_BH_P == " "] <- "ns"
L.sig_loc_F01$sig_BH_P[L.sig_loc_F01$sig_BH_P == " "] <- "ns"
L.sig_tot_M00$sig_BH_P[L.sig_tot_M00$sig_BH_P == " "] <- "ns"
L.sig_loc_M00$sig_BH_P[L.sig_loc_M00$sig_BH_P == " "] <- "ns"

# barplot ~ total ----

L.gpr_bar_tot_F01 <- ggbarplot(data = L.tot_F01, y = "LesionTotal", 
                               x = "ConsortiumAbrv", color = greydient[1],
                               fill = "ConsortiumAbrv", font.family = fnt, 
                               order = L.xord, width = L.sze_bar_wid, 
                               add = "mean_se", error.plot = "upper_errorbar") +
  stat_compare_means(inherit.aes = T, comparison = L.comps_grp, geom = "text",
                     method = "wilcox.test", size = 0, 
                     bracket.size = L.sze_bkt_lne, tip.length = L.sze_tip_len, 
                     label.y = L.sig_tot_F01$y) +
  geom_text(data = L.sig_tot_F01, family = fnt, color = greydient[2], 
            size = L.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + L.sclr_sts, label = sig_BH_P)) +
  geom_text(data = L.lab_tot, aes(x = x, y = y, label = l), family = fnt,
            color = greydient[2], size = L.sze_loc_txt, inherit.aes = F) +
  scale_fill_manual(values = L.bar_fil_F01) +
  scale_x_discrete(labels = L.xlab_bar_F01) +
  scale_y_continuous(limits = L.ylim, breaks = L.ybrk, labels = L.ylab) +
  theme(text = element_text(family = fnt, color = greydient[1]),
        line = element_line(lineend = "square", color = greydient[1]),
        axis.text = element_text(size = L.sze_axs_txt),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        panel.ontop = T,
        panel.background = element_rect(fill = NA))

L.gpr_bar_tot_M00 <- ggbarplot(data = L.tot_M00, y = "LesionTotal", 
                               x = "ConsortiumAbrv", color = greydient[1],
                               fill = "ConsortiumAbrv", font.family = fnt, 
                               order = L.xord, width = L.sze_bar_wid, 
                               add = "mean_se", error.plot = "upper_errorbar") +
  stat_compare_means(inherit.aes = T, comparison = L.comps_grp, geom = "text",
                     method = "wilcox.test", size = 0,
                     bracket.size = L.sze_bkt_lne, tip.length = L.sze_tip_len, 
                     label.y = L.sig_tot_M00$y) +
  geom_text(data = L.sig_tot_M00, family = fnt, color = greydient[2], 
            size = L.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + L.sclr_sts, label = sig_BH_P)) +
  geom_text(data = L.lab_tot, aes(x = x, y = y, label = l), family = fnt,
            color = greydient[2], size = L.sze_loc_txt, inherit.aes = F) +
  scale_fill_manual(values = L.bar_fil_M00) +
  scale_x_discrete(labels = L.xlab_bar_M00) +
  scale_y_continuous(limits = L.ylim, breaks = L.ybrk, labels = L.ylab) +
  theme(text = element_text(family = fnt, color = greydient[1]),
        line = element_line(lineend = "square", color = greydient[1]),
        axis.text = element_text(size = L.sze_axs_txt),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        panel.ontop = T,
        panel.background = element_rect(fill = NA))

# dotplot ----

L.gpr_dot_loc_F01 <- ggdotplot(data = L.loc_F01, y = "LesionTot", 
                               x = "ConsortiumAbrvAlt", 
                               color = "ConsortiumAbrvAlt", 
                               fill = "ConsortiumAbrvAlt", 
                               size = L.sze_pts_dot, font.family = fnt,
                               order = L.xord_alt, binwidth = L.sze_bin_wid) +
  stat_compare_means(inherit.aes = T, comparison = L.comps_alt, geom = "text",
                     method = "wilcox.test", size = 0, 
                     bracket.size = L.sze_bkt_lne, tip.length = L.sze_tip_len,
                     label.y = L.sig_loc_F01$y) +
  geom_text(data = L.sig_loc_F01, family = fnt, color = greydient[2],
            size = L.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + L.sclr_sts, label = sig_BH_P)) +
  geom_segment(inherit.aes = F, data = L.lne_loc,
               aes(x = x, xend = x, y = ylne, yend = 0), linetype = "dotted",
               color = greydient[6], size = L.sze_loc_lne) +
  geom_text(data = L.lab_loc, aes(x = x, y = y, label = l), family = fnt,
            color = greydient[2], size = L.sze_loc_txt, inherit.aes = F) +
  scale_color_manual(values = L.shp_hex_F01_loc) +
  scale_fill_manual(values = L.shp_fil_F01_loc) +
  scale_x_discrete(labels = L.xlab_dot_F01) +
  scale_y_continuous(limits = L.ylim, breaks = L.ybrk, labels = L.ylab) +
  theme(text = element_text(family = fnt, color = greydient[1]),
        line = element_line(lineend = "square", color = greydient[1]),
        axis.text = element_text(size = L.sze_axs_txt),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        panel.ontop = T,
        panel.background = element_rect(fill = NA))

L.gpr_dot_loc_M00 <- ggdotplot(data = L.loc_M00, y = "LesionTot", 
                               x = "ConsortiumAbrvAlt", 
                               color = "ConsortiumAbrvAlt", 
                               fill = "ConsortiumAbrvAlt", 
                               size = L.sze_pts_dot, font.family = fnt,
                               order = L.xord_alt, binwidth = L.sze_bin_wid) +
  stat_compare_means(inherit.aes = T, comparison = L.comps_alt, geom = "text",
                     method = "wilcox.test", size = 0, 
                     bracket.size = L.sze_bkt_lne, tip.length = L.sze_tip_len,
                     label.y = L.sig_loc_M00$y) +
  geom_text(data = L.sig_loc_M00, family = fnt, color = greydient[2],
            size = L.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + L.sclr_sts, label = sig_BH_P)) +
  geom_segment(inherit.aes = F, data = L.lne_loc,
               aes(x = x, xend = x, y = ylne, yend = 0), linetype = "dotted",
               color = greydient[6], size = L.sze_loc_lne) +
  geom_text(data = L.lab_loc, aes(x = x, y = y, label = l), family = fnt,
            color = greydient[2], size = L.sze_loc_txt, inherit.aes = F) +
  scale_color_manual(values = L.shp_hex_M00_loc) +
  scale_fill_manual(values = L.shp_fil_M00_loc) +
  scale_x_discrete(labels = L.xlab_dot_M00) +
  scale_y_continuous(limits = L.ylim, breaks = L.ybrk, labels = L.ylab) +
  theme(text = element_text(family = fnt, color = greydient[1]),
        line = element_line(lineend = "square", color = greydient[1]),
        axis.text = element_text(size = L.sze_axs_txt),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        panel.ontop = T,
        panel.background = element_rect(fill = NA))

# murine sex-based plots ----

# F01
L.gpr_bar_sex_F01 <- ggbarplot(data = L.tot_sex_F01, y = "LesionTotal",
                               x = "ConsortiumAnimalSex", color = greydient[1],
                               fill = "ConsortiumAbrv", font.family = fnt, 
                               order = L.xord_sex_F01, width = L.sze_sex_bar_wid,
                               add = "mean_se", error.plot = "upper_errorbar") +
  stat_compare_means(inherit.aes = T, comparison = L.comps_sex_F01, geom = "text",
                     method = "wilcox.test", 
                     size = 0, 
                     label.y = L.laby_sex_F01,
                     bracket.size = L.sze_sex_bkt_lne, 
                     tip.length = L.sze_sex_tip_len, 
                     hide.ns = T, 
                     na.rm = T) +
  geom_text(data = L.non_sig_sex_F01, aes(x = x, y = y + L.sclr_sts, label = label), 
            family = fnt, color = greydient[2], size = L.sze_sex_sts_lab, 
            inherit.aes = F) +
  scale_fill_manual(values = L.bar_fil_sex_F01) +
  scale_x_discrete(labels = L.xlab_sex_bar_F01) +
  scale_y_continuous(limits = L.ylim, breaks = L.ybrk, labels = L.ylab) +
  theme(text = element_text(family = fnt, color = greydient[1]),
        line = element_line(lineend = "square", color = greydient[1]),
        axis.text = element_text(size = L.sze_sex_axs_txt),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        panel.ontop = T,
        panel.background = element_rect(fill = NA))

# M00 (P = 0.019 for CMC comps)
L.gpr_bar_sex_M00 <- ggbarplot(data = L.tot_sex_M00, y = "LesionTotal",
                               x = "ConsortiumAnimalSex", color = greydient[1],
                               fill = "ConsortiumAbrv", font.family = fnt, 
                               order = L.xord_sex_M00, width = L.sze_sex_bar_wid,
                               add = "mean_se", error.plot = "upper_errorbar") +
  stat_compare_means(inherit.aes = T, comparison = L.comps_sex_M00, geom = "text",
                     method = "wilcox.test",
                     size = 0,
                     # size = 2, label = "p.format",
                     label.y = L.laby_sex_M00,
                     bracket.size = L.sze_sex_bkt_lne, 
                     tip.length = L.sze_sex_tip_len) +
  geom_text(data = L.val_sig_sex_M00, aes(x = x, y = y + L.sclr_sts, label = label),
            family = fnt, color = greydient[2], size = L.sze_sex_sts_lab,
            inherit.aes = F) +
  scale_fill_manual(values = L.bar_fil_sex_M00) +
  scale_x_discrete(labels = L.xlab_sex_bar_M00) +
  scale_y_continuous(limits = L.ylim, breaks = L.ybrk, labels = L.ylab) +
  theme(text = element_text(family = fnt, color = greydient[1]),
        line = element_line(lineend = "square", color = greydient[1]),
        axis.text = element_text(size = L.sze_sex_axs_txt),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        panel.ontop = T,
        panel.background = element_rect(fill = NA))
L.gpr_bar_sex_M00
### ************************************
### L - STEP #b - plotting: arrange ----
### ************************************

# arrange by sex of human donor
L.gga_F01_0 <- ggarrange(L.gpr_bar_tot_F01, L.gpr_dot_loc_F01, 
                         labels = c("B", ""), font.label = pan_fnt,
                         ncol = 2, nrow = 1, widths = c(1, 3))
L.gga_M00_0 <- ggarrange(L.gpr_bar_tot_M00, L.gpr_dot_loc_M00, 
                         labels = c("C", ""), font.label = pan_fnt,
                         ncol = 2, nrow = 1, widths = c(1, 3))

# annotate arrangement with a universal y-axis title
L.yttl <- text_grob("number of lesions/colon", color = greydient[1], size = 08,
                    family = fnt, rot = 90)

L.gga_F01 <- annotate_figure(L.gga_F01_0, left = L.yttl)
L.gga_M00 <- annotate_figure(L.gga_M00_0, left = L.yttl)

# arrange arrangements
L.gga <- ggarrange(L.gga_F01, L.gga_M00, labels = NULL, ncol = 1, nrow = 2)

# murine sex-based plots
# L.gga_sex <- ggarrange(L.gpr_bar_sex_F01, L.gpr_bar_sex_M00, 
#                        labels = c("A", "D"), font.label = pan_fnt,
#                        ncol = 1, nrow = 2)

L.gga_bar_F01_sex <- ggarrange(L.gpr_bar_sex_F01, 
                               labels = "A", font.label = pan_fnt,
                               ncol = 1, nrow = 1)
L.gga_bar_M00_sex <- ggarrange(L.gpr_bar_sex_M00, 
                               labels = "D", font.label = pan_fnt,
                               ncol = 1, nrow = 1)

### ************************************
### L - WRITE OUTPUTS ----
### ************************************

# total lesions F01
L.ofv_gpr_bar_tot_F01 <- "TransFaunation/vault/plot_bar_lesions_F01_total.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 55, height = 55,
       filename = L.ofv_gpr_bar_tot_F01, plot = L.gpr_bar_tot_F01)

# distribution of lesions F01
L.ofv_gpr_dot_loc_F01 <- "TransFaunation/vault/plot_bar_lesions_F01_distr.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 110, height = 55,
       filename = L.ofv_gpr_dot_loc_F01, plot = L.gpr_dot_loc_F01)

# arrangement of total and distribution of lesions F01
L.ofv_gga_F01 <- "TransFaunation/vault/plot_bar_lesions_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = L.ofv_gga_F01, plot = L.gga_F01)

# total lesions M00
L.ofv_gpr_bar_tot_M00 <- "TransFaunation/vault/plot_bar_lesions_M00_total.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 55, height = 55,
       filename = L.ofv_gpr_bar_tot_M00, plot = L.gpr_bar_tot_M00)

# distribution of lesions M00
L.ofv_gpr_dot_loc_M00 <- "TransFaunation/vault/plot_bar_lesions_M00_distr.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 110, height = 55,
       filename = L.ofv_gpr_dot_loc_M00, plot = L.gpr_dot_loc_M00)

# arrangement of total and distribution of lesions M00
L.ofv_gga_M00 <- "TransFaunation/vault/plot_bar_lesions_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 55,
       filename = L.ofv_gga_M00, plot = L.gga_M00)

# arranged arrangements
L.ofv_plot <- "TransFaunation/vault/plot_bar_lesions.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 110,
       filename = L.ofv_plot, plot = L.gga)

# new version of how to save these
L.ofv_plot_F01_sex <- "TransFaunation/vault/plot_bar_lesions_F01_animal_sex.pdf"
L.ofv_plot_M00_sex <- "TransFaunation/vault/plot_bar_lesions_M00_animal_sex.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 40, height = 45,
       filename = L.ofv_plot_F01_sex, plot = L.gga_bar_F01_sex)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 40, height = 45,
       filename = L.ofv_plot_M00_sex, plot = L.gga_bar_M00_sex)
# ^^REWRITE THIS NOTE: height for L.gga_sex must be equal to height of B.gga_sEPD_sex

L.ofv_stats_lesion <- "TransFaunation/vault/table_stats_lesion.txt"
write.table(sep = "\t", row.names = F, x = L.rslt, file = L.ofv_stats_lesion)

# save workspace
L.obj <- ls(pattern = "L.")
L.lst <- c(L.obj[grep(pattern = "L.", x = L.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS)
save(list = L.lst, file = L.ofv_wksp)

### ************************************
### BEGIN Section Y ----
### ************************************

# note for kdp: taken from TransFaunation/section_Y.R

### preface ----

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)

setwd("~/Desktop/")
load(file = "TransFaunation/vault/WS_Script02_PREFACE.RData")
load(file = ofv_COSMOS_wksp)

Y.ifp_mtb_cln <- "TransFaunation/mtbs/mtbs_colon.txt"
Y.ifp_mtb_plm <- "TransFaunation/mtbs/mtbs_plasma.txt"

Y.ofv_wksp <- paste(sep = "", path_vault, "/WS_", "SectionY.RData")

### ************************************
### Y - FUNCTION ----
### ************************************

# create a vector naming all Section A functions (used when saving workspace)
Y.function <- "iso_var_exp"

# ** note for KDP: version 0.2 ** #
# function takes an input of class "pcoa" or class "prcomp" and isolates...
# the amount of variation explained by the PC specified
# NOTE: isolates 1 value at a time as specified by pc.num
# NOTE: returns % values with 3 digits
# iso_var_exp = isolate_variation_explained
iso_var_exp <- function(data = c(pcoa, prcomp), pc.num, 
                        output = c("value", "label")) {
  # internal check to ensure correct input class
  if (!inherits(data, c("pcoa", "prcomp"))) {
    stop("data must be class 'pcoa' OR class 'prcomp'")
  }
  if (!inherits(pc.num, "numeric")) {
    stop("pc.num must be numeric")
  }
  if (!output == "value" && !output == "label") {
    stop("output must be one of: 'value' OR 'label'")
  }
  if (inherits(data, "pcoa")) {
    # isolate pcoa$values (create a data.frame)
    pcoa_vals <- data.frame(data$values)
    # isolate the Relative_eig values
    rel_eigs <- pcoa_vals$Relative_eig
    # convert rel eigs to percentage and round to 3 digits
    rel_eigs_rnd <- (round(rel_eigs, digits = 3)) * 100
    # isolate the specified PC
    pcoa_val <- rel_eigs_rnd[pc.num]
    # statements to determine what to return
    if (output == "value") {
      return(pcoa_val)
    }
    if (output == "label") {
      pcoa_lab <- paste("PC ", pc.num, " (", pcoa_val, "%)", sep = "")
      return(pcoa_lab)
    }
  }
  if (inherits(data, "prcomp")) {
    # isolate standard deviations for all pcs
    pcs_sdv <- data$sdev
    # square each sdev, sum the squares, and divide each sdev by the sum
    pcs_sdv_sqr <- pcs_sdv^2
    pcs_sdv_sqr_sum <- sum(pcs_sdv_sqr)
    pca_vals_raw <- pcs_sdv_sqr/pcs_sdv_sqr_sum
    # convert to percentage and round to 3 digits
    pca_vals_rnd <- (round(pca_vals_raw, digits = 3)) * 100
    # isolate the specified PC
    pca_val <- pca_vals_rnd[pc.num]
    # statements to determine what to return
    if (output == "value") {
      return(pca_val)
    }
    if (output == "label") {
      pca_lab <- paste("PC ", pc.num, " (", pca_val, "%)", sep = "")
      return(pca_lab)
    }
  }
}
# 
# # example usage:
# val <- iso_var_exp(data = pcoa, pc.num = 1, output = "value") # returns number
# lab <- iso_var_exp(data = pcoa, pc.num = 1, output = "label") # returns label

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

### ************************************
### Y - STEP  1 - format metabolite data for PCA & clustering ----
### ************************************

# this step formats the metabolite data in two parts:

# part 1 = replace NA's w/ 0 & subsets by sample type & sex of human donor:
# _cln = all colon samples for both sexes of human donors
# _prx_F01 = proximal colon; cmc & rmc from female
# _prx_M00 = proximal colon; cmc & rmc from male
# _dtl_F01 = distal colon; cmc & rmc from female
# _dtl_M00 = distal colon; cmc & rmc from male
# _plm_F01 = plasma; cmc & rmc from female
# _plm_M00 = plasma; cmc & rmc from male

# part 2 = remove any features with a total count of 0 across all samples ...
# ... and format to convert features to cols and samples to rows

# part 1 of this process occurs as follows:
# (1) read in the metabolite data
# (2) create a copy to avoid overwriting the originals
# (3) convert column FeatureID into row.names
# (4) remove unneeded columns
# (5;6) create a copy of the data.frame and replace any NA values with zero
# (7) filter to isolate appropriate samples
# NOTE: the mislabeled sample (T1_13_P_I) in the plm metabolite data is handled

Y.abs_mtb_cln_0 <- read.table(Y.ifp_mtb_cln, sep = "\t", header = T, as.is = T, 
                              stringsAsFactors = F, check.names = F, quote = "",
                              comment.char = "")
Y.abs_mtb_plm_0 <- read.table(Y.ifp_mtb_plm, sep = "\t", header = T, as.is = T, 
                              stringsAsFactors = F, check.names = F, quote = "",
                              comment.char = "")

Y.abs_mtb_cln_1 <- Y.abs_mtb_cln_0
Y.abs_mtb_plm_1 <- Y.abs_mtb_plm_0

row.names(Y.abs_mtb_cln_1) <- Y.abs_mtb_cln_1$FeatureID
row.names(Y.abs_mtb_plm_1) <- Y.abs_mtb_plm_1$FeatureID

Y.abs_mtb_cln_2 <- dplyr::select(Y.abs_mtb_cln_1, -FeatureID, -BIOCHEMICAL, 
                                 -SUPERPATHWAY, -SUBPATHWAY, -RI, -MASS, -CAS, 
                                 -PUBCHEM, -KEGG, -HMDB)
Y.abs_mtb_plm_2 <- dplyr::select(Y.abs_mtb_plm_1, -FeatureID, -BIOCHEMICAL, 
                                 -SUPERPATHWAY, -SUBPATHWAY, -RI, -MASS, -CAS, 
                                 -PUBCHEM, -KEGG, -HMDB)

Y.abs_mtb_cln_3 <- Y.abs_mtb_cln_2
Y.abs_mtb_plm_3 <- Y.abs_mtb_plm_2
Y.abs_mtb_cln_3[is.na(Y.abs_mtb_cln_3)] <- 0
Y.abs_mtb_plm_3[is.na(Y.abs_mtb_plm_3)] <- 0

Y.abs_mtb_cln_F01 <- dplyr::select(Y.abs_mtb_cln_3, 
                                   dplyr::contains("C1", ignore.case = F),
                                   dplyr::contains("T1", ignore.case = F))
Y.abs_mtb_cln_M00 <- dplyr::select(Y.abs_mtb_cln_3, 
                                   dplyr::contains("C2", ignore.case = F),
                                   dplyr::contains("T2", ignore.case = F))
Y.abs_mtb_prx_F01_0 <- dplyr::select(Y.abs_mtb_cln_F01, 
                                     dplyr::contains("CO_P", ignore.case = F))
Y.abs_mtb_prx_M00_0 <- dplyr::select(Y.abs_mtb_cln_M00, 
                                     dplyr::contains("CO_P", ignore.case = F))
Y.abs_mtb_dtl_F01_0 <- dplyr::select(Y.abs_mtb_cln_F01, 
                                     dplyr::contains("CO_D", ignore.case = F))
Y.abs_mtb_dtl_M00_0 <- dplyr::select(Y.abs_mtb_cln_M00, 
                                     dplyr::contains("CO_D", ignore.case = F))
Y.abs_mtb_plm_F01_0 <- dplyr::select(Y.abs_mtb_plm_3, 
                                     dplyr::contains("C1", ignore.case = F),
                                     dplyr::contains("T1", ignore.case = F),
                                     -T1_13_P_I)
Y.abs_mtb_plm_M00_0 <- dplyr::select(Y.abs_mtb_plm_3, 
                                     dplyr::contains("C2", ignore.case = F),
                                     dplyr::contains("T2", ignore.case = F),
                                     T1_13_P_I)

# part 2 of this process occurs as follows:
# (1) sum absolute counts for each feature across all samples
# (2) create a data.frame of features that will be removed
# (3) ensure that FeatureID is a character and not a factor
# (4;5) for original dfs, convert FeatureID back into a col
# (6) anti join the above to remove all features present in the rmv df
# (7;8;9) create copies; convert FeatureID into row.names; remove unneeded cols
# (10) transpose to convert features to cols and samples to rows

Y.mtb_prx_F01_sum <- data.frame("FeatureID" = row.names(Y.abs_mtb_prx_F01_0),
                                "FeatureTotal" = rowSums(Y.abs_mtb_prx_F01_0),
                                row.names = NULL)
Y.mtb_prx_M00_sum <- data.frame("FeatureID" = row.names(Y.abs_mtb_prx_M00_0),
                                "FeatureTotal" = rowSums(Y.abs_mtb_prx_M00_0),
                                row.names = NULL)
Y.mtb_dtl_F01_sum <- data.frame("FeatureID" = row.names(Y.abs_mtb_dtl_F01_0),
                                "FeatureTotal" = rowSums(Y.abs_mtb_dtl_F01_0),
                                row.names = NULL)
Y.mtb_dtl_M00_sum <- data.frame("FeatureID" = row.names(Y.abs_mtb_dtl_M00_0),
                                "FeatureTotal" = rowSums(Y.abs_mtb_dtl_M00_0),
                                row.names = NULL)
Y.mtb_plm_F01_sum <- data.frame("FeatureID" = row.names(Y.abs_mtb_plm_F01_0),
                                "FeatureTotal" = rowSums(Y.abs_mtb_plm_F01_0),
                                row.names = NULL)
Y.mtb_plm_M00_sum <- data.frame("FeatureID" = row.names(Y.abs_mtb_plm_M00_0),
                                "FeatureTotal" = rowSums(Y.abs_mtb_plm_M00_0),
                                row.names = NULL)

Y.mtb_prx_F01_rmv <- dplyr::filter(Y.mtb_prx_F01_sum, FeatureTotal < 1)
Y.mtb_prx_M00_rmv <- dplyr::filter(Y.mtb_prx_M00_sum, FeatureTotal < 1)
Y.mtb_dtl_F01_rmv <- dplyr::filter(Y.mtb_dtl_F01_sum, FeatureTotal < 1)
Y.mtb_dtl_M00_rmv <- dplyr::filter(Y.mtb_dtl_M00_sum, FeatureTotal < 1)
Y.mtb_plm_F01_rmv <- dplyr::filter(Y.mtb_plm_F01_sum, FeatureTotal < 1)
Y.mtb_plm_M00_rmv <- dplyr::filter(Y.mtb_plm_M00_sum, FeatureTotal < 1)

Y.mtb_prx_F01_rmv$FeatureID <- as.character(Y.mtb_prx_F01_rmv$FeatureID)
Y.mtb_prx_M00_rmv$FeatureID <- as.character(Y.mtb_prx_M00_rmv$FeatureID)
Y.mtb_dtl_F01_rmv$FeatureID <- as.character(Y.mtb_dtl_F01_rmv$FeatureID)
Y.mtb_dtl_M00_rmv$FeatureID <- as.character(Y.mtb_dtl_M00_rmv$FeatureID)
Y.mtb_plm_F01_rmv$FeatureID <- as.character(Y.mtb_plm_F01_rmv$FeatureID)
Y.mtb_plm_M00_rmv$FeatureID <- as.character(Y.mtb_plm_M00_rmv$FeatureID)

Y.abs_mtb_prx_F01_1 <- Y.abs_mtb_prx_F01_0
Y.abs_mtb_prx_M00_1 <- Y.abs_mtb_prx_M00_0
Y.abs_mtb_dtl_F01_1 <- Y.abs_mtb_dtl_F01_0
Y.abs_mtb_dtl_M00_1 <- Y.abs_mtb_dtl_M00_0
Y.abs_mtb_plm_F01_1 <- Y.abs_mtb_plm_F01_0
Y.abs_mtb_plm_M00_1 <- Y.abs_mtb_plm_M00_0
Y.abs_mtb_prx_F01_1$FeatureID <- row.names(Y.abs_mtb_prx_F01_1)
Y.abs_mtb_prx_M00_1$FeatureID <- row.names(Y.abs_mtb_prx_M00_1)
Y.abs_mtb_dtl_F01_1$FeatureID <- row.names(Y.abs_mtb_dtl_F01_1)
Y.abs_mtb_dtl_M00_1$FeatureID <- row.names(Y.abs_mtb_dtl_M00_1)
Y.abs_mtb_plm_F01_1$FeatureID <- row.names(Y.abs_mtb_plm_F01_1)
Y.abs_mtb_plm_M00_1$FeatureID <- row.names(Y.abs_mtb_plm_M00_1)

Y.abs_mtb_prx_F01_2 <- dplyr::anti_join(x = Y.abs_mtb_prx_F01_1, 
                                        y = Y.mtb_prx_F01_rmv, 
                                        by = "FeatureID")
Y.abs_mtb_prx_M00_2 <- dplyr::anti_join(x = Y.abs_mtb_prx_M00_1, 
                                        y = Y.mtb_prx_M00_rmv, 
                                        by = "FeatureID")
Y.abs_mtb_dtl_F01_2 <- dplyr::anti_join(x = Y.abs_mtb_dtl_F01_1, 
                                        y = Y.mtb_dtl_F01_rmv, 
                                        by = "FeatureID")
Y.abs_mtb_dtl_M00_2 <- dplyr::anti_join(x = Y.abs_mtb_dtl_M00_1, 
                                        y = Y.mtb_dtl_M00_rmv, 
                                        by = "FeatureID")
Y.abs_mtb_plm_F01_2 <- dplyr::anti_join(x = Y.abs_mtb_plm_F01_1, 
                                        y = Y.mtb_plm_F01_rmv, 
                                        by = "FeatureID")
Y.abs_mtb_plm_M00_2 <- dplyr::anti_join(x = Y.abs_mtb_plm_M00_1, 
                                        y = Y.mtb_plm_M00_rmv, 
                                        by = "FeatureID")

Y.abs_mtb_prx_F01_3 <- Y.abs_mtb_prx_F01_2
Y.abs_mtb_prx_M00_3 <- Y.abs_mtb_prx_M00_2
Y.abs_mtb_dtl_F01_3 <- Y.abs_mtb_dtl_F01_2
Y.abs_mtb_dtl_M00_3 <- Y.abs_mtb_dtl_M00_2
Y.abs_mtb_plm_F01_3 <- Y.abs_mtb_plm_F01_2
Y.abs_mtb_plm_M00_3 <- Y.abs_mtb_plm_M00_2
row.names(Y.abs_mtb_prx_F01_3) <- Y.abs_mtb_prx_F01_3$FeatureID
row.names(Y.abs_mtb_prx_M00_3) <- Y.abs_mtb_prx_M00_3$FeatureID
row.names(Y.abs_mtb_dtl_F01_3) <- Y.abs_mtb_dtl_F01_3$FeatureID
row.names(Y.abs_mtb_dtl_M00_3) <- Y.abs_mtb_dtl_M00_3$FeatureID
row.names(Y.abs_mtb_plm_F01_3) <- Y.abs_mtb_plm_F01_3$FeatureID
row.names(Y.abs_mtb_plm_M00_3) <- Y.abs_mtb_plm_M00_3$FeatureID
Y.abs_mtb_prx_F01_4 <- dplyr::select(Y.abs_mtb_prx_F01_3, -FeatureID)
Y.abs_mtb_prx_M00_4 <- dplyr::select(Y.abs_mtb_prx_M00_3, -FeatureID)
Y.abs_mtb_dtl_F01_4 <- dplyr::select(Y.abs_mtb_dtl_F01_3, -FeatureID)
Y.abs_mtb_dtl_M00_4 <- dplyr::select(Y.abs_mtb_dtl_M00_3, -FeatureID)
Y.abs_mtb_plm_F01_4 <- dplyr::select(Y.abs_mtb_plm_F01_3, -FeatureID)
Y.abs_mtb_plm_M00_4 <- dplyr::select(Y.abs_mtb_plm_M00_3, -FeatureID)

Y.abs_mtb_prx_F01_5 <- as.data.frame(t(Y.abs_mtb_prx_F01_4))
Y.abs_mtb_prx_M00_5 <- as.data.frame(t(Y.abs_mtb_prx_M00_4))
Y.abs_mtb_plm_F01_5 <- as.data.frame(t(Y.abs_mtb_plm_F01_4))
Y.abs_mtb_plm_M00_5 <- as.data.frame(t(Y.abs_mtb_plm_M00_4))
Y.abs_mtb_dtl_F01_5 <- as.data.frame(t(Y.abs_mtb_dtl_F01_4))
Y.abs_mtb_dtl_M00_5 <- as.data.frame(t(Y.abs_mtb_dtl_M00_4))

### ************************************
### Y - STEP  2 - compute CZM/clr/PCA/ ----
### ************************************

# for CZM and clr, the process occurs as follows:
# (1) replace zero counts
# (2) check proportion conversion
# (3) transpose to convert features to rows and samples to cols
# (4) clr transform
# NOTE: for cmultRepl(), input data = features as cols and samples as rows
# NOTE: for clr transform, input data = features as rows and samples as cols

Y.czm_mtb_prx_F01_0 <- zCompositions::cmultRepl(Y.abs_mtb_prx_F01_5, 
                                                method = "CZM", output = "prop")
Y.czm_mtb_prx_M00_0 <- zCompositions::cmultRepl(Y.abs_mtb_prx_M00_5, 
                                                method = "CZM", output = "prop")
Y.czm_mtb_dtl_F01_0 <- zCompositions::cmultRepl(Y.abs_mtb_dtl_F01_5, 
                                                method = "CZM", output = "prop")
Y.czm_mtb_dtl_M00_0 <- zCompositions::cmultRepl(Y.abs_mtb_dtl_M00_5, 
                                                method = "CZM", output = "prop")
Y.czm_mtb_plm_F01_0 <- zCompositions::cmultRepl(Y.abs_mtb_plm_F01_5, 
                                                method = "CZM", output = "prop")
Y.czm_mtb_plm_M00_0 <- zCompositions::cmultRepl(Y.abs_mtb_plm_M00_5, 
                                                method = "CZM", output = "prop")

print(rowSums(Y.czm_mtb_prx_F01_0)) # all == 1
print(rowSums(Y.czm_mtb_prx_M00_0)) # all == 1
print(rowSums(Y.czm_mtb_dtl_F01_0)) # all == 1
print(rowSums(Y.czm_mtb_dtl_M00_0)) # all == 1
print(rowSums(Y.czm_mtb_plm_F01_0)) # all == 1
print(rowSums(Y.czm_mtb_plm_M00_0)) # all == 1

Y.czm_mtb_prx_F01_1 <- as.data.frame(t(Y.czm_mtb_prx_F01_0))
Y.czm_mtb_prx_M00_1 <- as.data.frame(t(Y.czm_mtb_prx_M00_0))
Y.czm_mtb_dtl_F01_1 <- as.data.frame(t(Y.czm_mtb_dtl_F01_0))
Y.czm_mtb_dtl_M00_1 <- as.data.frame(t(Y.czm_mtb_dtl_M00_0))
Y.czm_mtb_plm_F01_1 <- as.data.frame(t(Y.czm_mtb_plm_F01_0))
Y.czm_mtb_plm_M00_1 <- as.data.frame(t(Y.czm_mtb_plm_M00_0))

Y.clr_mtb_prx_F01_0 <- as.data.frame(
  apply(Y.czm_mtb_prx_F01_1, 2, function(x) {log2(x) - mean(log2(x))})) 
Y.clr_mtb_prx_M00_0 <- as.data.frame(
  apply(Y.czm_mtb_prx_M00_1, 2, function(x) {log2(x) - mean(log2(x))})) 
Y.clr_mtb_dtl_F01_0 <- as.data.frame(
  apply(Y.czm_mtb_dtl_F01_1, 2, function(x) {log2(x) - mean(log2(x))})) 
Y.clr_mtb_dtl_M00_0 <- as.data.frame(
  apply(Y.czm_mtb_dtl_M00_1, 2, function(x) {log2(x) - mean(log2(x))}))
Y.clr_mtb_plm_F01_0 <- as.data.frame(
  apply(Y.czm_mtb_plm_F01_1, 2, function(x) {log2(x) - mean(log2(x))})) 
Y.clr_mtb_plm_M00_0 <- as.data.frame(
  apply(Y.czm_mtb_plm_M00_1, 2, function(x) {log2(x) - mean(log2(x))})) 

# for PCA, the process occurs as follows:
# (1) transpose to convert features to cols and samples to rows
# (2) compute PCA
# (3) create PCA plot axis labels (i.e. isolate variation explained by PC1;PC2)
# (4;5;6) create ggbiplots, extract x & y coordinates, merge with sample data
# NOTE: for prcomp(), input data = features as cols and samples as rows

Y.clr_mtb_prx_F01_1 <- as.data.frame(t(Y.clr_mtb_prx_F01_0))
Y.clr_mtb_prx_M00_1 <- as.data.frame(t(Y.clr_mtb_prx_M00_0))
Y.clr_mtb_dtl_F01_1 <- as.data.frame(t(Y.clr_mtb_dtl_F01_0))
Y.clr_mtb_dtl_M00_1 <- as.data.frame(t(Y.clr_mtb_dtl_M00_0))
Y.clr_mtb_plm_F01_1 <- as.data.frame(t(Y.clr_mtb_plm_F01_0))
Y.clr_mtb_plm_M00_1 <- as.data.frame(t(Y.clr_mtb_plm_M00_0))

Y.pca_mtb_prx_F01 <- prcomp(Y.clr_mtb_prx_F01_1)
Y.pca_mtb_prx_M00 <- prcomp(Y.clr_mtb_prx_M00_1)
Y.pca_mtb_dtl_F01 <- prcomp(Y.clr_mtb_dtl_F01_1)
Y.pca_mtb_dtl_M00 <- prcomp(Y.clr_mtb_dtl_M00_1)
Y.pca_mtb_plm_F01 <- prcomp(Y.clr_mtb_plm_F01_1)
Y.pca_mtb_plm_M00 <- prcomp(Y.clr_mtb_plm_M00_1)

Y.pc1_mtb_prx_F01 <- iso_var_exp(data = Y.pca_mtb_prx_F01, pc.num = 1, 
                                 output= "label")
Y.pc2_mtb_prx_F01 <- iso_var_exp(data = Y.pca_mtb_prx_F01, pc.num = 2, 
                                 output = "label")
Y.pc1_mtb_prx_M00 <- iso_var_exp(data = Y.pca_mtb_prx_M00, pc.num = 1, 
                                 output = "label")
Y.pc2_mtb_prx_M00 <- iso_var_exp(data = Y.pca_mtb_prx_M00, pc.num = 2, 
                                 output = "label")
Y.pc1_mtb_dtl_F01 <- iso_var_exp(data = Y.pca_mtb_dtl_F01, pc.num = 1, 
                                 output = "label")
Y.pc2_mtb_dtl_F01 <- iso_var_exp(data = Y.pca_mtb_dtl_F01, pc.num = 2, 
                                 output = "label")
Y.pc1_mtb_dtl_M00 <- iso_var_exp(data = Y.pca_mtb_dtl_M00, pc.num = 1, 
                                 output = "label")
Y.pc2_mtb_dtl_M00 <- iso_var_exp(data = Y.pca_mtb_dtl_M00, pc.num = 2, 
                                 output = "label")
Y.pc1_mtb_plm_F01 <- iso_var_exp(data = Y.pca_mtb_plm_F01, pc.num = 1, 
                                 output = "label")
Y.pc2_mtb_plm_F01 <- iso_var_exp(data = Y.pca_mtb_plm_F01, pc.num = 2, 
                                 output = "label")
Y.pc1_mtb_plm_M00 <- iso_var_exp(data = Y.pca_mtb_plm_M00, pc.num = 1,
                                 output = "label")
Y.pc2_mtb_plm_M00 <- iso_var_exp(data = Y.pca_mtb_plm_M00, pc.num = 2, 
                                 output = "label")

Y.gbp_mtb_prx_F01 <- ggbiplot::ggbiplot(Y.pca_mtb_prx_F01, scale = 0, 
                                        var.axes = F, 
                                        labels = row.names(Y.clr_mtb_prx_F01_1))
Y.gbp_mtb_prx_M00 <- ggbiplot::ggbiplot(Y.pca_mtb_prx_M00, scale = 0, 
                                        var.axes = F,
                                        labels = row.names(Y.clr_mtb_prx_M00_1))
Y.gbp_mtb_dtl_F01 <- ggbiplot::ggbiplot(Y.pca_mtb_dtl_F01, scale = 0, 
                                        var.axes = F,
                                        labels = row.names(Y.clr_mtb_dtl_F01_1))
Y.gbp_mtb_dtl_M00 <- ggbiplot::ggbiplot(Y.pca_mtb_dtl_M00 , scale = 0, 
                                        var.axes = F,
                                        labels = row.names(Y.clr_mtb_dtl_M00_1))
Y.gbp_mtb_plm_F01 <- ggbiplot::ggbiplot(Y.pca_mtb_plm_F01, scale = 0, 
                                        var.axes = F,
                                        labels = row.names(Y.clr_mtb_plm_F01_1))
Y.gbp_mtb_plm_M00 <- ggbiplot::ggbiplot(Y.pca_mtb_plm_M00 , scale = 0, 
                                        var.axes = F,
                                        labels = row.names(Y.clr_mtb_plm_M00_1))

Y.cord_pca_mtb_prx_F01 <- Y.gbp_mtb_prx_F01[["data"]]
Y.cord_pca_mtb_prx_M00 <- Y.gbp_mtb_prx_M00[["data"]]
Y.cord_pca_mtb_dtl_F01 <- Y.gbp_mtb_dtl_F01[["data"]]
Y.cord_pca_mtb_dtl_M00 <- Y.gbp_mtb_dtl_M00[["data"]]
Y.cord_pca_mtb_plm_F01 <- Y.gbp_mtb_plm_F01[["data"]]
Y.cord_pca_mtb_plm_M00 <- Y.gbp_mtb_plm_M00[["data"]]

Y.smp_cord_pca_mtb_prx_F01 <- merge(x = Y.cord_pca_mtb_prx_F01, y = smp_dat, 
                                    all = F, sort = F, by.x = "labels", 
                                    by.y = "MTBSampleID")
Y.smp_cord_pca_mtb_prx_M00 <- merge(x = Y.cord_pca_mtb_prx_M00, y = smp_dat, 
                                    all = F, sort = F, by.x = "labels", 
                                    by.y = "MTBSampleID")
Y.smp_cord_pca_mtb_dtl_F01 <- merge(x = Y.cord_pca_mtb_dtl_F01, y = smp_dat, 
                                    all = F, sort = F, by.x = "labels", 
                                    by.y = "MTBSampleID")
Y.smp_cord_pca_mtb_dtl_M00 <- merge(x = Y.cord_pca_mtb_dtl_M00, y = smp_dat, 
                                    all = F, sort = F, by.x = "labels", 
                                    by.y = "MTBSampleID")
Y.smp_cord_pca_mtb_plm_F01 <- merge(x = Y.cord_pca_mtb_plm_F01, y = smp_dat, 
                                    all = F, sort = F, by.x = "labels", 
                                    by.y = "MTBSampleID")
Y.smp_cord_pca_mtb_plm_M00 <- merge(x = Y.cord_pca_mtb_plm_M00, y = smp_dat, 
                                    all = F, sort = F, by.x = "labels", 
                                    by.y = "MTBSampleID")

### ************************************
### Y - STEP  3 - define plot parameters/customize plot aesthetics ----
### ***********************************

# NOTE: for consortium group aesthetics, order = Control - Rice.bran.modified

# vectors for proximal colon samples (shape outline and shape fill)
Y.hex_prx_F01 <- c(hex_F01_cmc[3], hex_F01_rmc[3])
Y.fil_prx_F01 <- c(hex_F01_cmc[5], hex_F01_rmc[5])
Y.hex_prx_M00 <- c(hex_M00_cmc[3], hex_M00_rmc[3])
Y.fil_prx_M00 <- c(hex_M00_cmc[5], hex_M00_rmc[5])

# vectors for distal colon samples (shape outline and shape fill)
Y.hex_dtl_F01 <- c(hex_F01_cmc[1], hex_F01_rmc[1])
Y.fil_dtl_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2])
Y.hex_dtl_M00 <- c(hex_M00_cmc[1], hex_M00_rmc[1])
Y.fil_dtl_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2])

# vectors for plasma samples (shape outline and shape fill)
Y.hex_plm_F01 <- c(hex_F01_cmc[3], hex_F01_rmc[3])
Y.fil_plm_F01 <- c(hex_F01_cmc[4], hex_F01_rmc[4])
Y.hex_plm_M00 <- c(hex_M00_cmc[3], hex_M00_rmc[3])
Y.fil_plm_M00 <- c(hex_M00_cmc[4], hex_M00_rmc[4])

# vector for shapes & colors by murine sex; plot order = Female mice - Male mice
Y.sex_shp <- c(shp_mur_fem, shp_mur_mal) 
Y.sex_hex_F01 <- c(hex_F01_cmc[2], hex_F01_rmc[2])
Y.sex_hex_M00 <- c(hex_M00_cmc[2], hex_M00_rmc[2])

# plot labeling
Y.ttl_prx <- "proximal colon"
Y.ttl_dtl <- "distal colon"
Y.ttl_plm <- "plasma"

# sizing parameters
Y.sze_smp_pts <- 1.01 # (pca) size for sample points
Y.sze_smp_txt <- 2.22 # (pca) size for sample text
Y.sze_int_lne <- 0.31 # (pca) size for origin intercept lines
Y.sze_int_crc <- 1.81 # (pca) size for origin circle
Y.sze_int_sqr <- 0.42 # (pca) size for origin square
Y.sze_smp_pts <- 1.01 # (pca) size for sample points
Y.sze_ptl <- 07 # size for plot titles
Y.sze_atl <- 07 # size for axis titles (x and y)
Y.sze_ltx <- 5.5 # size for legend text
Y.sze_lgn_pts <- 1.31 # size for legend points
Y.sze_lgn_stk <- 0.55 # size for legend point stroke

# PCoA and PCA custom theme parameters
Y.JacksonP <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(size = Y.sze_atl),
  plot.title = element_text(hjust = 0.5, size = Y.sze_ptl),
  legend.position = "none",
  panel.ontop = T,
  panel.background = element_rect(fill = NA))

Y.PurpleRain <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  legend.position = c(0.5, 0.5),
  legend.direction = "vertical",
  legend.background = element_rect(fill = greydient[8]),
  legend.key = element_rect(fill = NA, color = NA),
  legend.box.margin = margin(1, 1, 2, 1, "mm"),
  legend.key.size = unit(0, "mm"),
  legend.text = element_text(
    margin = margin(0.5, 0, 0.5, 0, "mm"), size = Y.sze_ltx, 
    color = greydient[1], hjust = 0))

### ************************************
### Y - STEP 4a - plotting: PCA ----
### ************************************

# pca with origin intercept line segments
Y.gpr_pca_mtb_prx_F01_seg <- ggscatter(data = Y.smp_cord_pca_mtb_prx_F01, 
                                       x = "xvar", y = "yvar",
                                       font.family = fnt, shape = 32, 
                                       color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = Y.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = Y.sze_int_sqr)

Y.gpr_pca_mtb_prx_M00_seg <- ggscatter(data = Y.smp_cord_pca_mtb_prx_M00, 
                                       x = "xvar", y = "yvar",
                                       font.family = fnt, shape = 32, 
                                       color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = Y.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = Y.sze_int_sqr)

Y.gpr_pca_mtb_dtl_F01_seg <- ggscatter(data = Y.smp_cord_pca_mtb_dtl_F01, 
                                       x = "xvar", y = "yvar",
                                       font.family = fnt, shape = 32, 
                                       color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = Y.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = Y.sze_int_sqr)

Y.gpr_pca_mtb_dtl_M00_seg <- ggscatter(data = Y.smp_cord_pca_mtb_dtl_M00, 
                                       x = "xvar", y = "yvar",
                                       font.family = fnt, shape = 32, 
                                       color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = Y.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = Y.sze_int_sqr)

Y.gpr_pca_mtb_plm_F01_seg <- ggscatter(data = Y.smp_cord_pca_mtb_plm_F01, 
                                       x = "xvar", y = "yvar",
                                       font.family = fnt, shape = 32, 
                                       color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = Y.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = Y.sze_int_sqr)

Y.gpr_pca_mtb_plm_M00_seg <- ggscatter(data = Y.smp_cord_pca_mtb_plm_M00, 
                                       x = "xvar", y = "yvar",
                                       font.family = fnt, shape = 32, 
                                       color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7], 
             size = Y.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8], 
             size = Y.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7], 
             size = Y.sze_int_sqr)

# pca for proximal colon
Y.gpr_pca_mtb_prx_F01 <- ggscatter(data = Y.smp_cord_pca_mtb_prx_F01,
                                   title = Y.ttl_prx, shape = shp_prx,
                                   xlab = Y.pc1_mtb_prx_F01, 
                                   ylab = Y.pc2_mtb_prx_F01,
                                   ggp = Y.gpr_pca_mtb_prx_F01_seg,
                                   x = "xvar", y = "yvar", 
                                   color = "ConsortiumAbrv", 
                                   fill = "ConsortiumAbrv", 
                                   font.family = fnt, size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.hex_prx_F01) +
  scale_fill_manual(values = Y.fil_prx_F01) +
  border(color = greydient[1]) +
  Y.JacksonP

Y.gpr_pca_mtb_prx_M00 <- ggscatter(data = Y.smp_cord_pca_mtb_prx_M00,
                                   title = Y.ttl_prx, shape = shp_prx,
                                   xlab = Y.pc1_mtb_prx_M00, 
                                   ylab = Y.pc2_mtb_prx_M00,
                                   ggp = Y.gpr_pca_mtb_prx_M00_seg,
                                   x = "xvar", y = "yvar",
                                   color = "ConsortiumAbrv", 
                                   fill = "ConsortiumAbrv", 
                                   font.family = fnt, size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.hex_prx_M00) +
  scale_fill_manual(values = Y.fil_prx_M00) +
  border(color = greydient[1]) +
  Y.JacksonP

# pca for distal colon
Y.gpr_pca_mtb_dtl_F01 <- ggscatter(data = Y.smp_cord_pca_mtb_dtl_F01,
                                   title = Y.ttl_dtl, shape = shp_dtl,
                                   xlab = Y.pc1_mtb_dtl_F01, 
                                   ylab = Y.pc2_mtb_dtl_F01,
                                   ggp = Y.gpr_pca_mtb_dtl_F01_seg,
                                   x = "xvar", y = "yvar",
                                   color = "ConsortiumAbrv", 
                                   fill = "ConsortiumAbrv", 
                                   font.family = fnt, size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.hex_dtl_F01) +
  scale_fill_manual(values = Y.fil_dtl_F01) +
  border(color = greydient[1]) +
  Y.JacksonP

Y.gpr_pca_mtb_dtl_M00 <- ggscatter(data = Y.smp_cord_pca_mtb_dtl_M00,
                                   title = Y.ttl_dtl, shape = shp_dtl,
                                   xlab = Y.pc1_mtb_dtl_M00, 
                                   ylab = Y.pc2_mtb_dtl_M00,
                                   ggp = Y.gpr_pca_mtb_dtl_M00_seg,
                                   x = "xvar", y = "yvar",
                                   color = "ConsortiumAbrv", 
                                   fill = "ConsortiumAbrv", 
                                   font.family = fnt, size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.hex_dtl_M00) +
  scale_fill_manual(values = Y.fil_dtl_M00) +
  border(color = greydient[1]) +
  Y.JacksonP

# pca for plasma
Y.gpr_pca_mtb_plm_F01 <- ggscatter(data = Y.smp_cord_pca_mtb_plm_F01,
                                   title = Y.ttl_plm, shape = shp_plm,
                                   xlab = Y.pc1_mtb_plm_F01, 
                                   ylab = Y.pc2_mtb_plm_F01,
                                   ggp = Y.gpr_pca_mtb_plm_F01_seg,
                                   x = "xvar", y = "yvar",
                                   color = "ConsortiumAbrv", 
                                   fill = "ConsortiumAbrv", 
                                   font.family = fnt, size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.hex_plm_F01) +
  scale_fill_manual(values = Y.fil_plm_F01) +
  border(color = greydient[1]) +
  Y.JacksonP

Y.gpr_pca_mtb_plm_M00 <- ggscatter(data = Y.smp_cord_pca_mtb_plm_M00,
                                   title = Y.ttl_plm, shape = shp_plm,
                                   xlab = Y.pc1_mtb_plm_M00, 
                                   ylab = Y.pc2_mtb_plm_M00,
                                   ggp = Y.gpr_pca_mtb_plm_M00_seg,
                                   x = "xvar", y = "yvar",
                                   color = "ConsortiumAbrv", 
                                   fill = "ConsortiumAbrv", 
                                   font.family = fnt, size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.hex_plm_M00) +
  scale_fill_manual(values = Y.fil_plm_M00) +
  border(color = greydient[1]) +
  Y.JacksonP

# pca for proximal colon with samples shaped by murine sex
Y.gpr_pca_mtb_prx_F01_sex <- ggscatter(data = Y.smp_cord_pca_mtb_prx_F01,
                                       title = Y.ttl_prx, 
                                       xlab = Y.pc1_mtb_prx_F01, 
                                       ylab = Y.pc2_mtb_prx_F01,
                                       ggp = Y.gpr_pca_mtb_prx_F01_seg,
                                       x = "xvar", y = "yvar",
                                       shape = "AnimalSex", 
                                       color = "ConsortiumAbrv",
                                       fill = greydient[8], font.family = fnt, 
                                       size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.sex_hex_F01) +
  scale_shape_manual(values = Y.sex_shp) +
  border(color = greydient[1]) +
  Y.JacksonP

Y.gpr_pca_mtb_prx_M00_sex <- ggscatter(data = Y.smp_cord_pca_mtb_prx_M00,
                                       title = Y.ttl_prx, 
                                       xlab = Y.pc1_mtb_prx_M00, 
                                       ylab = Y.pc2_mtb_prx_M00,
                                       ggp = Y.gpr_pca_mtb_prx_M00_seg,
                                       x = "xvar", y = "yvar",
                                       shape = "AnimalSex", 
                                       color = "ConsortiumAbrv",
                                       fill = greydient[8], font.family = fnt, 
                                       size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.sex_hex_M00) +
  scale_shape_manual(values = Y.sex_shp) +
  border(color = greydient[1]) +
  Y.JacksonP

# pca for distal colon with samples shaped by murine sex
Y.gpr_pca_mtb_dtl_F01_sex <- ggscatter(data = Y.smp_cord_pca_mtb_dtl_F01,
                                       title = Y.ttl_dtl, 
                                       xlab = Y.pc1_mtb_dtl_F01, 
                                       ylab = Y.pc2_mtb_dtl_F01,
                                       ggp = Y.gpr_pca_mtb_dtl_F01_seg,
                                       x = "xvar", y = "yvar",
                                       shape = "AnimalSex", 
                                       color = "ConsortiumAbrv",
                                       fill = greydient[8], font.family = fnt, 
                                       size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.sex_hex_F01) +
  scale_shape_manual(values = Y.sex_shp) +
  border(color = greydient[1]) +
  Y.JacksonP

Y.gpr_pca_mtb_dtl_M00_sex <- ggscatter(data = Y.smp_cord_pca_mtb_dtl_M00,
                                       title = Y.ttl_dtl, 
                                       xlab = Y.pc1_mtb_dtl_M00, 
                                       ylab = Y.pc2_mtb_dtl_M00,
                                       ggp = Y.gpr_pca_mtb_dtl_M00_seg,
                                       x = "xvar", y = "yvar",
                                       shape = "AnimalSex", 
                                       color = "ConsortiumAbrv",
                                       fill = greydient[8], font.family = fnt, 
                                       size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.sex_hex_M00) +
  scale_shape_manual(values = Y.sex_shp) +
  border(color = greydient[1]) +
  Y.JacksonP

# pca for plasma with samples shaped by murine sex
Y.gpr_pca_mtb_plm_F01_sex <- ggscatter(data = Y.smp_cord_pca_mtb_plm_F01,
                                       title = Y.ttl_plm,
                                       xlab = Y.pc1_mtb_plm_F01, 
                                       ylab = Y.pc2_mtb_plm_F01,
                                       ggp = Y.gpr_pca_mtb_plm_F01_seg,
                                       x = "xvar", y = "yvar",
                                       shape = "AnimalSex", 
                                       color = "ConsortiumAbrv",
                                       fill = greydient[8], font.family = fnt, 
                                       size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.sex_hex_F01) +
  scale_shape_manual(values = Y.sex_shp) +
  border(color = greydient[1]) +
  Y.JacksonP

Y.gpr_pca_mtb_plm_M00_sex <- ggscatter(data = Y.smp_cord_pca_mtb_plm_M00,
                                       title = Y.ttl_plm,
                                       xlab = Y.pc1_mtb_plm_M00, 
                                       ylab = Y.pc2_mtb_plm_M00,
                                       ggp = Y.gpr_pca_mtb_plm_M00_seg,
                                       x = "xvar", y = "yvar",
                                       shape = "AnimalSex", 
                                       color = "ConsortiumAbrv",
                                       fill = greydient[8], font.family = fnt, 
                                       size = Y.sze_smp_pts) +
  scale_color_manual(values = Y.sex_hex_M00) +
  scale_shape_manual(values = Y.sex_shp) +
  border(color = greydient[1]) +
  Y.JacksonP

### ************************************
### Y - STEP 4b - plotting: create custom legends ----
### ************************************

# consortium group-based legends (specific to donor groups)
# create data.frames with arbitrary x and y values that will be plotted ...
# ... and define legend breaks and labels

Y.lgn_grp <- data.frame(group = c("C0", "C1", "C2", 
                                  "e",
                                  "R0", "R1", "R2"), x = 0, y = 0)

Y.lgn_grp <- data.frame(group = c("A", "B", "C", "D", 
                                  "E",
                                  "F", "G", "H", 
                                  "I", 
                                  "J", "K", "L"), x = 0, y = 0,
                        stringsAsFactors = F)

Y.lgn_grp_brk <- c("A", "B", "C", "D", 
                   "E",
                   "F", "G", "H", 
                   "I", 
                   "J", "K", "L")

Y.lgn_grp_lab_F01 <- c("CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-f: proximal colon", 
                       "CMC-f: distal colon",
                       "CMC-f: plasma",
                       " ",
                       "RMC-f: proximal colon", 
                       "RMC-f: distal colon",
                       "RMC-f: plasma")

Y.lgn_grp_lab_M00 <- c("CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-m: proximal colon", 
                       "CMC-m: distal colon",
                       "CMC-m: plasma",
                       " ",
                       "RMC-m: proximal colon", 
                       "RMC-m: distal colon",
                       "RMC-m: plasma")

# vector for shape codes by sample Type; 
# order = P - D - plm for CMC then RMC (defined above)
Y.lgn_shp <- c(32, 32, 32, 32, 32,
               shp_prx, shp_dtl, shp_plm, 
               32, 
               shp_prx, shp_dtl, shp_plm)

# vector for shape colors by consortium group/sample type;
# vector for shape fill by consortium group/sample type;
Y.lgn_hex_F01 <- c(rep(greydient[8], times = 5),
                   hex_F01_cmc[3], hex_F01_cmc[1], hex_F01_cmc[3],
                   greydient[8],
                   hex_F01_rmc[3], hex_F01_rmc[1], hex_F01_rmc[3])
Y.lgn_fil_F01 <- c(rep(greydient[8], times = 5),
                   hex_F01_cmc[5], hex_F01_cmc[2], hex_F01_cmc[4],
                   greydient[8],
                   hex_F01_rmc[5], hex_F01_rmc[2], hex_F01_rmc[4])

Y.lgn_hex_M00 <- c(rep(greydient[8], times = 5),
                   hex_M00_cmc[3], hex_M00_cmc[1], hex_M00_cmc[3], 
                   greydient[8],
                   hex_M00_rmc[3], hex_M00_rmc[1], hex_M00_rmc[3])
Y.lgn_fil_M00 <- c(rep(greydient[8], times = 5),
                   hex_M00_cmc[5], hex_M00_cmc[2], hex_M00_cmc[4], 
                   greydient[8],
                   hex_M00_rmc[5], hex_M00_rmc[2], hex_M00_rmc[5])

# (1) legend for CMC-f and RMC-f:
# (2) legend for CMC-m and RMC-m:

Y.gtb_lgn_F01 <- ggplot(data = Y.lgn_grp, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group, fill = group),
             size = Y.sze_lgn_pts, stroke = Y.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = Y.lgn_shp, 
                     breaks = Y.lgn_grp_brk, labels = Y.lgn_grp_lab_F01) +
  scale_color_manual(name = NULL, values = Y.lgn_hex_F01, 
                     breaks = Y.lgn_grp_brk, labels = Y.lgn_grp_lab_F01) +
  scale_fill_manual(name = NULL, values = Y.lgn_fil_F01, 
                    breaks = Y.lgn_grp_brk, labels = Y.lgn_grp_lab_F01) +
  theme_void(base_family = fnt) +
  Y.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

Y.gtb_lgn_M00 <- ggplot(data = Y.lgn_grp, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group, fill = group),
             size = Y.sze_lgn_pts, stroke = Y.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = Y.lgn_shp, 
                     breaks = Y.lgn_grp_brk, labels = Y.lgn_grp_lab_M00) +
  scale_color_manual(name = NULL, values = Y.lgn_hex_M00, 
                     breaks = Y.lgn_grp_brk, labels = Y.lgn_grp_lab_M00) +
  scale_fill_manual(name = NULL, values = Y.lgn_fil_M00, 
                    breaks = Y.lgn_grp_brk, labels = Y.lgn_grp_lab_M00) +
  theme_void(base_family = fnt) +
  Y.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))


# extract legends and convert from class 'gtable' to 'ggplot'
Y.ggp_lgn_F01 <- get_legend(Y.gtb_lgn_F01)
Y.ggp_lgn_M00 <- get_legend(Y.gtb_lgn_M00)
Y.gpr_lgn_F01 <- as_ggplot(Y.ggp_lgn_F01)
Y.gpr_lgn_M00 <- as_ggplot(Y.ggp_lgn_M00)

# murine sex-based legends
# create data.frames with arbitrary x and y values that will be plotted ...
# ... and define legend breaks and labels

Y.lgn_sex_F01 <- data.frame(group = c("A", "B", "C", "D", 
                                      "E",
                                      "F",
                                      "G", 
                                      "H", "I"), x = 0, y = 0,
                            stringsAsFactors = F)
Y.lgn_sex_brk_F01 <- c("A", "B", "C", "D", 
                       "E",
                       "F",
                       "G", 
                       "H", "I")
Y.lgn_sex_lab_F01 <- c("CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-f: Female mice",
                       " ",
                       "RMC-f: Female mice", 
                       "RMC-f: Male mice")

Y.lgn_sex_M00 <- data.frame(group = c("A", "B", "C", "D", 
                                      "E",
                                      "F", "G", 
                                      "H", 
                                      "I", "J"), x = 0, y = 0,
                            stringsAsFactors = F)
Y.lgn_sex_brk_M00 <- c("A", "B", "C", "D", 
                       "E",
                       "F", "G", 
                       "H", 
                       "I", "J")
Y.lgn_sex_lab_M00 <- c("CMC: Control Microbial",
                       "     Consortium",
                       "RMC: Rice bran modified",
                       "     Microbial Consortium",
                       "",
                       "CMC-f: Female mice", 
                       "CMC-f: Male mice", 
                       " ",
                       "RMC-f: Female mice",
                       "RMC-f: Male mice")

# vector for colors by consortium group; vectors for shapes by murine sex
# order for F01 = C.Female - e - R.s - R.Female - R.male
# order for M00 = C.s - C.Female - C.Male - e - R.s - R.Female - R.male
Y.sex_hex_F01 <- c(rep(greydient[8], times = 5),
                   hex_F01_cmc[2],
                   greydient[8],
                   hex_F01_rmc[2], hex_F01_rmc[2])
Y.sex_hex_M00 <- c(rep(greydient[8], times = 5),
                   hex_M00_cmc[2], hex_M00_cmc[2],
                   greydient[8],
                   hex_M00_rmc[2], hex_M00_rmc[2])
Y.sex_shp_F01 <- c(32, 32, 32, 32, 
                   32, 
                   shp_mur_fem, 
                   32,
                   shp_mur_fem, shp_mur_mal)
Y.sex_shp_M00 <- c(32, 32, 32, 32, 
                   32, 
                   shp_mur_fem, shp_mur_mal, 
                   32, 
                   shp_mur_fem, shp_mur_mal)

# (1) legend for CMC-f and RMC-f:
# (2) legend for CMC-m and RMC-m:

Y.gtb_lgn_sex_F01 <- ggplot(data = Y.lgn_sex_F01, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group), fill = greydient[8],
             size = Y.sze_lgn_pts, stroke = Y.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = Y.sex_shp_F01,
                     breaks = Y.lgn_sex_brk_F01, labels = Y.lgn_sex_lab_F01) +
  scale_color_manual(name = NULL, values = Y.sex_hex_F01,
                     breaks = Y.lgn_sex_brk_F01, labels = Y.lgn_sex_lab_F01) +
  theme_void(base_family = fnt) +
  Y.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

Y.gtb_lgn_sex_M00 <- ggplot(data = Y.lgn_sex_M00, aes(x = x, y = y)) +
  geom_point(aes(shape = group, color = group), fill = greydient[8],
             size = Y.sze_lgn_pts, stroke = Y.sze_lgn_stk) +
  scale_shape_manual(name = NULL, values = Y.sex_shp_M00,
                     breaks = Y.lgn_sex_brk_M00, labels = Y.lgn_sex_lab_M00) +
  scale_color_manual(name = NULL, values = Y.sex_hex_M00,
                     breaks = Y.lgn_sex_brk_M00, labels = Y.lgn_sex_lab_M00) +
  theme_void(base_family = fnt) +
  Y.PurpleRain +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# extract legends and convert from class 'gtable' to 'ggplot'
Y.ggp_lgn_sex_F01 <- get_legend(Y.gtb_lgn_sex_F01)
Y.ggp_lgn_sex_M00 <- get_legend(Y.gtb_lgn_sex_M00)
Y.gpr_lgn_sex_F01 <- as_ggplot(Y.ggp_lgn_sex_F01)
Y.gpr_lgn_sex_M00 <- as_ggplot(Y.ggp_lgn_sex_M00)

### ************************************
### Y - STEP 4c - plotting: arrange ----
### ************************************

# arrange plots by sex of human donor
Y.gga_mtb_pca_F01 <- ggarrange(Y.gpr_pca_mtb_prx_F01, Y.gpr_pca_mtb_dtl_F01, 
                               Y.gpr_pca_mtb_plm_F01,
                               labels = c("A", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = c(1, 1, 1), 
                               align = "hv")

Y.gga_mtb_pca_M00 <- ggarrange(Y.gpr_pca_mtb_prx_M00, Y.gpr_pca_mtb_dtl_M00, 
                               Y.gpr_pca_mtb_plm_M00,
                               labels = c("B", "", ""), font.label = pan_fnt,
                               ncol = 3, nrow = 1, widths = c(1, 1, 1), 
                               align = "hv")

# arrange plots and legends by sex of human donor
Y.gga_plt_lgn_F01 <- ggarrange(Y.gga_mtb_pca_F01, Y.gpr_lgn_F01, labels = NULL,
                               ncol = 2, nrow = 1, widths = c(3, 0.95))
Y.gga_plt_lgn_M00 <- ggarrange(Y.gga_mtb_pca_M00, Y.gpr_lgn_M00, labels = NULL,
                               ncol = 2, nrow = 1, widths = c(3, 0.95))

# arrange arrangements
Y.gga_mtb_pca <- ggarrange(Y.gga_plt_lgn_F01, Y.gga_plt_lgn_M00, labels = NULL,
                           ncol = 1, nrow = 2, align = "hv")


# murine sex-based plots - arrange plots by sex of human donor
Y.gga_mtb_pca_F01_sex <- ggarrange(Y.gpr_pca_mtb_prx_F01_sex, 
                                   Y.gpr_pca_mtb_dtl_F01_sex, 
                                   Y.gpr_pca_mtb_plm_F01_sex,
                                   labels = c("C", "", ""), 
                                   font.label = pan_fnt, ncol = 3, nrow = 1, 
                                   widths = c(1, 1, 1), align = "hv")

Y.gga_mtb_pca_M00_sex <- ggarrange(Y.gpr_pca_mtb_prx_M00_sex, 
                                   Y.gpr_pca_mtb_dtl_M00_sex, 
                                   Y.gpr_pca_mtb_plm_M00_sex,
                                   labels = c("F", "", ""), 
                                   font.label = pan_fnt, ncol = 3, nrow = 1, 
                                   widths = c(1, 1, 1), align = "hv")

# murine sex-based plots - arrange plots and legends by sex of human donor
Y.gga_plt_sex_lgn_F01 <- ggarrange(Y.gga_mtb_pca_F01_sex, Y.gpr_lgn_sex_F01,
                                   labels = NULL, ncol = 2, nrow = 1,
                                   widths = c(3, 0.95))
Y.gga_plt_sex_lgn_M00 <- ggarrange(Y.gga_mtb_pca_M00_sex, Y.ggp_lgn_sex_M00,
                                   labels = NULL, ncol = 2, nrow = 1, 
                                   widths = c(3, 0.95))

### ************************************
### Y - WRITE OUTPUTS ----
### ************************************

# consortium-based plots
Y.ofv_gga_plt_lgn_F01 <- "TransFaunation/vault/plot_pca_mtb_F01.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 45,
       filename = Y.ofv_gga_plt_lgn_F01, plot = Y.gga_plt_lgn_F01)

Y.ofv_gga_plt_lgn_M00 <- "TransFaunation/vault/plot_pca_mtb_M00.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 45,
       filename = Y.ofv_gga_plt_lgn_M00, plot = Y.gga_plt_lgn_M00)

Y.ofv_gga_mtb_pca <- "TransFaunation/vault/plot_pca_mtb.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 90,
       filename = Y.ofv_gga_mtb_pca, plot = Y.gga_mtb_pca)

# murine sex-based plots
Y.ofv_gga_F01_sex <- "TransFaunation/vault/plot_pca_mtb_F01_animal_sex.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 45,
       filename = Y.ofv_gga_F01_sex, plot = Y.gga_plt_sex_lgn_F01)

Y.ofv_gga_M00_sex <- "TransFaunation/vault/plot_pca_mtb_M00_animal_sex.pdf"
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 45,
       filename = Y.ofv_gga_M00_sex, plot = Y.gga_plt_sex_lgn_M00)

# save workspace
Y.obj <- ls(pattern = "Y.")
Y.lst <- c(Y.obj[grep(pattern = "Y.", x = Y.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS)
save(list = Y.lst, file = Y.ofv_wksp)
