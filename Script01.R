### ************************************
### PREFACE ----
### ************************************

# this file houses the entirety of code for analysis performed in QIIME2^
# ^after creating a manifest in R
# ** note for KDP: PREFACE version 0.5 ** #

# NOTE: my machine runs MacOS; paths have not been tested on other platforms
# use objects created in '# configure the script' when planning for paths
# for QIIME2, the working directory will ALWAYS be: print(path_projd)
# paths will ALWAYS be relative to that working directory
# raw fastq files are located in: path_projd/raw_fastq/
# database files for tax classifiers are located in the home directory: ~/
# outputs are stored in the central vault located at: print(path_vault)

# configure the script
setwd("~/Desktop/") # *gasp*
name_scrpt <- "Script01" # script filename
name_projd <- "TransFaunation" # name of project directory
path_projd <- paste(sep = "", getwd(), "/", name_projd) # path to project dir
path_vault <- paste(sep = "", path_projd, "/vault") # path to storage vault
PREFACE <- c("name_scrpt", "name_projd", "path_projd", "path_vault")

# paths to universal inputs stored in the project directory (path_projd)
ifp_smp_dat <- paste(sep = "", path_projd, "/sampledata_TF.txt")

# paths to section specific inputs stored in the project directory (path_projd)
M.ifp_fsq <- paste(sep = "", path_projd, "/raw_fastq/")

# paths for outputs stored in the central vault (path_vault)
I.ofv_info <- paste(sep = "", path_vault, "/info_", name_scrpt, "_SectionI.txt")
I.ofv_wksp <- paste(sep = "", path_vault, "/WS_", name_scrpt, "_SectionI.RData")

M.ofv_mft <- paste(sep = "", path_vault, "/manifest_R1.csv")
M.ofv_info <- paste(sep = "", path_vault, "/info_", name_scrpt, "_SectionM.txt")
M.ofv_prov <- paste(sep = "", path_vault, "/prov_", name_scrpt,"_SectionM.txt")
M.ofv_wksp <- paste(sep = "", path_vault, "/WS_", name_scrpt, "_SectionM.RData")

Q.ofv_info <- paste(sep = "", path_vault, "/info_", name_scrpt, "_SectionQ.txt")
Q.ofv_prov <- paste(sep = "", path_vault, "/prov_", name_scrpt,"_SectionQ.txt")
Q.ofv_wksp <- paste(sep = "", path_vault, "/WS_", name_scrpt, "_SectionQ.RData")

# save PREFACE workspace
PREFACE.lst <- c(ls(pattern = "ifp"), ls(pattern = "ifv"), ls(pattern = "ofv"), 
                 PREFACE, "PREFACE", "PREFACE.lst")
save(list = PREFACE.lst,
     file = paste(sep = "", path_vault, "/WS_", name_scrpt, "_PREFACE.RData"))

### ************************************
### I - MACHINE/PACKAGE/VERSION Info ----
### ************************************

# NOTE: section I requires objects from the PREFACE to be in the environment
# ** note for KDP: section I version 0.2 ** #

# capture QIIME2-related information:
# NOTE: Q2 version info curated manually from Terminal print out: $ qiime info
# MOTE: Greengenes/SILVA information curated manually
I.qcli_ctg <- "QIIME2"
I.qcli_sec <- "Q"

I.QIIME2t_a <- data.frame(info = "QIIME2 version", value = "2019.7.0",
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_b <- data.frame(info = "Python version", value = "3.6.7",
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_c <- data.frame(info = "Greengenes database version", value = "13_8",
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_d <- data.frame(info = "Greengenes OTU file used", 
                          value = "gg_13_8_otus/rep_set/99_otus.fasta",
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_e <- data.frame(info = "Greengenes taxonomy file used", 
                          value = "gg_13_8_otus/taxonomy/99_otu_taxonomy.txt",
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_f <- data.frame(info = "SILVA database version", value = "132",
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_g <- data.frame(info = "SILVA OTU file used", 
                          value = 
                            paste(sep = "", "SILVA_132_QIIME_release/rep_set/",
                                  "rep_set_16S_only/99/silva_132_99_16S.fna"),
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_h <- data.frame(info = "SILVA taxonomy file used", 
                          value = 
                            paste(sep = "", "SILVA_132_QIIME_release/taxonomy/",
                                  "16S_only/99/taxonomy_7_levels.txt"),
                          category = I.qcli_ctg, script = name_scrpt,
                          section = I.qcli_sec, stringsAsFactors = F)
I.QIIME2t_0 <- rbind(I.QIIME2t_a, I.QIIME2t_b, I.QIIME2t_c, I.QIIME2t_d,
                     I.QIIME2t_e, I.QIIME2t_f, I.QIIME2t_g, I.QIIME2t_h)
I.QIIME2t <- dplyr::select(I.QIIME2t_0, category, info, value, script, section)

# R packages accessed via namespace:
# benchmarkme
# dplyr

# capture R package-related information:
I.Rpac_ctg <- "R package version"
I.Rpackge_a <- data.frame(info = "benchmarkme", section = "I", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_b <- data.frame(info = "dplyr", section = "I", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_0 <- rbind(I.Rpackge_a, I.Rpackge_b)

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

I.pathsto_a <- data.frame(info = "working directory for R (Console)",
                          value = paste(getwd(), "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_b <- data.frame(info = "working directory for QIIME2 (Terminal)",
                          value = path_projd,
                          category = I.path_ctg, script = name_scrpt,
                          section = "Q", stringsAsFactors = F)
I.pathsto_c <- data.frame(info = "path to project directory",
                          value = paste(path_projd, "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_d <- data.frame(info = "path to project's central '/vault/'",
                          value = paste(path_vault, "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_0 <- rbind(I.pathsto_a, I.pathsto_b, I.pathsto_c, I.pathsto_d)
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
I.info <- rbind(I.project, I.machine, I.pathsto, I.baseRpl, I.QIIME2t,
                I.RStudio, I.Rpackge)

# outputs to the vault
write.table(sep = "\t", row.names = F, x = I.info, file = I.ofv_info)

I.obj <- ls(pattern = "I.")
I.lst <- c(I.obj[grep(pattern = "I.", x = I.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, "captureVersions")
save(list = I.lst, file = I.ofv_wksp)

### ************************************
### DESCRIPTIONS FOR EACH SECTION ----
### ************************************

## SECTION M - create manifest for QIIME2
# purpose:
# this section creates a comma-separated "Fastq manifest" file (.csv) that ...
# ... is used to import individual (demultiplexed) fastq files into QIIME2
# NOTE: the fastq files are Casava 1.8 formatted, nullifying the need for ...
# ... a manifest file; however, this section also serves to validate ...
# ... that SampleIDs in the sampledata file match the fastq file names ...
# ... meaning that it ensures the correct samples are imported into QIIME2

# description of steps:
# STEP 1 = creates three columns (sample-id, absolute-filepath, direction)
# STEP 2 = creates a manifest and checks the manifest/sampledata for errors

## SECTION Q - QIIME2 analysis

# purpose:
# this section stores all QIIME2 commands
# unless otherwise noted, commands were exectued using the R Studio Terminal
# if viewing this script in R Studio, it is recommended to uncheck the box here: 
# Tools --> Global Options --> Code --> Diagnostics --> Show diagnostics for R

# NOTE: some R code is present within a few subsections of section Q
# that code stores information obtained from .qzv files ... 
# ... and also provides some level of provenance outside of this script
# i.e. from reads.qzv: 'Demultiplexed sequence counts summary - Total:'
# as a result of needing to run R code, make sure to run the PREFACE section ...
# ... in order to have the appropiate paths and variables in the environment

### ************************************
### M - STEP 1 - create columns ----
### ************************************

# ** note for KDP: MS create manifest version 0.5 ** #

# for column 'sample-id', the process occurs as follows:
# (1) read in the filenames and simultaneously replace the file extension
# (2) successively truncate the extension-less filenames to obtain the SampleID
# NOTE: this process uses the Illumina generated sample number to truncate ...
# .. the extension-less filenames (this also removes the Cassava 1.8 formatting)
# PRACTICAL EXAMPLE:
# 1Swab_S51_L001_R1_001.fastq.gz = raw filename
# 1Swab_S51_L001_R1_001 = extension-less filename; obtained after (1)
# 1Swab = SampleID; obtained after (2)
# S51 = Illumina sample number
# _L001_R1_001 = Cassava 1.8 formatting
# .fastq.gz = file extension
# WARNING: if a sampleID is named with _S (i.e. Control_Swab) then ... 
# this process will not occur correctly

M.fname <- gsub(pattern = ".fastq.gz", replacement = "",
                list.files(path = M.ifp_fsq, full.names = F))

M.sID_0 <- sapply(M.fname, function(x) {strsplit(x, split = "_S0")[[1]][1]})
M.sID_1 <- sapply(M.sID_0, function(x) {strsplit(x, split = "_S1")[[1]][1]})
M.sID_2 <- sapply(M.sID_1, function(x) {strsplit(x, split = "_S2")[[1]][1]})
M.sID_3 <- sapply(M.sID_2, function(x) {strsplit(x, split = "_S3")[[1]][1]})
M.sID_4 <- sapply(M.sID_3, function(x) {strsplit(x, split = "_S4")[[1]][1]})
M.sID_5 <- sapply(M.sID_4, function(x) {strsplit(x, split = "_S5")[[1]][1]})
M.sID_6 <- sapply(M.sID_5, function(x) {strsplit(x, split = "_S6")[[1]][1]})
M.sID_7 <- sapply(M.sID_6, function(x) {strsplit(x, split = "_S7")[[1]][1]})
M.sID_8 <- sapply(M.sID_7, function(x) {strsplit(x, split = "_S8")[[1]][1]})
M.sID_9 <- sapply(M.sID_8, function(x) {strsplit(x, split = "_S9")[[1]][1]})

# for column 'absolute-filepath':
# obtain file names with extensions and prepend names with the filepath
M.afp <- gsub(pattern = M.ifp_fsq, replacement = "$PWD/raw_fastq",
              list.files(path = M.ifp_fsq, full.names = T))

# for column 'direction':
# match pattern in filenames to specify the read direction
M.drc <- sapply(M.afp, function(x) {
  ifelse(grepl(pattern = "R1", x, ignore.case = F), yes = "forward", 
         no = ifelse(grepl(pattern = "R2", x, ignore.case = F),
                     yes = "reverse", no = "index"))})

### ************************************
### M - STEP 2 - create manifest ----
### ************************************

# ** note for KDP: MS create manifest version 0.3 ** #

# create manifest df containing info from the vetors created above
M.mft_0 <- data.frame("sample_id" = M.sID_9,
                      "absolute_filepath" = M.afp,
                      "direction" = M.drc)

# check to confirm manifest SampleIDs match sampledata SampleIDs:
# this process occurs as follows:
# (1;2) remove duplicate rows from manifest and create sorted vector
# (3;4;5) read in sample data, retain metataxonomics, create sorted vector
# (6;7) run logical check

M.mft_unq_0 <- data.frame("SampleID" = unique(M.mft_0$sample_id))
M.mft_unq_1 <- as.character(sort(M.mft_unq_0$SampleID))

M.smp_0 <- read.table(ifp_smp_dat, header = T, sep = "\t", as.is = T)
M.smp_1 <- dplyr::filter(M.smp_0, MTXFilter == "yay")
M.smp_2 <- as.character(sort(M.smp_1$SampleID))

LOGICAL_M.sID <- all(M.mft_unq_1 == M.smp_2)

if (!isTRUE(LOGICAL_M.sID)) {
  stop("SampleIDs in manifest do not match SampleIDs in sampledata file")
}

# format manifest prior to output, this process occurs as follows:
# (1) filter to isolate forward read SampleIDs only
# (2) create a copy to avoid overwriting data
# (3) format column names to match requirements of QIIME2 manifest
# (4) copy to remove variable number

M.mft_1 <- dplyr::filter(M.mft_0, direction == "forward")
M.mft_2 <- M.mft_1
names(M.mft_2) <- gsub(pattern = "\\_", replacement = "-", x = names(M.mft_2))
M.mft <- M.mft_2

### ************************************
### M - STEP 3 - information gathering and provenance ----
### ************************************

# count the number of samples in the manifest
M.info <- data.frame("info" = "number of samples in manifest", 
                     value = length(M.mft$`sample-id`), "object" = "M.mft", 
                     "script" = name_scrpt, "section" = "Section M - STEP 2",
                     "heading" = "create manifest", stringsAsFactors = F)

# provenance for outputs to the vault
M.prov <- data.frame("info" = "provenance for output",
                     "path" = M.ofv_mft,
                     "object" = "M.mft",
                     "script" = name_scrpt,
                     "section" = "Section M - STEP 2",
                     "heading" = "create manifest",
                     stringsAsFactors = F)

### ************************************
### M - WRITE OUTPUTS ----
### ************************************

# outputs to the vault
write.csv(row.names = F, quote = F, x = M.mft, file = M.ofv_mft)

write.table(sep = "\t", row.names = F, x = M.info, file = M.ofv_info)
write.table(sep = "\t", row.names = F, x = M.prov, file = M.ofv_prov)

M.obj <- ls(pattern = "M.")
M.lst <- c(M.obj[grep(pattern = "M.", x = M.obj, ignore.case = F, fixed = T)],
           PREFACE.lst)
save(list = M.lst, file = M.ofv_wksp)

### ************************************
### Q - STEP 1 - import raw fastq seqs ----
### ************************************

# NOTE: section Q requires objects from the PREFACE to be in the environment

# ** note for KDP: QS import raw fastq seqs version 0.5 ** #

# make new directory to store sequence reads
mkdir reads/
  
# import FASTQ sequences into QIIME2
qiime tools import \
--input-path vault/manifest_R1.csv \
--type 'SampleData[SequencesWithQuality]' \
--input-format SingleEndFastqManifestPhred33 \
--output-path reads/reads.qza

# visualize summary statistics
qiime demux summarize \
--i-data reads/reads.qza \
--o-visualization reads/reads.qzv

# view
qiime tools view \
reads/reads.qzv

# quit viewer
q

# in reads.qzv:
# Per-sample sequence counts - Total Samples: 147
# Demultiplexed sequence counts summary - Total: 3,464,689

# in R Console run:
Q.S1_TotalSamples <- 147
Q.S1_TotalSequences <- 3464689

Q.info_S1a <- data.frame("info" = "number samples post import", 
                         value = Q.S1_TotalSamples, "object" = "reads.qza", 
                         "script" = name_scrpt, 
                         "section" = "SECTION Q - STEP 1",
                         "heading" = "import raw fastq seqs",
                         stringsAsFactors = F)
Q.info_S1b <- data.frame("info" = "number reads post import", 
                         value = Q.S1_TotalSequences, "object" = "reads.qza", 
                         "script" = name_scrpt, 
                         "section" = "SECTION Q - STEP 1",
                         "heading" = "import raw fastq seqs",
                         stringsAsFactors = F)
Q.info_S1 <- rbind(Q.info_S1a, Q.info_S1b)

### ************************************
### Q - STEP 2 - remove adapter sequences ----
### ************************************

# ** note for KDP: QS remove adapter sequences version 0.4 ** #

# remove non-biological sequences (e.g. adapters, 515F (Parada) primer, etc.)
# NOTE: this is a precautionary step as primers/adapters should not be present
# runtime: about 2 minutes
qiime cutadapt trim-single \
--i-demultiplexed-sequences reads/reads.qza \
--p-cores 3 \
--p-front GTGYCAGCMGCCGCGGTAA \
--p-error-rate 0.1 \
--p-match-adapter-wildcards \
--o-trimmed-sequences reads/reads_trimd.qza

# visualize summary statistics
qiime demux summarize \
--i-data reads/reads_trimd.qza \
--o-visualization reads/reads_trimd.qzv

# view
qiime tools view \
reads/reads_trimd.qzv

# quit viewer
q

# in reads_trimd.qzv:
# Per-sample sequence counts - Total Samples: 147
# Demultiplexed sequence counts summary - Total: 3,464,689

# in R Console run:
Q.S2_TotalSamples <- 147
Q.S2_TotalSequences <- 3464689

Q.info_S2a <- data.frame("info" = "number samples post cutadapt", 
                         value = Q.S2_TotalSamples, 
                         "object" = "reads_trimd.qzv", "script" = name_scrpt, 
                         "section" = "SECTION Q - STEP 2",
                         "heading" = "remove adapter sequences",
                         stringsAsFactors = F)
Q.info_S2b <- data.frame("info" = "number reads post cutadapt", 
                         value = Q.S2_TotalSequences, 
                         "object" = "reads_trimd.qzv", "script" = name_scrpt, 
                         "section" = "SECTION Q - STEP 2",
                         "heading" = "remove adapter sequences",
                         stringsAsFactors = F)
Q.info_S2 <- rbind(Q.info_S2a, Q.info_S2b)

### ************************************
### Q - STEP 3 - dada2 ----
### ************************************

# ** note for KDP: QS dada2 version 0.4 ** #

# NOTE: this step was performed using the normal mac terminal

# make new directory to store outputs
mkdir dada2/
  
# runtime: ~20 minutes
# --p-n-reads-learn 1000000 change from default = 4000000
qiime dada2 denoise-single \
--p-n-threads 0 \
--i-demultiplexed-seqs reads/reads_trimd.qza \
--p-trunc-len 248 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 1 \
--p-n-reads-learn 4000000 \
--o-table dada2/fts_tab.qza \
--o-representative-sequences dada2/rep_seq.qza \
--o-denoising-stats dada2/stats.qza

# visualize summary statistics
qiime feature-table summarize \
--i-table dada2/fts_tab.qza \
--o-visualization dada2/fts_tab.qzv

# view
qiime tools view \
dada2/fts_tab.qzv

# quit viewer
q

# in fts_tab.qzv: 'Table summary'
# Number of samples 147
# Number of features 1,219
# Total frequency 2,873,924

# NOTE: protocol control samples PCRnegC2 and NC1 not retained following dada2

# in R Console run:
Q.S3_NumberSamples <- 147
Q.S3_NumberFeatures <- 1219
Q.S3_TotalFrequency <- 2873924
Q.prov_secstep_QS3 <- "SECTION Q - STEP 3"
Q.prov_heading_QS3 <- "dada2"

Q.info_S3a <- data.frame("info" = "Number of samples", 
                         value = Q.S3_NumberSamples, "object" = "fts_tab.qza", 
                         "script" = name_scrpt, 
                         "section" = Q.prov_secstep_QS3,
                         "heading" = Q.prov_heading_QS3,
                         stringsAsFactors = F)
Q.info_S3b <- data.frame("info" = "Number of features", 
                         value = Q.S3_NumberFeatures, "object" = "fts_tab.qza", 
                         "script" = name_scrpt, 
                         "section" = Q.prov_secstep_QS3,
                         "heading" = Q.prov_heading_QS3,
                         stringsAsFactors = F)
Q.info_S3c <- data.frame("info" = "Total frequency", 
                         value = Q.S3_TotalFrequency, "object" = "fts_tab.qza", 
                         "script" = name_scrpt, 
                         "section" = Q.prov_secstep_QS3,
                         "heading" = Q.prov_heading_QS3,
                         stringsAsFactors = F)
Q.info_S3 <- rbind(Q.info_S3a, Q.info_S3b, Q.info_S3c)

### ************************************
### Q - STEP 4 - tabulate rep seqs ----
### ************************************

# ** note for KDP: QS tabulate rep seqs version 0.1 ** #

# tabulate seqs for easy 'BLASTing'
qiime feature-table tabulate-seqs \
--i-data dada2/rep_seq.qza \
--o-visualization dada2/rep_seq.qzv

### ************************************
### Q - STEP 5 - train 'in-house' taxonomy classifiers ----
### ************************************

# ** note for KDP: QS train 'in-house' taxonomy classifiers version Q2019.7 ** #

# NOTE: these classifiers use the 515F/806R (Parada/Apprill) primer pair

# make new directory to store outputs
mkdir class/
  
# import relevant files
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ../../gg_13_8_otus/rep_set/99_otus.fasta \
--output-path class/ggs_otu_99.qza

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ../../SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna \
--output-path class/slv_otu_99.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ../../gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
--output-path class/ggs_tax_99.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ../../SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt \
--output-path class/slv_tax_99.qza

# extract reference reads
# runtime: ~13 minutes
qiime feature-classifier extract-reads \
--i-sequences class/ggs_otu_99.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--p-trunc-len 248 \
--p-min-length 50 \
--p-max-length 0 \
--o-reads class/ggs_13_8_99_515Par_806App_248bp_ref_seqs_Q2019_7.qza

# runtime: ~25 minutes
qiime feature-classifier extract-reads \
--i-sequences class/slv_otu_99.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--p-trunc-len 248 \
--p-min-length 50 \
--p-max-length 0 \
--o-reads class/slv_132_99_515Par_806App_248bp_ref_seqs_Q2019_7.qza

# train the classifiers
# runtime: ~4 minutes
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads class/ggs_13_8_99_515Par_806App_248bp_ref_seqs_Q2019_7.qza \
--i-reference-taxonomy class/ggs_tax_99.qza \
--o-classifier class/ggs_13_8_99_515Par_806App_248bp_nb_classifier_Q2019_7.qza

# runtime: ~2 hours
# NOTE: this step was performed using the normal mac terminal
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads class/slv_132_99_515Par_806App_248bp_ref_seqs_Q2019_7.qza \
--i-reference-taxonomy class/slv_tax_99.qza \
--o-classifier class/slv_132_99_515Par_806App_248bp_nb_classifier_Q2019_7.qza

### ************************************
### Q - STEP 6 - assign taxonomy ----
### ************************************

# ** note for KDP: QS assign taxonomy version Q2019.7 ** #

# make new directory to store outputs
mkdir taxa/

# assign taxonomy with in-house (int) classifiers
# --p-confidence 0.7 change from default = 0.75
qiime feature-classifier classify-sklearn \
--p-n-jobs -1 \
--p-confidence 0.75 \
--i-reads dada2/rep_seq.qza \
--i-classifier class/ggs_13_8_99_515Par_806App_248bp_nb_classifier_Q2019_7.qza \
--o-classification taxa/int_ggs.qza

# runtime: ~20 minutes
# --p-confidence 0.7 change from default = 0.75
qiime feature-classifier classify-sklearn \
--p-n-jobs -1 \
--p-confidence 0.75 \
--i-reads dada2/rep_seq.qza \
--i-classifier class/slv_132_99_515Par_806App_248bp_nb_classifier_Q2019_7.qza \
--o-classification taxa/int_slv.qza

### ************************************
### Q - STEP 7 - generate phylogenetic tree ----
### ************************************

# ** note for KDP: QS generate phylogenetic tree version 0.1 ** #

# make new directory to store outputs
mkdir tree/
  
# perform multiple sequence alignment of rep seqs
qiime alignment mafft \
--i-sequences dada2/rep_seq.qza \
--o-alignment tree/rep_seq_aln.qza

# filter the alignment to remove any highly variable positions
qiime alignment mask \
--i-alignment tree/rep_seq_aln.qza \
--o-masked-alignment tree/rep_seq_aln_flt.qza

# create an unrooted phylogenetic tree
qiime phylogeny fasttree \
--i-alignment tree/rep_seq_aln_flt.qza \
--o-tree tree/unrt_tree.qza

# create a rooted phylogenetic tree
qiime phylogeny midpoint-root \
--i-tree tree/unrt_tree.qza \
--o-rooted-tree tree/root_tree.qza

### ************************************
### Q - STEP 8 - alpha diversity analysis ----
### ************************************

# ** note for KDP: QS alpha diversity analysis version 0.1 ** #

# make new directory
mkdir dvrs/
  
# alpha diversity via the chao 1 index
qiime diversity alpha \
--i-table dada2/fts_tab.qza \
--p-metric chao1 \
--o-alpha-diversity dvrs/chao.qza

# alpha diversity via the shannon's index
qiime diversity alpha \
--i-table dada2/fts_tab.qza \
--p-metric shannon \
--o-alpha-diversity dvrs/shan.qza

### ************************************
### Q - EXPORT ARTIFACTS from QIIME2 ----
### ************************************

# NOTE:: all ouputs to the "vault/" should have "q2_" appended to the file name

# feature table
qiime tools export \
--input-path dada2/fts_tab.qza \
--output-path dada2/
  
# convert from biom format to classic format (tab-delimited) for use with R
biom convert -i dada2/feature-table.biom -o dada2/fts_tab.tsv --to-tsv

# copy feature table to the "vault/"
cp dada2/fts_tab.tsv vault/q2_fts_tab.tsv

# representative sequences
qiime tools export \
--input-path dada2/rep_seq.qza \
--output-path dada2/
  
# copy representative sequences to the "vault/" (and rename)
cp dada2/dna-sequences.fasta vault/q2_rep_seq.fasta

# in-house trained greengenes taxonomy table
qiime tools export \
--input-path taxa/int_ggs.qza \
--output-path taxa/

# copy to rename, and then copy greengenes taxonomy table to the "vault/"
cp taxa/taxonomy.tsv taxa/int_ggs.tsv
cp taxa/int_ggs.tsv vault/q2_int_ggs.tsv

# in-house trained silva taxonomy table
qiime tools export \
--input-path taxa/int_slv.qza \
--output-path taxa/

# copy to rename, then copy silva taxonomy table to the "vault/"
cp taxa/taxonomy.tsv taxa/int_slv.tsv
cp taxa/int_slv.tsv vault/q2_int_slv.tsv

# rooted phylogenetic tree
qiime tools export \
--input-path tree/root_tree.qza \
--output-path tree/

# copy rooted phylogenetic tree to the "vault/" (and rename)
cp tree/tree.nwk vault/q2_tree.nwk

# alpha diversity via the chao 1 index
qiime tools export \
--input-path dvrs/chao.qza \
--output-path dvrs/
  
# copy to rename, then copy table to the "vault/"
cp dvrs/alpha-diversity.tsv dvrs/chao.tsv
cp dvrs/chao.tsv vault/q2_chao.tsv

# alpha diversity via via shannon's index
qiime tools export \
--input-path dvrs/shan.qza \
--output-path dvrs/
  
# copy to rename, then copy table to the "vault/"
cp dvrs/alpha-diversity.tsv dvrs/shan.tsv
cp dvrs/shan.tsv vault/q2_shan.tsv

### ************************************
### Q - PROVENANCE FOR EXPORTED QIIME ARTIFACTS ----
### ************************************

# in R Console run:
Q.prov_output1 <- data.frame("info" = "provenance for output",
                             "path" = paste(sep = "", path_vault, 
                                            "/q2_fts_tab.tsv"),
                             "object" = "q2_fts_tab.tsv",
                             "lineage" = paste(sep = " ", 
                                               "q2_fts_tab.tsv",
                                               "<- $ cp", 
                                               "fts_tab.tsv", 
                                               "<- $ biom convert", 
                                               "feature-table.biom", 
                                               "<- $ export",
                                               "fts_tab.qza",
                                               "<- $ dada2 denoise-single",
                                               "reads_trimd.qza",
                                               "<- $ cutadapt trim-single",
                                               "reads.qza", 
                                               "<- $ import",
                                               "raw_fastq/"),
                             "script" = name_scrpt,
                             "section" = Q.prov_secstep_QS3,
                             "heading" = Q.prov_heading_QS3,
                             stringsAsFactors = F)

Q.prov_output2 <- data.frame("info" = "provenance for output",
                             "path" = paste(sep = "", path_vault, 
                                            "/q2_rep_seq.fasta"),
                             "object" = "q2_rep_seq.fasta",
                             "lineage" = paste(sep = " ", 
                                               "q2_rep_seq.fasta",
                                               "<- $ cp", 
                                               "dna-sequences.fasta", 
                                               "<- $ export",
                                               "rep_seq.qza",
                                               "<- $ dada2 denoise-single",
                                               "reads_trimd.qza",
                                               "<- $ cutadapt trim-single",
                                               "reads.qza", 
                                               "<- $ import",
                                               "raw_fastq/"),
                             "script" = name_scrpt,
                             "section" = Q.prov_secstep_QS3,
                             "heading" = Q.prov_heading_QS3,
                             stringsAsFactors = F)

Q.prov_output3 <- data.frame("info" = "provenance for output",
                             "path" = paste(sep = "", path_vault, 
                                            "/q2_int_ggs.tsv"),
                             "object" = "q2_int_ggs.tsv",
                             "lineage" = paste(sep = " ", 
                                               "q2_int_ggs.tsv",
                                               "<- $ cp",
                                               "int_ggs.tsv",
                                               "<- $ cp",
                                               "taxonomy.tsv",
                                               "<- $ export", 
                                               "int_ggs.qza",
                                               "<- $ feature-classifier",
                                               "classify-sklearn",
                                               "rep_seq.qza",
                                               "<- $ dada2 denoise-single",
                                               "reads_trimd.qza",
                                               "<- $ cutadapt trim-single",
                                               "reads.qza", 
                                               "<- $ import",
                                               "raw_fastq/"),
                             "script" = name_scrpt,
                             "section" = "SECTION Q - STEP 6",
                             "heading" = "assign taxonomy",
                             stringsAsFactors = F)

Q.prov_output4 <- data.frame("info" = "provenance for output",
                             "path" = paste(sep = "", path_vault, 
                                            "/q2_int_slv.tsv"),
                             "object" = "q2_int_slv.tsv",
                             "lineage" = paste(sep = " ", 
                                               "q2_int_slv.tsv",
                                               "<- $ cp",
                                               "int_slv.tsv",
                                               "<- $ cp",
                                               "taxonomy.tsv",
                                               "<- $ export", 
                                               "int_slv.qza",
                                               "<- $ feature-classifier",
                                               "classify-sklearn",
                                               "rep_seq.qza",
                                               "<- $ dada2 denoise-single",
                                               "reads_trimd.qza",
                                               "<- $ cutadapt trim-single",
                                               "reads.qza", 
                                               "<- $ import",
                                               "raw_fastq/"),
                             "script" = name_scrpt,
                             "section" = "SECTION Q - STEP 6",
                             "heading" = "assign taxonomy",
                             stringsAsFactors = F)

Q.prov_output5 <- data.frame("info" = "provenance for output",
                             "path" = paste(sep = "", path_vault, 
                                            "/q2_tree.nwk"),
                             "object" = "q2_tree.nwk",
                             "lineage" = paste(sep = " ", 
                                               "q2_tree.nwk",
                                               "<- $ cp",
                                               "tree.nwk",
                                               "<- $ export", 
                                               "root_tree.qza",
                                               "<- $ phylogeny midpoint-root",
                                               "unrt_tree.qza",
                                               "<- $ phylogeny fasttree",
                                               "rep_seq_aln_flt.qza",
                                               "<- $ alignment mask",
                                               "rep_seq_aln.qza",
                                               "<- $ alignment mafft",
                                               "rep_seq.qza",
                                               "<- $ dada2 denoise-single",
                                               "reads_trimd.qza",
                                               "<- $ cutadapt trim-single",
                                               "reads.qza", 
                                               "<- $ import",
                                               "raw_fastq/"),
                             "script" = name_scrpt,
                             "section" = "SECTION Q - STEP 7",
                             "heading" = "generate phylogenetic tree",
                             stringsAsFactors = F)

Q.prov_output6 <- data.frame("info" = "provenance for output",
                             "path" = paste(sep = "", path_vault, 
                                            "/q2_chao.tsv"),
                             "object" = "q2_chao.tsv",
                             "lineage" = paste(sep = " ", 
                                               "q2_chao.tsv",
                                               "<- $ cp", 
                                               "chao.tsv",
                                               "<- $ cp", 
                                               "alpha-diversity.tsv",
                                               "<- $ export",
                                               "chao.qza",
                                               "<- $ diversity alpha",
                                               "fts_tab.qza",
                                               "<- $ dada2 denoise-single",
                                               "reads_trimd.qza",
                                               "<- $ cutadapt trim-single",
                                               "reads.qza", 
                                               "<- $ import",
                                               "raw_fastq/"),
                             "script" = name_scrpt,
                             "section" = "SECTION Q - STEP 8",
                             "heading" = "alpha diversity analysis",
                             stringsAsFactors = F)

Q.prov_output7 <- data.frame("info" = "provenance for output",
                             "path" = paste(sep = "", path_vault, 
                                            "/q2_shan.tsv"),
                             "object" = "q2_shan.tsv",
                             "lineage" = paste(sep = " ", 
                                               "q2_shan.tsv",
                                               "<- $ cp", 
                                               "shan.tsv",
                                               "<- $ cp", 
                                               "alpha-diversity.tsv",
                                               "<- $ export",
                                               "shan.qza",
                                               "<- $ diversity alpha",
                                               "fts_tab.qza",
                                               "<- $ dada2 denoise-single",
                                               "reads_trimd.qza",
                                               "<- $ cutadapt trim-single",
                                               "reads.qza", 
                                               "<- $ import",
                                               "raw_fastq/"),
                             "script" = name_scrpt,
                             "section" = "SECTION Q - STEP 8",
                             "heading" = "alpha diversity analysis",
                             stringsAsFactors = F)

### ************************************
### Q - WRITE OUTPUTS from R ----
### ************************************

# rbind Q.info data.frames together
Q.info <- rbind(Q.info_S1, Q.info_S2, Q.info_S3)

# rbind Q.prov data.frames together
Q.prov <- rbind(Q.prov_output1, Q.prov_output2, Q.prov_output3, Q.prov_output4,
                Q.prov_output5, Q.prov_output6, Q.prov_output7)

# outputs to the vault
write.table(sep = "\t", row.names = F, x = Q.info, file = Q.ofv_info)
write.table(sep = "\t", row.names = F, x = Q.prov, file = Q.ofv_prov)

Q.obj <- ls(pattern = "Q.")
Q.lst <- c(Q.obj[grep(pattern = "Q.", x = Q.obj, ignore.case = F, fixed = T)],
           PREFACE.lst)
save(list = Q.lst, file = Q.ofv_wksp)

# proceed to next script
