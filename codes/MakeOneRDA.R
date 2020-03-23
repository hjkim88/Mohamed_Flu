###
#   File name : MakeOneRDA.R
#   Author    : Hyunjin Kim
#   Date      : Mar 20, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Mohamed has the cytokine level data and patients' sample info.
#               This script will load the data and make them into R objects,
#               then save all of them in one RDA file.
#
#   Instruction
#               1. Source("MakeOneRDA.R")
#               2. Run the function "makeRDA" - specify the necessary input paths and the output RDA file path
#               3. The RDA file will be generated in the output path
#
#   Example
#               > source("The_directory_of_MakeOneRDA.R/MakeOneRDA.R")
#               > makeRDA(cytokine_file_path="./data/Copy of FLU09 Cytokines_season 2009-2014_17July2014.xlsx",
#                         sample_info_file_path="./data/Copy of FLU09_Demographics_21July2014.xlsx",
#                         output_rda_path="./data/flu09_cytokine.rda")
###

makeRDA <- function(cytokine_file_path="./data/Copy of FLU09 Cytokines_season 2009-2014_17July2014.xlsx",
                    sample_info_file_path="./data/Copy of FLU09_Demographics_21July2014.xlsx",
                    output_rda_path="./data/flu09_cytokine.rda") {
  
  ### load library
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### load data
  cyto_nw <- read.xlsx2(file = cytokine_file_path, sheetName = "NW Cytokines",
                        colClasses = c(rep("character", 5), rep("numeric", 42)),
                        stringsAsFactors = FALSE, check.names = FALSE)
  cyto_plasma <- read.xlsx2(file = cytokine_file_path, sheetName = "Plasma Cytokines",
                            colClasses = c(rep("character", 6), rep("numeric", 39)),
                            stringsAsFactors = FALSE, check.names = FALSE)
  cyto_sample <- read.xlsx2(file = sample_info_file_path, sheetIndex = 1,
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  ### change IL17A in cyto_nw to IL17
  ### they are the same thing but for consistency
  colnames(cyto_nw)[which(colnames(cyto_nw) == "IL17A")] <- "IL17"
  
  ### remove empty rows in the sample info
  cyto_sample <- cyto_sample[which(cyto_sample$ID != ""),]
  
  ### set ID as the row names
  rownames(cyto_sample) <- cyto_sample$ID
  
  ### set README function
  README <- function() {
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"cyto_nw\" are cytokine levels data from nasal wash.")
    writeLines("The \"cyto_plasma\" are cytokine levels data from PBMC.")
    writeLines("Each patient has some samples that were taken at different time points.")
    writeLines("The \"cyto_sample\" has the patient information.")
    writeLines("")
    writeLines("The original data is from Mohamed Ghonim (Mohamed.Ghonim@STJUDE.ORG),")
    writeLines("and is combined into this one RDA file by Hyunjin Kim (Hyunjin.Kim@STJUDE.ORG).")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save the objects into one RDA file
  save(list = c("cyto_nw", "cyto_plasma", "cyto_sample", "README"),
       file = output_rda_path)
  
}
