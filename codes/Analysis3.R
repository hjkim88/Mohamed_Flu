###
#   File name : Analysis3.R
#   Author    : Hyunjin Kim
#   Date      : Apr 9, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : There may be a concern that age can affect the cytokine analysis.
#               e.g., we see cytokine level difference between IV+ and IV- but
#               it can happen due to age difference.
#               Therefore, I would like to make pie charts to see age differences
#               among the groups and also generate age-corrected cytokine level data,
#               then will produce line graphs with the age-corrected data.
#
#   Instruction
#               1. Source("Analysis3.R")
#               2. Run the function "flu09_analysis3" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Analysis3.R/Analysis3.R")
#               > flu09_analysis3(data_path="./data/flu09_cytokine.rda",
#                                 output_dir="./results/analysis3/")
###

flu09_analysis3 <- function(data_path="./data/flu09_cytokine.rda",
                            output_dir="./results/analysis3/") {
  
  ### load libraries
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  
  ### load the data
  load(data_path)
  
  ### set types of cytokines
  ### defined by the Mohamed's presentation file
  th1_cytokines <- c("IL10",
                     "GRO",
                     "IL8",
                     "GCSF",
                     "IFNa2",
                     "FKN",
                     "MCP1",
                     "MCP3",
                     "IL1a",
                     "IL1Ra2",
                     "TNFa",
                     "VEGF")
  th2_cytokines <- c("IL4",
                     "IL5",
                     "IL6",
                     "IL9",
                     "IL13",
                     "IL15",
                     "IL17")
  
  ### clinical factors to be tested
  ### must be in colnames(cyto_sample)
  factor_list <- c("FLU.Strain.Designation",
                   "Age.at.Enrollment",
                   "Race",
                   "Gender",
                   "Respiratory.Disease",
                   "Tobacco.Use")
  
  ### data to be tested
  plot_df <- cyto_sample[,factor_list]
  
  ### annotate additional info - IV+/Resp+, IV+/Resp-, IV-/Resp+, IV-/Resp-
  plot_df$IV.Resp <- sapply(1:nrow(plot_df), function(x) {
    if(plot_df[x,"FLU.Strain.Designation"] == "Negative" && plot_df[x,"Respiratory.Disease"] == "Yes") {
      return("IV-/Resp+")
    } else if(plot_df[x,"FLU.Strain.Designation"] == "Negative" && plot_df[x,"Respiratory.Disease"] == "No") {
      return("IV-/Resp-")
    } else if(plot_df[x,"FLU.Strain.Designation"] != "Negative" && plot_df[x,"Respiratory.Disease"] == "Yes") {
      return("IV+/Resp+")
    } else if(plot_df[x,"FLU.Strain.Designation"] != "Negative" && plot_df[x,"Respiratory.Disease"] == "No"){
      return("IV+/Resp-")
    } else {
      return(NA)
    }
  })
  
  ### numerize the numeric column
  plot_df[,"Age.at.Enrollment"] <- as.numeric(plot_df[,"Age.at.Enrollment"])
  cyto_nw[,"Study.Day"] <- as.numeric(cyto_nw[,"Study.Day"])
  cyto_plasma[,"Study.Day"] <- as.numeric(cyto_plasma[,"Study.Day"])
  
  ### N/A -> NA in Tobacco.Use
  plot_df[which(plot_df[,"Tobacco.Use"] == "N/A"),"Tobacco.Use"] <- "NA"
  
  ### annotate clinical info to the cytokine level data
  plot_df_nw <- merge(cyto_nw[,c("ID", "Study.Day", th1_cytokines, th2_cytokines)], plot_df,
                      by.x = "ID", by.y = "row.names")
  plot_df_ps <- merge(cyto_plasma[,c("ID", "Study.Day", th1_cytokines, th2_cytokines)], plot_df,
                      by.x = "ID", by.y = "row.names")
  
  ### remove samples that are categorized into a group with less than 10 samples
  for(factor in setdiff(c(factor_list, "Study.Day"), "Age.at.Enrollment")) {
    ### NW
    unique_vals <- unique(plot_df_nw[,factor])
    remove_ind <- NULL
    for(val in unique_vals) {
      target_ind <- which(plot_df_nw[,factor] == val)
      if(length(target_ind) < 10) {
        remove_ind <- c(remove_ind, target_ind)
      }
    }
    if(length(remove_ind) > 0) {
      plot_df_nw <- plot_df_nw[-remove_ind,]
    }
    
    ### plasma
    unique_vals <- unique(plot_df_ps[,factor])
    remove_ind <- NULL
    for(val in unique_vals) {
      target_ind <- which(plot_df_ps[,factor] == val)
      if(length(target_ind) < 10) {
        remove_ind <- c(remove_ind, target_ind)
      }
    }
    if(length(remove_ind) > 0) {
      plot_df_ps <- plot_df_ps[-remove_ind,]
    }
  }
  
  ### add new column - Age.over.50
  plot_df_nw$Age.over.50 <- "FALSE"
  plot_df_nw$Age.over.50[plot_df_nw$Age.at.Enrollment >= 50] <- "TRUE"
  plot_df_ps$Age.over.50 <- "FALSE"
  plot_df_ps$Age.over.50[plot_df_ps$Age.at.Enrollment >= 50] <- "TRUE"
  
  
  ### 1. pie charts
  
  ### create result directory
  outDir <- paste0(output_dir, "1_Pie_Chart/")
  dir.create(path = outDir, showWarnings = FALSE, recursive = TRUE)
  
  ### NW
  ### create an analysis-ready data for the pie charts
  pie_data <- vector("list", length = length(unique(plot_df_nw$IV.Resp)))
  names(pie_data) <- unique(plot_df_nw$IV.Resp)
  for(group in unique(plot_df_nw$IV.Resp)) {
    pie_data[[group]] <- data.frame(Age.over.50=c("FALSE", "TRUE"),
                                    n=c(length(intersect(which(plot_df_nw$IV.Resp == group),
                                                         which(plot_df_nw$Age.over.50 == "FALSE"))),
                                        length(intersect(which(plot_df_nw$IV.Resp == group),
                                                         which(plot_df_nw$Age.over.50 == "TRUE")))),
                                    stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  ### get pie charts
  p <- vector("list", length = length(unique(plot_df_nw$IV.Resp)))
  names(p) <- unique(plot_df_nw$IV.Resp)
  for(group in unique(plot_df_nw$IV.Resp)) {
    total_sample_num <- sum(pie_data[[group]]$n)
    pct <- sapply(pie_data[[group]]$n, function(x) signif(x*100/total_sample_num, digits = 3))
    pct <- paste0(pie_data[[group]]$n, "(", pct, "%)")
    p[[group]] <- ggplot(data = pie_data[[group]], aes(x = "", y = n, fill = Age.over.50)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta="y") +
      geom_text(label = pct, position = position_stack(vjust = 0.5), size = 6) +
      labs(x = NULL, y = NULL, title = group) +
      theme_classic(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, color = "black"),
            axis.line = element_blank(),
            axis.ticks = element_blank())
  }
  
  ### arrange the plots and print out
  fName <- paste0("Pie_Chart_NW_IV.Resp_Age")
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2,
                   top = fName)
  ggsave(file = paste0(outDir, fName, ".png"), g, width = 15, height = 9, dpi = 300)
  
  ### PS
  ### create an analysis-ready data for the pie charts
  pie_data <- vector("list", length = length(unique(plot_df_ps$IV.Resp)))
  names(pie_data) <- unique(plot_df_ps$IV.Resp)
  for(group in unique(plot_df_ps$IV.Resp)) {
    pie_data[[group]] <- data.frame(Age.over.50=c("FALSE", "TRUE"),
                                    n=c(length(intersect(which(plot_df_ps$IV.Resp == group),
                                                         which(plot_df_ps$Age.over.50 == "FALSE"))),
                                        length(intersect(which(plot_df_ps$IV.Resp == group),
                                                         which(plot_df_ps$Age.over.50 == "TRUE")))),
                                    stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  ### get pie charts
  p <- vector("list", length = length(unique(plot_df_ps$IV.Resp)))
  names(p) <- unique(plot_df_ps$IV.Resp)
  for(group in unique(plot_df_ps$IV.Resp)) {
    total_sample_num <- sum(pie_data[[group]]$n)
    pct <- sapply(pie_data[[group]]$n, function(x) signif(x*100/total_sample_num, digits = 3))
    pct <- paste0(pie_data[[group]]$n, "(", pct, "%)")
    p[[group]] <- ggplot(data = pie_data[[group]], aes(x = "", y = n, fill = Age.over.50)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta="y") +
      geom_text(label = pct, position = position_stack(vjust = 0.5), size = 6) +
      labs(x = NULL, y = NULL, title = group) +
      theme_classic(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, color = "black"),
            axis.line = element_blank(),
            axis.ticks = element_blank())
  }
  
  ### arrange the plots and print out
  fName <- paste0("Pie_Chart_PS_IV.Resp_Age")
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2,
                   top = fName)
  ggsave(file = paste0(outDir, fName, ".png"), g, width = 15, height = 9, dpi = 300)
  
}
