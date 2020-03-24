###
#   File name : Analysis1.R
#   Author    : Hyunjin Kim
#   Date      : Mar 23, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : First analysis of the Mohamed's flu09 data.
#               It has 4 tasks to do:
#               1. Explore the data set. Find any association between cytokine levels
#                  and various clinical factors. This includes the plot per cytokine to put
#                  all the samples to compare Asthma/IV and IV. Different time points are
#                  painted with different colored dots.
#               2. PCA. Mark with Asthma info and with time point info to see
#                  if the samples are clustered well.
#               3. a) Line graph per each cytokine. The x-axis is time points and the y-axis is
#                     cytokine level. Color the graphs based on Asthma/IV vs IV, so that we can
#                     distinguish them in the plot.
#                  b) Line graph per each cytokine. Take mean of each time point and compare
#                     Asthma/IV vs IV. There would be only two lines in each graph.
#               4. Statistics table. For each time point, calculate mean/median difference
#                  between Asthma/IV vs IV, and perform t-test and one-way ANOVA to get p-values.
#
#               * Here, it is not about difference of trends but to know the sense of actual
#                 difference in cytokine levels in each time point or as a whole. So we
#                 are not doing a time-series analysis.
#
#   Instruction
#               1. Source("Analysis1.R")
#               2. Run the function "flu09_analysis1" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Analysis1.R/Analysis1.R")
#               > flu09_analysis1(data_path="./data/flu09_cytokine.rda",
#                                 output_dir="./results/analysis1/")
###

flu09_analysis1 <- function(data_path="./data/flu09_cytokine.rda",
                            output_dir="./results/analysis1/") {
  
  ### load libraries
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
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
  
  
  ### 1. Explore the data set.
  
  ### create result directory
  outDir <- paste0(output_dir, "1_Explore/")
  dir.create(path = outDir, showWarnings = FALSE)
  
  ### clinical factors to be tested
  ### must be in colnames(cyto_sample)
  factor_list <- c("FLU.Strain.Designation",
                   "Enrollment.Status",
                   "Age.at.Enrollment",
                   "Race",
                   "Gender",
                   "Respiratory.Disease",
                   "Tobacco.Use")
  
  ### data to be tested
  plot_df <- cyto_sample[,factor_list]
  
  ### numerize the numeric column
  plot_df[,"Age.at.Enrollment"] <- as.numeric(plot_df[,"Age.at.Enrollment"])
  cyto_nw[,"Study.Day"] <- as.numeric(cyto_nw[,"Study.Day"])
  cyto_plasma[,"Study.Day"] <- as.numeric(cyto_plasma[,"Study.Day"])
  
  ### annotate clinical info to the cytokine level data
  plot_df_nw <- merge(cyto_nw[,c("ID", "Study.Day", th1_cytokines, th2_cytokines)], plot_df,
                      by.x = "ID", by.y = "row.names")
  plot_df_ps <- merge(cyto_plasma[,c("ID", "Study.Day", th1_cytokines, th2_cytokines)], plot_df,
                      by.x = "ID", by.y = "row.names")
  
  ### if discrete: make beeswarm plots
  ### if continous: make correlation plots
  for(factor in factor_list) {
    
    ### correlation
    if(class(plot_df[,factor]) == "numeric") {
      
      ### NW
      
      ### TH1 cytokines plot
      p <- vector("list", length = length(th1_cytokines))
      names(p) <- th1_cytokines
      
      ### correlation plot per cytokine
      for(cytokine in th1_cytokines) {
        p[[cytokine]] <- ggplot(data = plot_df_nw, aes_string(x=factor, y=cytokine)) +
          geom_point(aes_string(col="FLU.Strain.Designation"), size = 2) +
          labs(subtitle=paste("Pearson Correlation = ", round(cor(plot_df_nw[,factor],
                                                                 plot_df_nw[,cytokine],
                                                                 use = "pairwise.complete.obs"), 5),
                              "\nP-value = ", signif(cor.test(plot_df_nw[,factor], plot_df_nw[,cytokine])$p.value, 5))) +
          xlab(factor) +
          ylab(cytokine) +
          geom_smooth(method = lm, color="black", se=FALSE) +
          theme_classic(base_size = 16)
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_MCP-Counter_Result_", group)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4, 5), c(6, 7, 8, 9, 10)),
                       top = fName)
      ggsave(file = paste0(outputDir, fName, ".png"), g, width = 22, height = 12)
      
    
    ### beeswarm
    } else if(class(plot_df[,factor]) == "character") {
      
      ggplot(cyto_sample, aes_string(x="Survival_Status", y="logFDR")) +
        theme_classic(base_size = 16) +
        geom_boxplot() +
        geom_beeswarm(aes_string(color="Survival_Status"), na.rm = TRUE) +
        stat_compare_means(label.y = max_y*1.1) +
        ggtitle(paste0("ECHs on Viper(", tissue, ") -log10(GSEA FDR) between Good and Bad Survival")) +
        labs(x = paste0("Survival Status Top-Bottom ", params[[2]], "%"), y = "-log10(GSEA FDR)") +
        theme(legend.position = "None")
      ggsave(filename = paste0(result_dir, "beeswarm_plot_FDR_", tissue, ".png"), width = 12, height = 10)
      
    } else {
      writeLines(paste0("ERROR: class(plot_df[,", factor, "]) is not character nor numeric- ", class(plot_df[,factor])))
    }
    
    
  }
  
  
  
  
  
  
    
}
