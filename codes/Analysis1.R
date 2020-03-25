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
  
  # ### remove flu- samples since we are not interested in them in this analysis
  # plot_df <- plot_df[which(plot_df[,"FLU.Strain.Designation"] != "Negative"),]
  
  ### numerize the numeric column
  plot_df[,"Age.at.Enrollment"] <- as.numeric(plot_df[,"Age.at.Enrollment"])
  cyto_nw[,"Study.Day"] <- as.numeric(cyto_nw[,"Study.Day"])
  cyto_plasma[,"Study.Day"] <- as.numeric(cyto_plasma[,"Study.Day"])
  
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
  
  ### if discrete: make beeswarm plots
  ### if continous: make correlation plots
  for(factor in factor_list) {
    
    ### create result directory
    outDir <- paste0(output_dir, "1_Explore/", factor, "/")
    dir.create(path = outDir, showWarnings = FALSE, recursive = TRUE)
    
    ### correlation
    if(class(plot_df[,factor]) == "numeric") {
      
      ### NW
      
      ### TH1 cytokines plot
      p <- vector("list", length = length(th1_cytokines))
      names(p) <- th1_cytokines
      
      ### correlation plot per cytokine
      for(cytokine in th1_cytokines) {
        p[[cytokine]] <- ggplot(data = plot_df_nw, aes_string(x=factor, y=cytokine)) +
          geom_point(aes_string(col="Respiratory.Disease"), size = 2) +
          labs(title=paste("Pearson Correlation = ", round(cor(plot_df_nw[,factor],
                                                               plot_df_nw[,cytokine],
                                                               use = "pairwise.complete.obs"), 5),
                              "\nP-value = ", signif(cor.test(plot_df_nw[,factor], plot_df_nw[,cytokine])$p.value, 5))) +
          xlab(factor) +
          ylab(cytokine) +
          geom_smooth(method = lm, color="black", se=FALSE) +
          theme_classic(base_size = 16) +
          theme(legend.position = c(0.5, 1),
                legend.justification = c("top"),
                legend.direction = "horizontal",
                legend.title = element_text(size=12),
                legend.text = element_text(size=10),
                plot.title=element_text(size=12, hjust=0.5, margin=margin(b=-8)))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Correlation_Plot_NW_TH1_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### TH2 cytokines plot
      p <- vector("list", length = length(th2_cytokines))
      names(p) <- th2_cytokines
      
      ### correlation plot per cytokine
      for(cytokine in th2_cytokines) {
        p[[cytokine]] <- ggplot(data = plot_df_nw, aes_string(x=factor, y=cytokine)) +
          geom_point(aes_string(col="Respiratory.Disease"), size = 2) +
          labs(title=paste("Pearson Correlation = ", round(cor(plot_df_nw[,factor],
                                                               plot_df_nw[,cytokine],
                                                               use = "pairwise.complete.obs"), 5),
                           "\nP-value = ", signif(cor.test(plot_df_nw[,factor], plot_df_nw[,cytokine])$p.value, 5))) +
          xlab(factor) +
          ylab(cytokine) +
          geom_smooth(method = lm, color="black", se=FALSE) +
          theme_classic(base_size = 16) +
          theme(legend.position = c(0.5, 1),
                legend.justification = c("top"),
                legend.direction = "horizontal",
                legend.title = element_text(size=12),
                legend.text = element_text(size=10),
                plot.title=element_text(size=12, hjust=0.5, margin=margin(b=-8)))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Correlation_Plot_NW_TH2_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      
      ### Plasma
      
      ### TH1 cytokines plot
      p <- vector("list", length = length(th1_cytokines))
      names(p) <- th1_cytokines
      
      ### correlation plot per cytokine
      for(cytokine in th1_cytokines) {
        p[[cytokine]] <- ggplot(data = plot_df_ps, aes_string(x=factor, y=cytokine)) +
          geom_point(aes_string(col="Respiratory.Disease"), size = 2) +
          labs(title=paste("Pearson Correlation = ", round(cor(plot_df_ps[,factor],
                                                               plot_df_ps[,cytokine],
                                                               use = "pairwise.complete.obs"), 5),
                           "\nP-value = ", signif(cor.test(plot_df_ps[,factor], plot_df_ps[,cytokine])$p.value, 5))) +
          xlab(factor) +
          ylab(cytokine) +
          geom_smooth(method = lm, color="black", se=FALSE) +
          theme_classic(base_size = 16) +
          theme(legend.position = c(0.5, 1),
                legend.justification = c("top"),
                legend.direction = "horizontal",
                legend.title = element_text(size=12),
                legend.text = element_text(size=10),
                plot.title=element_text(size=12, hjust=0.5, margin=margin(b=-8)))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Correlation_Plot_PS_TH1_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### TH2 cytokines plot
      p <- vector("list", length = length(th2_cytokines))
      names(p) <- th2_cytokines
      
      ### correlation plot per cytokine
      for(cytokine in th2_cytokines) {
        p[[cytokine]] <- ggplot(data = plot_df_ps, aes_string(x=factor, y=cytokine)) +
          geom_point(aes_string(col="Respiratory.Disease"), size = 2) +
          labs(title=paste("Pearson Correlation = ", round(cor(plot_df_ps[,factor],
                                                               plot_df_ps[,cytokine],
                                                               use = "pairwise.complete.obs"), 5),
                           "\nP-value = ", signif(cor.test(plot_df_ps[,factor], plot_df_ps[,cytokine])$p.value, 5))) +
          xlab(factor) +
          ylab(cytokine) +
          geom_smooth(method = lm, color="black", se=FALSE) +
          theme_classic(base_size = 16) +
          theme(legend.position = c(0.5, 1),
                legend.justification = c("top"),
                legend.direction = "horizontal",
                legend.title = element_text(size=12),
                legend.text = element_text(size=10),
                plot.title=element_text(size=12, hjust=0.5, margin=margin(b=-8)))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Correlation_Plot_PS_TH2_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
    ### beeswarm
    } else if(class(plot_df[,factor]) == "character") {
      
      ### NW
      
      ### TH1 cytokines plot
      p <- vector("list", length = length(th1_cytokines))
      names(p) <- th1_cytokines
      
      ### beeswarm plot per cytokine
      for(cytokine in th1_cytokines) {
        p[[cytokine]] <- ggplot(plot_df_nw, aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_NW_TH1_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### because of the outliers, the plot can look messy
      ### therefore, only using q1-q3 data
      ### beeswarm plot per cytokine
      for(cytokine in th1_cytokines) {
        q <- quantile(plot_df_nw[,cytokine], na.rm = TRUE)
        
        p[[cytokine]] <- ggplot(plot_df_nw[intersect(which(plot_df_nw[,cytokine] >= q[2]),
                                                     which(plot_df_nw[,cytokine] <= q[3])),],
                                aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_NW_TH1_Q1-Q3_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### TH2 cytokines plot
      p <- vector("list", length = length(th2_cytokines))
      names(p) <- th2_cytokines
      
      ### beeswarm plot per cytokine
      for(cytokine in th2_cytokines) {
        p[[cytokine]] <- ggplot(plot_df_nw, aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_NW_TH2_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### because of the outliers, the plot can look messy
      ### therefore, only using q1-q3 data
      ### beeswarm plot per cytokine
      for(cytokine in th2_cytokines) {
        q <- quantile(plot_df_nw[,cytokine], na.rm = TRUE)
        
        p[[cytokine]] <- ggplot(plot_df_nw[intersect(which(plot_df_nw[,cytokine] >= q[2]),
                                                     which(plot_df_nw[,cytokine] <= q[3])),],
                                aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_NW_TH2_Q1-Q3_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      
      ### Plasma
      
      ### TH1 cytokines plot
      p <- vector("list", length = length(th1_cytokines))
      names(p) <- th1_cytokines
      
      ### beeswarm plot per cytokine
      for(cytokine in th1_cytokines) {
        p[[cytokine]] <- ggplot(plot_df_ps, aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_PS_TH1_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### because of the outliers, the plot can look messy
      ### therefore, only using q1-q3 data
      ### beeswarm plot per cytokine
      for(cytokine in th1_cytokines) {
        q <- quantile(plot_df_ps[,cytokine], na.rm = TRUE)
        
        p[[cytokine]] <- ggplot(plot_df_ps[intersect(which(plot_df_ps[,cytokine] >= q[2]),
                                                     which(plot_df_ps[,cytokine] <= q[3])),],
                                aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_PS_TH1_Q1-Q3_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### TH2 cytokines plot
      p <- vector("list", length = length(th2_cytokines))
      names(p) <- th2_cytokines
      
      ### beeswarm plot per cytokine
      for(cytokine in th2_cytokines) {
        p[[cytokine]] <- ggplot(plot_df_ps, aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_PS_TH2_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
      ### because of the outliers, the plot can look messy
      ### therefore, only using q1-q3 data
      ### beeswarm plot per cytokine
      for(cytokine in th2_cytokines) {
        q <- quantile(plot_df_ps[,cytokine], na.rm = TRUE)
        
        p[[cytokine]] <- ggplot(plot_df_ps[intersect(which(plot_df_ps[,cytokine] >= q[2]),
                                                     which(plot_df_ps[,cytokine] <= q[3])),],
                                aes_string(x=factor, y=cytokine)) +
          geom_boxplot() +
          geom_beeswarm(aes_string(color="Study.Day"), na.rm = TRUE) +
          stat_compare_means() +
          theme_classic(base_size = 16) +
          theme(legend.title = element_text(size=10),
                legend.text = element_text(size=8))
      }
      
      ### arrange the plots and print out
      fName <- paste0("Beeswarm_Plot_PS_TH2_Q1-Q3_", factor)
      g <- arrangeGrob(grobs = p,
                       layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)),
                       top = fName)
      ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
      
    } else {
      writeLines(paste0("ERROR: class(plot_df[,", factor, "]) is not character nor numeric- ", class(plot_df[,factor])))
    }
    
  }
  
  
  ### 2. PCA
  
  ### A function to perform 2D PCA and save a plot
  ### mat: rows are genes and columns are samples
  ### grp: group information of the samples
  ### num: the number of top genes to be used based on variance (-1 [default]: use all the genes)
  ### isScaled: whether the given matrix is scaled or not
  ### sampleName: whether to present sample names or not
  ### component: to draw a plot with PC1 & PC2 or PC2 & PC3
  ### title: title of the plot
  ### legendTitle: title of the legend
  ### suppliment: if it is TRUE, this function additionally generates figures and tables
  ###             related to contributions. For example, which genes are highly contributed to
  ###             the PC1, PC2, etc.
  ### outDir: output directory for the plot
  pca_plot <- function(mat,
                       grp,
                       num = -1,
                       isScaled = FALSE,
                       sampleName = FALSE,
                       component=c("PC1&PC2", "PC2&PC3"),
                       title="PCA_Plot",
                       legendTitle="group",
                       suppliment=FALSE,
                       outDir="./") {
    ### load library
    if(!require(ggfortify, quietly = TRUE)) {
      install.packages("ggfortify")
      library(ggfortify, quietly = TRUE)
    }
    if(!require(FactoMineR, quietly = TRUE)) {
      install.packages("FactoMineR")
      library(FactoMineR, quietly = TRUE)
    }
    if(!require(factoextra, quietly = TRUE)) {
      install.packages("factoextra")
      library(factoextra, quietly = TRUE)
    }
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### select the top genes based on variance
    if(num >= 0 && num <= nrow(mat)) {
      v <- apply(mat, 1, var)
      v <- v[order(-v)]
      top_genes <- names(v)[1:num]
    } else {
      top_genes <- rownames(mat)
    }
    
    ### PCA
    grp <- as.character(grp)
    if(isScaled) {
      pca_result <- PCA(t(mat[top_genes,]), graph = FALSE, scale.unit = FALSE) 
    } else {
      pca_result <- PCA(t(mat[top_genes,]), graph = FALSE, scale.unit = TRUE)
    }
    colnames(pca_result$ind$coord) <- paste0("PC", 1:ncol(pca_result$ind$coord))
    colnames(pca_result$var$contrib) <- paste0("PC", 1:ncol(pca_result$var$contrib))
    pca_group <- data.frame(pca_result$ind$coord, group=grp)
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as png
    if(component[1] == "PC1&PC2") {
      if(sampleName) {
        ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
          labs(title=paste0(title, "_PC1-2")) +
          geom_text(aes(label=colnames(mat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16) +
          guides(color=guide_legend(legendTitle))
      } else {
        ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
          labs(title=paste0(title, "_PC1-2")) +
          geom_point() +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16) +
          guides(color=guide_legend(legendTitle))
      }
      ggsave(filename = paste0(outDir, title, "_PC1-2", ".png"), width = 10, height = 8)
      
      if(suppliment) {
        fviz_contrib(pca_result, choice = "var", axes = 1, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC1_contribution.png"), width = 12, height = 8)
        
        fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
      }
    } else if(component[1] == "PC2&PC3") {
      if(sampleName) {
        ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
          labs(title=paste0(title, "_PC2-3")) +
          geom_text(aes(label=colnames(mat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16) +
          guides(color=guide_legend(legendTitle))
      } else {
        ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
          labs(title=paste0(title, "_PC2-3")) +
          geom_point() +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16) +
          guides(color=guide_legend(legendTitle))
      }
      ggsave(filename = paste0(outDir, title, "_PC2-3", ".png"), width = 10, height = 8)
      
      if(suppliment) {
        fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
        
        fviz_contrib(pca_result, choice = "var", axes = 3, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC3_contribution.png"), width = 12, height = 8)
      }
    } else {
      stop("\"component\" parameter should be \"PC1&PC2\" or \"PC2&PC3\"")
    }
    
    if(suppliment) {
      write.xlsx2(data.frame(Gene_Symbol=rownames(pca_result$var$contrib), pca_result$var$contrib,
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(outDir, title, "_PC_contribution.xlsx"),
                  sheetName = "PCA_contribution", row.names = FALSE)
    }
  }
  
  ### make PCA plots
  for(factor in setdiff(c(factor_list, "Study.Day"), "Age.at.Enrollment")) {
    
    ### create result directory
    outDir <- paste0(output_dir, "2_PCA/", factor, "/")
    dir.create(path = outDir, showWarnings = FALSE, recursive = TRUE)
    
    ### PCA with all the cytokine levels
    pca_plot(t(plot_df_nw[,c(th1_cytokines, th2_cytokines)]), grp = plot_df_nw[,factor],
             title = paste0("PCA_Plot_NW_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(plot_df_ps[,c(th1_cytokines, th2_cytokines)]), grp = plot_df_ps[,factor],
             title = paste0("PCA_Plot_PS_", factor),
             legendTitle = factor,
             outDir = outDir)
    
    ### PCA with TH1 cytokines only
    pca_plot(t(plot_df_nw[,th1_cytokines]), grp = plot_df_nw[,factor],
             title = paste0("PCA_Plot_NW_TH1_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(plot_df_ps[,th1_cytokines]), grp = plot_df_ps[,factor],
             title = paste0("PCA_Plot_PS_TH1_", factor),
             legendTitle = factor,
             outDir = outDir)
    
    ### PCA with TH2 cytokines only
    pca_plot(t(plot_df_nw[,th2_cytokines]), grp = plot_df_nw[,factor],
             title = paste0("PCA_Plot_NW_TH2_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(plot_df_ps[,th2_cytokines]), grp = plot_df_ps[,factor],
             title = paste0("PCA_Plot_PS_TH2_", factor),
             legendTitle = factor,
             outDir = outDir)
    
  }
  
  
  ### 3. a) Line graph with all the samples
  
  
    
}
