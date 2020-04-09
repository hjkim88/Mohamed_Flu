###
#   File name : Analysis2.R
#   Author    : Hyunjin Kim
#   Date      : Apr 6, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Second analysis of the Mohamed's flu09 data.
#               This analysis is a follow-up analysis after the "Analysis1.R".
#               1. PCA. When performing PCA, try with log-transformed data and
#                  turn on/off the embedded scaling in the PCA function.
#               2. Line graph per each cytokine. But this time, draw unique color dots
#                  for the significantly different time points and specify p-values in it.
#                  And also add standard errors.
#               3. Correlation between cytokine levels and age.
#                  Continuous and discrete (>=50 vs <50)
#
#   Instruction
#               1. Source("Analysis2.R")
#               2. Run the function "flu09_analysis2" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Analysis2.R/Analysis2.R")
#               > flu09_analysis2(data_path="./data/flu09_cytokine.rda",
#                                 output_dir="./results/analysis2/")
###

flu09_analysis2 <- function(data_path="./data/flu09_cytokine.rda",
                            output_dir="./results/analysis2/") {
  
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
  
  
  ### 1. PCA
  
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
  
  ### log-transform the data
  log_df_nw <- plot_df_nw
  log_df_nw[,c(th1_cytokines, th2_cytokines)] <- log(log_df_nw[,c(th1_cytokines, th2_cytokines)]+1)
  log_df_ps <- plot_df_ps
  log_df_ps[,c(th1_cytokines, th2_cytokines)] <- log(log_df_ps[,c(th1_cytokines, th2_cytokines)]+1)
  
  ### make PCA plots
  for(factor in setdiff(c(factor_list, "Study.Day", "IV.Resp", "Age.over.50"), "Age.at.Enrollment")) {
    
    ### A. with log-transformed and without embedded scaling 
    
    ### create result directory
    outDir <- paste0(output_dir, "1_PCA/", factor, "/log_transformed_and_no_scaling/")
    dir.create(path = outDir, showWarnings = FALSE, recursive = TRUE)
    
    ### PCA with all the cytokine levels
    pca_plot(t(log_df_nw[,c(th1_cytokines, th2_cytokines)]), grp = log_df_nw[,factor],
             isScaled = TRUE,
             title = paste0("PCA_Plot_log_transformed_NW_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(log_df_ps[,c(th1_cytokines, th2_cytokines)]), grp = log_df_ps[,factor],
             isScaled = TRUE,
             title = paste0("PCA_Plot_log_transformed_PS_", factor),
             legendTitle = factor,
             outDir = outDir)
    
    ### PCA with TH1 cytokines only
    pca_plot(t(log_df_nw[,th1_cytokines]), grp = log_df_nw[,factor],
             isScaled = TRUE,
             title = paste0("PCA_Plot_log_transformed_NW_TH1_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(log_df_ps[,th1_cytokines]), grp = log_df_ps[,factor],
             isScaled = TRUE,
             title = paste0("PCA_Plot_log_transformed_PS_TH1_", factor),
             legendTitle = factor,
             outDir = outDir)
    
    ### PCA with TH2 cytokines only
    pca_plot(t(log_df_nw[,th2_cytokines]), grp = log_df_nw[,factor],
             isScaled = TRUE,
             title = paste0("PCA_Plot_log_transformed_NW_TH2_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(log_df_ps[,th2_cytokines]), grp = log_df_ps[,factor],
             isScaled = TRUE,
             title = paste0("PCA_Plot_log_transformed_PS_TH2_", factor),
             legendTitle = factor,
             outDir = outDir)
    
    ### B. with log-transformed and with embedded scaling
    
    ### create result directory
    outDir <- paste0(output_dir, "1_PCA/", factor, "/log_transformed_and_scaling/")
    dir.create(path = outDir, showWarnings = FALSE, recursive = TRUE)
    
    ### PCA with all the cytokine levels
    pca_plot(t(log_df_nw[,c(th1_cytokines, th2_cytokines)]), grp = log_df_nw[,factor],
             title = paste0("PCA_Plot_log_transformed_NW_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(log_df_ps[,c(th1_cytokines, th2_cytokines)]), grp = log_df_ps[,factor],
             title = paste0("PCA_Plot_log_transformed_PS_", factor),
             legendTitle = factor,
             outDir = outDir)
    
    ### PCA with TH1 cytokines only
    pca_plot(t(log_df_nw[,th1_cytokines]), grp = log_df_nw[,factor],
             title = paste0("PCA_Plot_log_transformed_NW_TH1_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(log_df_ps[,th1_cytokines]), grp = log_df_ps[,factor],
             title = paste0("PCA_Plot_log_transformed_PS_TH1_", factor),
             legendTitle = factor,
             outDir = outDir)
    
    ### PCA with TH2 cytokines only
    pca_plot(t(log_df_nw[,th2_cytokines]), grp = log_df_nw[,factor],
             title = paste0("PCA_Plot_log_transformed_NW_TH2_", factor),
             legendTitle = factor,
             outDir = outDir)
    pca_plot(t(log_df_ps[,th2_cytokines]), grp = log_df_ps[,factor],
             title = paste0("PCA_Plot_log_transformed_PS_TH2_", factor),
             legendTitle = factor,
             outDir = outDir)
    
  }
  
  
  ### 2. Line graph per each cytokine + p-value + standard error
  
  ### option for 95% confidence level
  is_conf_lvl <- FALSE
  
  ### make the line graphs
  factors <- setdiff(c(factor_list, "IV.Resp", "Age.over.50"), "Age.at.Enrollment")
  plot_df_nw_mean <- vector("list", length = length(factors))
  names(plot_df_nw_mean) <- factors
  plot_df_nw_sem <- vector("list", length = length(factors))
  names(plot_df_nw_sem) <- factors
  plot_df_ps_mean <- vector("list", length = length(factors))
  names(plot_df_ps_mean) <- factors
  plot_df_ps_sem <- vector("list", length = length(factors))
  names(plot_df_ps_sem) <- factors
  stat_table_nw <- vector("list", length = length(factors))
  names(stat_table_nw) <- factors
  stat_table_ps <- vector("list", length = length(factors))
  names(stat_table_ps) <- factors
  for(factor in factors) {
    
    ### create result directory
    outDir <- paste0(output_dir, "2_Line_Graph/", factor, "/")
    dir.create(path = outDir, showWarnings = FALSE, recursive = TRUE)
    
    ### create mean data
    ### set the threshold for "Study.Day" selection
    ### if there are less than X samples in the category, we will not include that in the analysis
    sampleNum_threshold <- 5
    ### NW
    v1 <- NULL
    v2 <- NULL
    for(x in unique(plot_df_nw[,factor])) {
      y <- unique(plot_df_nw[which(plot_df_nw[,factor] == x),"Study.Day"])
      for(z in y) {
        if(length(intersect(which(plot_df_nw[,factor] == x),
                            which(plot_df_nw[,"Study.Day"] == z))) >= sampleNum_threshold) {
          v1 <- c(v1, x)
          v2 <- c(v2, z)
        }
      }
    }
    row_names <- paste(v1, v2, sep = "_")
    plot_df_nw_mean[[factor]] <- data.frame(matrix(0, length(row_names), length(c(th1_cytokines, th2_cytokines))))
    rownames(plot_df_nw_mean[[factor]]) <- row_names
    colnames(plot_df_nw_mean[[factor]]) <- c(th1_cytokines, th2_cytokines)
    plot_df_nw_sem[[factor]] <- data.frame(matrix(0, length(row_names), length(c(th1_cytokines, th2_cytokines))))
    rownames(plot_df_nw_sem[[factor]]) <- row_names
    colnames(plot_df_nw_sem[[factor]]) <- c(th1_cytokines, th2_cytokines)
    for(x in unique(plot_df_nw[,factor])) {
      y <- unique(plot_df_nw[which(plot_df_nw[,factor] == x),"Study.Day"])
      for(z in y) {
        w <- intersect(which(plot_df_nw[,factor] == x),
                       which(plot_df_nw[,"Study.Day"] == z))
        if(length(w) >= sampleNum_threshold) {
          plot_df_nw_mean[[factor]][paste0(x, "_", z),] <- apply(plot_df_nw[w,
                                                                            c(th1_cytokines, th2_cytokines)], 2, mean, na.rm=TRUE)
          plot_df_nw_sem[[factor]][paste0(x, "_", z),] <- apply(plot_df_nw[w,
                                                                           c(th1_cytokines, th2_cytokines)], 2, function(a) {
                                                                             b <- a[!is.na(a)]
                                                                             return(sd(b)/sqrt(length(b)))
                                                                           })
        }
      }
    }
    plot_df_nw_mean[[factor]] <- data.frame(plot_df_nw_mean[[factor]],v1, v2, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(plot_df_nw_mean[[factor]])[(ncol(plot_df_nw_mean[[factor]])-1):ncol(plot_df_nw_mean[[factor]])] <- c(factor, "Study.Day")
    colnames(plot_df_nw_sem[[factor]]) <- paste0(colnames(plot_df_nw_sem[[factor]]), "_sem")
    ### Plasma
    v1 <- NULL
    v2 <- NULL
    for(x in unique(plot_df_ps[,factor])) {
      y <- unique(plot_df_ps[which(plot_df_ps[,factor] == x),"Study.Day"])
      for(z in y) {
        if(length(intersect(which(plot_df_ps[,factor] == x),
                            which(plot_df_ps[,"Study.Day"] == z))) >= sampleNum_threshold) {
          v1 <- c(v1, x)
          v2 <- c(v2, z)
        }
      }
    }
    row_names <- paste(v1, v2, sep = "_")
    plot_df_ps_mean[[factor]] <- data.frame(matrix(0, length(row_names), length(c(th1_cytokines, th2_cytokines))))
    rownames(plot_df_ps_mean[[factor]]) <- row_names
    colnames(plot_df_ps_mean[[factor]]) <- c(th1_cytokines, th2_cytokines)
    plot_df_ps_sem[[factor]] <- data.frame(matrix(0, length(row_names), length(c(th1_cytokines, th2_cytokines))))
    rownames(plot_df_ps_sem[[factor]]) <- row_names
    colnames(plot_df_ps_sem[[factor]]) <- c(th1_cytokines, th2_cytokines)
    for(x in unique(plot_df_ps[,factor])) {
      y <- unique(plot_df_ps[which(plot_df_ps[,factor] == x),"Study.Day"])
      for(z in y) {
        w <- intersect(which(plot_df_ps[,factor] == x),
                       which(plot_df_ps[,"Study.Day"] == z))
        if(length(w) >= sampleNum_threshold) {
          plot_df_ps_mean[[factor]][paste0(x, "_", z),] <- apply(plot_df_ps[w,
                                                                            c(th1_cytokines, th2_cytokines)], 2, mean, na.rm=TRUE)
          plot_df_ps_sem[[factor]][paste0(x, "_", z),] <- apply(plot_df_ps[w,
                                                                           c(th1_cytokines, th2_cytokines)], 2, function(a) {
                                                                             b <- a[!is.na(a)]
                                                                             return(sd(b)/sqrt(length(b)))
                                                                           })
        }
      }
    }
    plot_df_ps_mean[[factor]] <- data.frame(plot_df_ps_mean[[factor]],v1, v2, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(plot_df_ps_mean[[factor]])[(ncol(plot_df_ps_mean[[factor]])-1):ncol(plot_df_ps_mean[[factor]])] <- c(factor, "Study.Day")
    colnames(plot_df_ps_sem[[factor]]) <- paste0(colnames(plot_df_ps_sem[[factor]]), "_sem")
    
    ### fill out the statistics table
    ### This is needed because we would like to specify t-test and ANOVA p-value in the line graphs
    ### NW
    
    ### make all possible pairs in the factor group
    unique_vals <- unique(plot_df_nw_mean[[factor]][,factor])
    factor_pairs <- matrix("", choose(length(unique_vals), 2), 2)
    cnt <- 1
    for(i in 1:(length(unique_vals)-1)) {
      for(j in (i+1):length(unique_vals)) {
        factor_pairs[cnt,1] <- unique_vals[i]
        factor_pairs[cnt,2] <- unique_vals[j]
        cnt <- cnt + 1
      }
    }
    
    ### for each pair
    stat_table_nw[[factor]] <- vector("list", length = nrow(factor_pairs))
    names(stat_table_nw[[factor]]) <- paste0(factor_pairs[,1], "_vs_", factor_pairs[,2])
    for(i in 1:nrow(factor_pairs)) {
      
      ### create an empty matrix
      specific_study_days <- intersect(unique(plot_df_nw_mean[[factor]][which(plot_df_nw_mean[[factor]][,factor] == factor_pairs[i,1]),"Study.Day"]),
                                       unique(plot_df_nw_mean[[factor]][which(plot_df_nw_mean[[factor]][,factor] == factor_pairs[i,2]),"Study.Day"]))
      specific_study_days <- as.character(specific_study_days[order(specific_study_days)])
      mat <- data.frame(matrix(NA, length(c(th1_cytokines, th2_cytokines)), length(specific_study_days)),
                        stringsAsFactors = FALSE, check.names = FALSE)
      rownames(mat) <- c(th1_cytokines, th2_cytokines)
      colnames(mat) <- specific_study_days
      
      for(cytokine in c(th1_cytokines, th2_cytokines)) {
        for(study_day in specific_study_days) {
          ### ANOVA
          a <- plot_df_nw[intersect(which(plot_df_nw[,factor] == factor_pairs[i,1]),
                                    which(plot_df_nw[,"Study.Day"] == as.numeric(study_day))),cytokine]
          b <- plot_df_nw[intersect(which(plot_df_nw[,factor] == factor_pairs[i,2]),
                                    which(plot_df_nw[,"Study.Day"] == as.numeric(study_day))),cytokine]
          c <- data.frame(cytokine_level=c(a, b), group=c(rep("a", length(a)), rep("b", length(b))),
                          stringsAsFactors = FALSE, check.names = FALSE)
          mat[cytokine,study_day] <- anova(lm(cytokine_level~group,data=c))$"Pr(>F)"[1]
        }
      }
      
      ### save the matrix
      stat_table_nw[[factor]][[i]] <- mat
      
    }
    
    ### Plasma
    
    ### make all possible pairs in the factor group
    unique_vals <- unique(plot_df_ps_mean[[factor]][,factor])
    factor_pairs <- matrix("", choose(length(unique_vals), 2), 2)
    cnt <- 1
    for(i in 1:(length(unique_vals)-1)) {
      for(j in (i+1):length(unique_vals)) {
        factor_pairs[cnt,1] <- unique_vals[i]
        factor_pairs[cnt,2] <- unique_vals[j]
        cnt <- cnt + 1
      }
    }
    
    ### for each pair
    stat_table_ps[[factor]] <- vector("list", length = nrow(factor_pairs))
    names(stat_table_ps[[factor]]) <- paste0(factor_pairs[,1], "_vs_", factor_pairs[,2])
    for(i in 1:nrow(factor_pairs)) {
      
      ### create an empty matrix
      specific_study_days <- intersect(unique(plot_df_ps_mean[[factor]][which(plot_df_ps_mean[[factor]][,factor] == factor_pairs[i,1]),"Study.Day"]),
                                       unique(plot_df_ps_mean[[factor]][which(plot_df_ps_mean[[factor]][,factor] == factor_pairs[i,2]),"Study.Day"]))
      specific_study_days <- as.character(specific_study_days[order(specific_study_days)])
      mat <- data.frame(matrix(NA, length(c(th1_cytokines, th2_cytokines)), length(specific_study_days)),
                        stringsAsFactors = FALSE, check.names = FALSE)
      rownames(mat) <- c(th1_cytokines, th2_cytokines)
      colnames(mat) <- specific_study_days
      
      for(cytokine in c(th1_cytokines, th2_cytokines)) {
        for(study_day in specific_study_days) {
          ### ANOVA
          a <- plot_df_ps[intersect(which(plot_df_ps[,factor] == factor_pairs[i,1]),
                                    which(plot_df_ps[,"Study.Day"] == as.numeric(study_day))),cytokine]
          b <- plot_df_ps[intersect(which(plot_df_ps[,factor] == factor_pairs[i,2]),
                                    which(plot_df_ps[,"Study.Day"] == as.numeric(study_day))),cytokine]
          c <- data.frame(cytokine_level=c(a, b), group=c(rep("a", length(a)), rep("b", length(b))),
                          stringsAsFactors = FALSE, check.names = FALSE)
          mat[cytokine,study_day] <- anova(lm(cytokine_level~group,data=c))$"Pr(>F)"[1]
        }
      }
      
      ### save the matrix
      stat_table_ps[[factor]][[i]] <- mat
      
    }
    
    ### add standard error
    plot_df_nw_mean[[factor]] <- cbind(plot_df_nw_mean[[factor]], plot_df_nw_sem[[factor]])
    plot_df_ps_mean[[factor]] <- cbind(plot_df_ps_mean[[factor]], plot_df_ps_sem[[factor]])
    
    
    ### Start drawing the plots!
    ### NW
    
    ### TH1 cytokines only
    p <- vector("list", length = length(th1_cytokines))
    names(p) <- th1_cytokines
    
    ### line graph per cytokine - mean
    for(cytokine in th1_cytokines) {
      top <- ""
      p[[cytokine]] <- ggplot(plot_df_nw_mean[[factor]], aes_string(x="Study.Day", y=cytokine, group=factor)) +
        geom_point(aes_string(color=factor), size=4, alpha=0.5) +
        geom_line(aes_string(color=factor), size=2, alpha=0.5) +
        theme_classic(base_size = 16)
      
      ### if there are comparisons that have significant p-values show them in the plot
      for(comparison in names(stat_table_nw[[factor]])) {
        for(col in colnames(stat_table_nw[[factor]][[comparison]])) {
          pv <- as.numeric(stat_table_nw[[factor]][[comparison]][cytokine,col])
          if(pv < 0.05) {
            temp <- strsplit(comparison, split = "_vs_", fixed = TRUE)[[1]]
            temp <- paste0(temp, "_", col)
            temp <- plot_df_nw_mean[[factor]][temp,cytokine]
            ### black circles
            p[[cytokine]] <- p[[cytokine]] +
              geom_point(x = as.numeric(col), y = temp[1], col = "black", size = 4, shape = 1) +
              geom_point(x = as.numeric(col), y = temp[2], col = "black", size = 4, shape = 1)
            ### add p-value to the top
            top <- paste(top, comparison, "|", col, "Day |", signif(pv, digits = 2), "\n")
          }
        }
      }
      
      ### attach top label to the ggplot
      p[[cytokine]] <- p[[cytokine]] +
        ggtitle(label = top) +
        theme(plot.title = element_text(size = 10, hjust=0.5))
      
      ### if 95% confidence level is TRUE,
      if(is_conf_lvl) {
        p[[cytokine]] <- p[[cytokine]] +
          geom_errorbar(aes_string(ymin=paste0(cytokine, "-1.96*", cytokine, "_sem"),
                                   ymax=paste0(cytokine, "+1.96*", cytokine, "_sem"),
                                   color=factor),
                        width=1,
                        alpha=1,
                        position=position_dodge(0.3))
      }
    }
    
    ### arrange the plots and print out
    fName <- paste0("Line_Plot_NW_TH1_Mean_", factor)
    g <- arrangeGrob(grobs = p,
                     layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                     top = fName)
    ggsave(file = paste0(outDir, fName, ".png"), g, width = 24, height = 12)
    
    ### TH2 cytokines only
    p <- vector("list", length = length(th2_cytokines))
    names(p) <- th2_cytokines
    
    ### line graph per cytokine - mean
    for(cytokine in th2_cytokines) {
      top <- ""
      p[[cytokine]] <- ggplot(plot_df_nw_mean[[factor]], aes_string(x="Study.Day", y=cytokine, group=factor)) +
        geom_point(aes_string(color=factor), size=4, alpha=0.5) +
        geom_line(aes_string(color=factor), size=2, alpha=0.5) +
        theme_classic(base_size = 16)
      
      ### if there are comparisons that have significant p-values show them in the plot
      for(comparison in names(stat_table_nw[[factor]])) {
        for(col in colnames(stat_table_nw[[factor]][[comparison]])) {
          pv <- as.numeric(stat_table_nw[[factor]][[comparison]][cytokine,col])
          if(pv < 0.05) {
            temp <- strsplit(comparison, split = "_vs_", fixed = TRUE)[[1]]
            temp <- paste0(temp, "_", col)
            temp <- plot_df_nw_mean[[factor]][temp,cytokine]
            ### black circles
            p[[cytokine]] <- p[[cytokine]] +
              geom_point(x = as.numeric(col), y = temp[1], col = "black", size = 4, shape = 1) +
              geom_point(x = as.numeric(col), y = temp[2], col = "black", size = 4, shape = 1)
            ### add p-value to the top
            top <- paste(top, comparison, "|", col, "Day |", signif(pv, digits = 2), "\n")
          }
        }
      }
      
      ### attach top label to the ggplot
      p[[cytokine]] <- p[[cytokine]] +
        ggtitle(label = top) +
        theme(plot.title = element_text(size = 10, hjust=0.5))
      
      ### if 95% confidence level is TRUE,
      if(is_conf_lvl) {
        p[[cytokine]] <- p[[cytokine]] +
          geom_errorbar(aes_string(ymin=paste0(cytokine, "-1.96*", cytokine, "_sem"),
                                   ymax=paste0(cytokine, "+1.96*", cytokine, "_sem"),
                                   color=factor),
                        width=1,
                        alpha=1,
                        position=position_dodge(0.3))
      }
    }
    
    ### arrange the plots and print out
    fName <- paste0("Line_Plot_NW_TH2_Mean_", factor)
    g <- arrangeGrob(grobs = p,
                     layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)),
                     top = fName)
    ggsave(file = paste0(outDir, fName, ".png"), g, width = 24, height = 12)
    
    
    ### Plasma
    
    ### TH1 cytokines only
    p <- vector("list", length = length(th1_cytokines))
    names(p) <- th1_cytokines
    
    ### line graph per cytokine - mean
    for(cytokine in th1_cytokines) {
      top <- ""
      p[[cytokine]] <- ggplot(plot_df_ps_mean[[factor]], aes_string(x="Study.Day", y=cytokine, group=factor)) +
        geom_point(aes_string(color=factor), size=4, alpha=0.5) +
        geom_line(aes_string(color=factor), size=2, alpha=0.5) +
        theme_classic(base_size = 16)
      
      ### if there are comparisons that have significant p-values show them in the plot
      for(comparison in names(stat_table_ps[[factor]])) {
        for(col in colnames(stat_table_ps[[factor]][[comparison]])) {
          pv <- as.numeric(stat_table_ps[[factor]][[comparison]][cytokine,col])
          if(pv < 0.05) {
            temp <- strsplit(comparison, split = "_vs_", fixed = TRUE)[[1]]
            temp <- paste0(temp, "_", col)
            temp <- plot_df_ps_mean[[factor]][temp,cytokine]
            ### black circles
            p[[cytokine]] <- p[[cytokine]] +
              geom_point(x = as.numeric(col), y = temp[1], col = "black", size = 4, shape = 1) +
              geom_point(x = as.numeric(col), y = temp[2], col = "black", size = 4, shape = 1)
            ### add p-value to the top
            top <- paste(top, comparison, "|", col, "Day |", signif(pv, digits = 2), "\n")
          }
        }
      }
      
      ### attach top label to the ggplot
      p[[cytokine]] <- p[[cytokine]] +
        ggtitle(label = top) +
        theme(plot.title = element_text(size = 10, hjust=0.5))
      
      ### if 95% confidence level is TRUE,
      if(is_conf_lvl) {
        p[[cytokine]] <- p[[cytokine]] +
          geom_errorbar(aes_string(ymin=paste0(cytokine, "-1.96*", cytokine, "_sem"),
                                   ymax=paste0(cytokine, "+1.96*", cytokine, "_sem"),
                                   color=factor),
                        width=1,
                        alpha=1,
                        position=position_dodge(0.3))
      }
    }
    
    ### arrange the plots and print out
    fName <- paste0("Line_Plot_PS_TH1_Mean_", factor)
    g <- arrangeGrob(grobs = p,
                     layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                     top = fName)
    ggsave(file = paste0(outDir, fName, ".png"), g, width = 24, height = 12)
    
    ### TH2 cytokines only
    p <- vector("list", length = length(th2_cytokines))
    names(p) <- th2_cytokines
    
    ### line graph per cytokine - mean
    for(cytokine in th2_cytokines) {
      top <- ""
      p[[cytokine]] <- ggplot(plot_df_ps_mean[[factor]], aes_string(x="Study.Day", y=cytokine, group=factor)) +
        geom_point(aes_string(color=factor), size=4, alpha=0.5) +
        geom_line(aes_string(color=factor), size=2, alpha=0.5) +
        theme_classic(base_size = 16)
      
      ### if there are comparisons that have significant p-values show them in the plot
      for(comparison in names(stat_table_ps[[factor]])) {
        for(col in colnames(stat_table_ps[[factor]][[comparison]])) {
          pv <- as.numeric(stat_table_ps[[factor]][[comparison]][cytokine,col])
          if(pv < 0.05) {
            temp <- strsplit(comparison, split = "_vs_", fixed = TRUE)[[1]]
            temp <- paste0(temp, "_", col)
            temp <- plot_df_ps_mean[[factor]][temp,cytokine]
            ### black circles
            p[[cytokine]] <- p[[cytokine]] +
              geom_point(x = as.numeric(col), y = temp[1], col = "black", size = 4, shape = 1) +
              geom_point(x = as.numeric(col), y = temp[2], col = "black", size = 4, shape = 1)
            ### add p-value to the top
            top <- paste(top, comparison, "|", col, "Day |", signif(pv, digits = 2), "\n")
          }
        }
      }
      
      ### attach top label to the ggplot
      p[[cytokine]] <- p[[cytokine]] +
        ggtitle(label = top) +
        theme(plot.title = element_text(size = 10, hjust=0.5))
      
      ### if 95% confidence level is TRUE,
      if(is_conf_lvl) {
        p[[cytokine]] <- p[[cytokine]] +
          geom_errorbar(aes_string(ymin=paste0(cytokine, "-1.96*", cytokine, "_sem"),
                                   ymax=paste0(cytokine, "+1.96*", cytokine, "_sem"),
                                   color=factor),
                        width=1,
                        alpha=1,
                        position=position_dodge(0.3))
      }
    }
    
    ### arrange the plots and print out
    fName <- paste0("Line_Plot_PS_TH2_Mean_", factor)
    g <- arrangeGrob(grobs = p,
                     layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)),
                     top = fName)
    ggsave(file = paste0(outDir, fName, ".png"), g, width = 24, height = 12)
    
  }
  
  
  ### 3. Comparing cytokine levels between Age >= 50 vs Age < 50
  
  ### create result directory
  outDir <- paste0(output_dir, "3_Age_Cytokine_level/")
  dir.create(path = outDir, showWarnings = FALSE, recursive = TRUE)
  
  ### Function to produce summary statistics (mean and +/- sd)
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  ### NW
  ### TH1 cytokines plot
  p <- vector("list", length = length(th1_cytokines))
  names(p) <- th1_cytokines
  
  ### box plot per cytokine
  for(cytokine in th1_cytokines) {
    means <- aggregate(as.formula(paste0(cytokine, "~", "Age.over.50")), plot_df_nw, mean)
    min_v <- min(plot_df_nw[,cytokine], na.rm = TRUE)
    max_v <- max(plot_df_nw[,cytokine], na.rm = TRUE)
    means[,cytokine] <- paste0("Mean: ", signif(as.numeric(means[,cytokine]), digits = 5))
    lbl <- sapply(levels(as.factor(plot_df_nw$Age.over.50)), function(x) {
      return(length(intersect(which(plot_df_nw[,"Age.over.50"] == x),
                              which(!is.na(plot_df_nw[,cytokine])))))
    })
    lbl <- paste0(levels(as.factor(plot_df_nw$Age.over.50)),
                  paste0("(n=", lbl, ")"))
    p[[cytokine]] <- ggplot(plot_df_nw, aes_string(x="Age.over.50", y=cytokine)) +
      geom_boxplot(na.rm = TRUE) +
      stat_compare_means(na.rm = TRUE) +
      geom_text(data = means, aes_string(label = cytokine, y = min_v-(max_v*0.05))) +
      scale_x_discrete(labels = lbl) +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(size=10),
            legend.text = element_text(size=8))
  }
  
  ### arrange the plots and print out
  fName <- paste0("Box_Plot_NW_TH1_", "Age.over.50")
  g <- arrangeGrob(grobs = p,
                   layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                   top = fName)
  ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
  
  ### TH2 cytokines plot
  p <- vector("list", length = length(th2_cytokines))
  names(p) <- th2_cytokines
  
  ### box plot per cytokine
  for(cytokine in th2_cytokines) {
    means <- aggregate(as.formula(paste0(cytokine, "~", "Age.over.50")), plot_df_nw, mean)
    min_v <- min(plot_df_nw[,cytokine], na.rm = TRUE)
    max_v <- max(plot_df_nw[,cytokine], na.rm = TRUE)
    means[,cytokine] <- paste0("Mean: ", signif(as.numeric(means[,cytokine]), digits = 5))
    lbl <- sapply(levels(as.factor(plot_df_nw$Age.over.50)), function(x) {
      return(length(intersect(which(plot_df_nw[,"Age.over.50"] == x),
                              which(!is.na(plot_df_nw[,cytokine])))))
    })
    lbl <- paste0(levels(as.factor(plot_df_nw$Age.over.50)),
                  paste0("(n=", lbl, ")"))
    p[[cytokine]] <- ggplot(plot_df_nw, aes_string(x="Age.over.50", y=cytokine)) +
      geom_boxplot(na.rm = TRUE) +
      stat_compare_means(na.rm = TRUE) +
      geom_text(data = means, aes_string(label = cytokine, y = min_v-(max_v*0.05))) +
      scale_x_discrete(labels = lbl) +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(size=10),
            legend.text = element_text(size=8))
  }
  
  ### arrange the plots and print out
  fName <- paste0("Box_Plot_NW_TH2_", "Age.over.50")
  g <- arrangeGrob(grobs = p,
                   layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)),
                   top = fName)
  ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
  
  ### Plasma
  ### TH1 cytokines plot
  p <- vector("list", length = length(th1_cytokines))
  names(p) <- th1_cytokines
  
  ### box plot per cytokine
  for(cytokine in th1_cytokines) {
    means <- aggregate(as.formula(paste0(cytokine, "~", "Age.over.50")), plot_df_ps, mean)
    min_v <- min(plot_df_ps[,cytokine], na.rm = TRUE)
    max_v <- max(plot_df_ps[,cytokine], na.rm = TRUE)
    means[,cytokine] <- paste0("Mean: ", signif(as.numeric(means[,cytokine]), digits = 5))
    lbl <- sapply(levels(as.factor(plot_df_ps$Age.over.50)), function(x) {
      return(length(intersect(which(plot_df_ps[,"Age.over.50"] == x),
                              which(!is.na(plot_df_ps[,cytokine])))))
    })
    lbl <- paste0(levels(as.factor(plot_df_ps$Age.over.50)),
                  paste0("(n=", lbl, ")"))
    p[[cytokine]] <- ggplot(plot_df_ps, aes_string(x="Age.over.50", y=cytokine)) +
      geom_boxplot(na.rm = TRUE) +
      stat_compare_means(na.rm = TRUE) +
      geom_text(data = means, aes_string(label = cytokine, y = min_v-(max_v*0.05))) +
      scale_x_discrete(labels = lbl) +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(size=10),
            legend.text = element_text(size=8))
  }
  
  ### arrange the plots and print out
  fName <- paste0("Box_Plot_PS_TH1_", "Age.over.50")
  g <- arrangeGrob(grobs = p,
                   layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)),
                   top = fName)
  ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
  
  ### TH2 cytokines plot
  p <- vector("list", length = length(th2_cytokines))
  names(p) <- th2_cytokines
  
  ### box plot per cytokine
  for(cytokine in th2_cytokines) {
    means <- aggregate(as.formula(paste0(cytokine, "~", "Age.over.50")), plot_df_ps, mean)
    min_v <- min(plot_df_ps[,cytokine], na.rm = TRUE)
    max_v <- max(plot_df_ps[,cytokine], na.rm = TRUE)
    means[,cytokine] <- paste0("Mean: ", signif(as.numeric(means[,cytokine]), digits = 5))
    lbl <- sapply(levels(as.factor(plot_df_ps$Age.over.50)), function(x) {
      return(length(intersect(which(plot_df_ps[,"Age.over.50"] == x),
                              which(!is.na(plot_df_ps[,cytokine])))))
    })
    lbl <- paste0(levels(as.factor(plot_df_ps$Age.over.50)),
                  paste0("(n=", lbl, ")"))
    p[[cytokine]] <- ggplot(plot_df_ps, aes_string(x="Age.over.50", y=cytokine)) +
      geom_boxplot(na.rm = TRUE) +
      stat_compare_means(na.rm = TRUE) +
      geom_text(data = means, aes_string(label = cytokine, y = min_v-(max_v*0.05))) +
      scale_x_discrete(labels = lbl) +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(size=10),
            legend.text = element_text(size=8))
  }
  
  ### arrange the plots and print out
  fName <- paste0("Box_Plot_PS_TH2_", "Age.over.50")
  g <- arrangeGrob(grobs = p,
                   layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)),
                   top = fName)
  ggsave(file = paste0(outDir, fName, ".png"), g, width = 22, height = 12)
  
}
