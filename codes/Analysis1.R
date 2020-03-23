###
#   File name : Analysis1.R
#   Author    : Hyunjin Kim
#   Date      : Mar 23, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : First analysis of the Mohamed's flu09 data.
#               It has 5 tasks to do:
#               1. Explore the data set. Find any association between cytokine levels
#                  and various clinical factors.
#               2. PCA. Mark with Asthma info and with time point info to see
#                  if the samples are clustered well.
#               3. Beeswarm plots. Make the plot per cytokine and put all the samples
#                  to compare Asthma/IV and IV. Different time points are painted with
#                  different colored dots.
#               4. a) Line graph per each cytokine. The x-axis is time points and the y-axis is
#                     cytokine level. Color the graphs based on Asthma/IV vs IV, so that we can
#                     distinguish them in the plot.
#                  b) Line graph per each cytokine. Take mean of each time point and compare
#                     Asthma/IV vs IV. There would be only two lines in each graph.
#               5. Statistics table. For each time point, calculate mean/median difference
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
  
  
    
}
