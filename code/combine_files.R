#Script to pull in and combine/summarize metrics from the different GWAS methods
setClass("config",
         slots = c(ncores = "numeric",
                   cytokine = "character",
                   patient_metadata = "character",
                   group = "character",
                   tree = "character",
                   genome = "character",
                   ml_methods = "character"))

config <- new("config",
              ncores = c(4),
              cytokine = c("serum_IL.2Ra"),
              patient_metadata = c("data/patient_metadata.csv"),
              group = c("raw",
                        "adjusted"),
              tree = c("data/cytokine_rooted_tree.tree"),
              genome = c("core",
                         "pan"),
              ml_methods = c("glmnet",
                             "rf",
                             "svmRadial"))

setClass("snakemake", slots = c(input = "character",
                                output = "character",
                                wildcards = "character",
                                log = "character",
                                resources = "numeric",
                                params = "character",
                                ml_methods = "character"))
snakemake <- new("snakemake",
                 input = c(R = "code/prepro_overall.R",
                           patient_metadata = config@patient_metadata[1][[1]],
                           pheno = c(paste0("data/pheno/", config@cytokine[1][[1]], "/", config@group[1][[1]],".tsv"))),
                 wildcards = c(cytokine = config@cytokine[1][[1]],
                               genome = c("pan")),
                 log = c(paste0("log/", config@cytokine[1][[1]],"_prepro_overall.txt")),
                 resources = c(ncores = config@ncores[1][[1]]),
                 params = c(tree = config@tree[1][[1]]))

# source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
# doFuture::registerDoFuture()
# future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#RESULT PATH ----
result_path <- paste0("../snakemake-cytokine-analysis/results/",
                      snakemake@wildcards[['cytokine']])

#MIKROPML ----
data_proc_path <- list.files(paste0("../snakemake-cytokine-analysis/data/mikropml/",
                                    snakemake@wildcards[['cytokine']]),
                             pattern = "*.dat_proc.Rds")

print(paste0("Reading in preprocessed data frame from ", data_proc_path[1]))

mikropml_frame <-  readRDS(paste0("../snakemake-cytokine-analysis/data/mikropml/",
                                  snakemake@wildcards[['cytokine']], "/", data_proc_path[1]))

mikropml_names <- names(mikropml_frame[['grp_feats']])
mikropml_names <- gsub("_1$", "", mikropml_names)
mikropml_names <- gsub("`", "", mikropml_names)

perf_list <- list.files(path = result_path,
                        pattern = "performance_results.csv")
perf_list <- perf_list[grepl(snakemake@wildcards[['genome']], perf_list)]
perf_combine <- c()

for(i in 1:length(perf_list)){
  a <- read.csv(paste0(result_path, "/", perf_list[i]))
  file_name <- str_split(perf_list[i], pattern = "\\.", simplify = TRUE)
  alt <- a %>%
    mutate("group" = file_name[1],
           "genome" = file_name[2])
  perf_combine <- rbind(perf_combine,
                        alt)
}

bench_list <-list.files(path = result_path,
                        pattern = "benchmarks_results.csv")
bench_list <- bench_list[grepl(snakemake@wildcards[['genome']], bench_list)]
benchmark_combine <- c()

for (i in 1:length(bench_list)){
  a <- read.csv(paste0(result_path, "/", bench_list[i]))
  benchmark_combine <- rbind(benchmark_combine,
                             a)
}

feat_path <- paste0(result_path,
                    "/runs")
feat_import_list <- list.files(path = feat_path,
                               pattern = "feature-importance.csv")
feat_import_list <- feat_import_list[grepl(snakemake@wildcards[['genome']], feat_import_list)]
feat_import_combine <- c()

# for (i in 1:4){
for (i in 1:length(feat_import_list)){
  a <- read.csv(paste0(feat_path, "/", feat_import_list[i]))
  file_name <- str_split(feat_import_list[i], pattern = "\\.", simplify = TRUE)
  alt <- a %>%
    mutate("GWAS_method" = "mikropml",
           "group" = file_name[1],
           "genome" = file_name[2]) %>%
    rename("GWAS_test" = "method")
  
  feat_import_combine <- rbind(feat_import_combine,
                               alt)
}

feat_import_combine <- feat_import_combine %>%
  mutate(names = gsub("`_1", "", names),
         names = gsub("`", "", names),
         ) # %>%
  # group_by(names, GWAS_test, group, genome, GWAS_method, perf_metric_name) %>%
  # summarize(mean_perf_metric_diff = mean(perf_metric_diff),
  #           mean_perf_metric = mean(perf_metric),
  #           sd_perf_metric = sd(perf_metric))

rm(a)
rm(alt)
rm(file_name)

#summarize feature importance metrics
feat_importance_collapsed <- tibble("perf_metric" = c(),
                                    "perf_metric_diff" = c(),
                                    "names" = c(),
                                    "GWAS_test" = c(),
                                    "perf_metric_name" = c(),
                                    "seed" = c(),
                                    "GWAS_method" = c(),
                                    "group" = c(),
                                    "genome" = c())
table_out <- c()
string_clean <- c()
table_out_clean <- c()

# i = 5

for(i in 1:10){
# for(i in 1:length(mikropml_names)){
  
  table_out <- feat_import_combine[grepl(mikropml_names[i], feat_import_combine[,"names"]),]
  
  string_clean <- paste("^", mikropml_names[i], "$|\\Q|\\E", mikropml_names[i], "$|^", mikropml_names[i], "\\Q|\\E|\\Q|\\E", mikropml_names[i], "\\|", sep = "")
  
  table_out_clean <- table_out[grepl(string_clean, table_out[,"names"]),] %>%
    mutate(names = mikropml_names[i])
  
  feat_importance_collapsed <- rbind(feat_importance_collapsed,
                                     table_out_clean)
  
}

feat_import_collapsed_summary <- feat_importance_collapsed %>%
  group_by(names, GWAS_test, group, genome, GWAS_method, perf_metric_name) %>%
  summarize(mean_perf_metric_diff = mean(perf_metric_diff),
            sd_perf_metric_diff = sd(perf_metric_diff))

#summarize model performances
mikropml_perf_summary <- perf_combine %>%
  group_by(method, group, genome) %>%
  summarize(mean_cv_metric_RMSE = mean(cv_metric_RMSE),
            sd_cv_metric_RMSE = sd(cv_metric_RMSE),
            mean_RMSE = mean(RMSE),
            sd_RMSE = sd(RMSE),
            mean_Rsquared = mean(Rsquared),
            sd_Rsquared = sd(Rsquared),
            mean_MAE = mean(MAE),
            sd_MAE = sd(MAE)) %>%
  rename("GWAS_test" = "method")

#combine mikropml data
mikropml_col <- intersect(colnames(feat_import_collapsed_summary), colnames(mikropml_perf_summary))

mikropml_combine <- feat_import_collapsed_summary %>%
  full_join(mikropml_perf_summary, by = all_of(mikropml_col))

#HOGWASH ----
#Create the data frames from the hogwash data
create_df <- function(x){
  a <- data.frame(x$hi_confidence_transition_edge)
  
  conf <- tibble("names" = colnames(a),
                 "High_Confidence_Transition" = 1) %>%
    mutate(names = gsub("\\.\\.\\.", "~~~", names))
  
  df_1 <- cbind(x$hit_pvals,
                x$convergence$geno_beta,
                x$convergence$epsilon,
                x$convergence$N) %>%
    rownames_to_column(var = "names") %>%
    left_join(conf, by = "names") %>%
    mutate("Hogwash_Sig_Pval" = if_else(names %in% rownames(x$sig_hits), 1, 0))
  
  colnames(df_1) <- c("names",
                      "Hit_FDR_corrected_pvalues",
                      "Geno_Beta",
                      "Epsilon",
                      "N",
                      "High_Confidence_Transition",
                      "Hogwash_Sig_Pval")
  
  return(df_1)
}

hogwash_list <- list.files(path = result_path,
                           pattern = "hogwash")
hogwash_list <- hogwash_list[grepl(".rda", hogwash_list)]
hogwash_list <- hogwash_list[grepl(snakemake@wildcards[['genome']], hogwash_list)]
hogwash_combine <- c()

for (i in 1:length(hogwash_list)){
  load(file=paste0(result_path, "/" , hogwash_list[i]))
  
  a <- load(file=paste0(result_path, "/" , hogwash_list[i]))
  
  assign(hogwash_list[i],get(a))
  
  file_name <- gsub("hogwash_continuous_", "", hogwash_list[i])
  
  file_name <- str_split(file_name, pattern = "\\.", simplify = TRUE)
  
  df <- create_df(get(hogwash_list[i]))
  
  df_add <- df %>%
    mutate(GWAS_method = "hogwash",
           GWAS_test = "ungrouped",
           group = file_name[1],
           genome = file_name[2])
  
  hogwash_combine <- rbind(hogwash_combine,
                            df_add)
}

rm(a)
rm(df)
rm(df_add)
rm(hogwash_continuous)
rm(list = hogwash_list)

#TREEWAS ----
treeWAS_list <- list.files(path = result_path,
                           pattern = "treeWAS.RData")
treeWAS_list <- treeWAS_list[grepl(snakemake@wildcards[['genome']], treeWAS_list)]
treeWAS_combine <- c()

treeWAS_test <- c("simultaneous", "terminal", "subsequent")

int <- c()

for (i in 1:length(treeWAS_list)){
  load(file=paste0(result_path, "/", treeWAS_list[i]))
  a<-load(file=paste0(result_path, "/", treeWAS_list[i]))
  assign(treeWAS_list[i],get(a))
  
  file_name <- str_split(treeWAS_list[i], pattern = "\\.", simplify = TRUE)
  
  int <- c()
  
  for(j in 1:length(treeWAS_test)){
    
    test_out <- get(treeWAS_list[i])[[treeWAS_test[j]]][["sig.snps"]]
    
    if(is.data.frame(test_out) == TRUE){
    
      df <- test_out %>%
        rownames_to_column(var = "names") %>%
        mutate(GWAS_method = "treeWAS",
               GWAS_test = treeWAS_test[j],
               group = file_name[1],
               genome = file_name[2])
      
      int <- rbind(int,
                   df)
      
    }
    
  }
  
  treeWAS_combine <- rbind(treeWAS_combine,
                           int)
}

rm(a)
rm(df)
rm(df_add)
rm(file_name)
rm(int)
# rm(treeWAS.out)
# rm(list = treeWAS_list)

#COMBINE ----
all_col <- unique(c(intersect(colnames(treeWAS_combine), colnames(hogwash_combine)),
                    intersect(colnames(treeWAS_combine), colnames(mikropml_combine)),
                    intersect(colnames(hogwash_combine), colnames(mikropml_combine))))

combined_outputs <- mikropml_combine %>%
  full_join(hogwash_combine, by = all_of(all_col)) %>%
  full_join(treeWAS_combine, by = all_of(all_col))


