library(dplyr)
library(tibble)
library(readxl)
library(writexl)
library(drc)
library(rstatix)
library(phytools)
library(nlme)
library(factoextra)
library(ggpubr)
library(phytools)
library(ape)
library(qualpalr)
library(ggrepel)
library(purrr)
library(patchwork)



# Normalization function using ((B2-B7)-B8)/(B9-B8)
normalize_dataset <- function(df) {
  if (nrow(df) < 8) stop("Data frame must have at least 8 rows for normalization.")
  
  df %>% mutate(across(everything(), ~ {
    col_vec <- .x
    B7 <- col_vec[6]
    B8 <- col_vec[7]
    B9 <- col_vec[8]
    (col_vec - B7 - B8) / (B9 - B8)
  }))
}

# function to process one excel file
process_file <- function(filepath) {
  file_id <- tools::file_path_sans_ext(basename(filepath))
  
  PlateTest <- read_xlsx(filepath, col_names = TRUE)
  
  # Remove first and 14th columns
  Plate <- PlateTest[, -c(1, 14)] |> as.data.frame()
  
  # Assign row names
  rownames(Plate) <- c("100%", "50%", "20%", "10%", "5%", "MB", "Pos", "Neg")
  
  # Assign column names
  colnames(Plate) <- c("I60_A","I60_B","I60_C","I60_D",
                       "I39_A","I39_B","I39_C","I39_D",
                       "I50_A","I50_B","I50_C","I50_D")
  
  # Apply normalization
  Corrected <- normalize_dataset(Plate)
  
  # Remove reference rows (MB, Pos, Neg)
  Cleaned <- Corrected[-c(6:8), ]
  
  # Add row names as a column
  df_out <- Cleaned |> rownames_to_column(var = "Row_Name")
  
  # Prefix all data columns with file ID
  df_out <- df_out |> rename_with(~ paste0(file_id, "_", .x), .cols = -Row_Name)
  
  return(df_out)
}

# function to process all excel files in a folder
process_folder_combined <- function(folder_path, output_csv = NULL) {
  files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) == 0) stop("No Excel files found in folder.")
  
  all_data <- lapply(files, process_file)
  
  # Combine all datasets by Row_Name
  combined_df <- Reduce(function(x, y) full_join(x, y, by = "Row_Name"), all_data)
  
  # Save as CSV if output path is specified
  if (!is.null(output_csv)) {
    write.csv(combined_df, output_csv, row.names = FALSE)
  }
  
  return(combined_df)
}

################################################################################
## Process all files, each folder should be separated by time point
## (one for 24 hours, the other for 72 hours)
################################################################################

folder1 <- ""
folder2 <- ""

## data processing for 24 hours
usefuldata1 <- process_folder_combined(folder1, "24hr_Absorbances.csv")

long_df <- usefuldata %>%
  pivot_longer(
    cols = -Row_Name,         # Pivot all columns except the one we just created
    names_to = "Column_Name",  # Original column names go here
    values_to = "Absorbance" # New variable name "Species_123"
  )

Doses <- long_df %>%
  filter(Absorbance >= 0)

write.csv(Doses, "Corrected_Absorbances_24hr_FINAL.csv")



## data processing for 72 hours
usefuldata2 <- process_folder_combined(folder2, "72hr_Absorbances.csv")

long_df <- usefuldata2 %>%
  pivot_longer(
    cols = -Row_Name,         # Pivot all columns except the one we just created
    names_to = "Column_Name",  # Original column names go here
    values_to = "Absorbance" # New variable name "Species_123"
  )

Doses <- long_df %>%
  filter(Absorbance >= 0)

write.csv(Doses, "Corrected_Absorbances_72hr_FINAL.csv")

################################################################################
## Curve fitting and parameter extractionf for data extracted at 24 hours
################################################################################

DoseResp <- read.csv("Corrected_Absorbances_24hr_FINAL.csv", header = T)

out_dir <- file.path(getwd(), "dose_response_plots_svg_perRep_24hr_FINAL")
dir.create(out_dir, showWarnings = FALSE)

all_params <- data.frame(
  curveID = character(),
  b = numeric(),
  c = numeric(),
  d = numeric(),
  e = numeric(),
  var_b = numeric(),
  se_b = numeric(),
  var_e = numeric(),
  se_e = numeric(),
  stringsAsFactors = FALSE
)


for (cid in curves) { 
  message("Fitting curve: ", cid) 
  
  sub <- droplevels(DoseResp_NoZero[DoseResp_NoZero$curveID == cid, ]) 
  
  model <- try(
    drm(
      Absorbance ~ LogConc, 
      data = sub, 
      fct = LL2.5()
    ),
    silent = TRUE
  ) 
  
  if (inherits(model, "try-error")) { 
    warning("Model failed for curve ", cid) 
    next 
  } 
  
  ## ---- Extract parameters ----
  cf <- coef(model)
  
  ## ---- Extract variance–covariance matrix safely ----
  vc <- try(vcov(model), silent = TRUE)
  
  if (inherits(vc, "try-error")) {
    warning("vcov failed for curve ", cid)
    var_b <- var_e <- se_b <- se_e <- NA
  } else {
    
    par_names <- names(coef(model))
    names(par_names) <- par_names  # ensure mapping
    
    idx_b <- par_names["b:(Intercept)"]
    idx_e <- par_names["e:(Intercept)"]
    
    idx_b <- which(par_names == "b:(Intercept)")
    idx_e <- which(par_names == "e:(Intercept)")
    
    var_b <- vc[idx_b, idx_b]
    var_e <- vc[idx_e, idx_e]
    
    se_b  <- sqrt(var_b)
    se_e  <- sqrt(var_e)
  }
  
  ## ---- Store results ----
  all_params <- rbind(
    all_params,
    data.frame(
      curveID = cid,
      b = cf["b:(Intercept)"],      # slope
      c = cf["c:(Intercept)"],
      d = cf["d:(Intercept)"],
      e = cf["e:(Intercept)"],      # EC50
      var_b = var_b,
      se_b  = se_b,
      var_e = var_e,
      se_e  = se_e
    )
  )
  
  ## ---- Generate prediction grid ----
  newdata <- data.frame(
    LogConc = seq(min(sub$LogConc), max(sub$LogConc), length.out = 200)
  )
  newdata$Absorbance <- predict(model, newdata)
  
  ## ---- ggplot ----
  p <- ggplot(sub, aes(LogConc, Absorbance)) + 
    geom_point(size = 2) + 
    geom_line(data = newdata, aes(LogConc, Absorbance), linewidth = 1) +
    labs(
      title = paste("Dose–Response Curve:", cid),
      x = "LogConc",
      y = "Absorbance"
    ) + 
    theme_bw()
  
  ## ---- Save as SVG ----
  ggsave(
    filename = file.path(out_dir, paste0("curve_", cid, ".svg")),
    plot = p,
    width = 6,
    height = 4,
    units = "in"
  )
}


all_params

write.csv(all_params, "dose_response_parameters_perRep_24hr_FINAL.csv", row.names = FALSE)

UsefulDoses <- read.csv("dose_response_parameters_perRep_72hr_FINAL.csv", header = T)

SpeciesList <- read.csv("Codes.csv", header = TRUE)  ## used to convert the JAD code to the species name, provided in the data file


NewData <- left_join(SpeciesList, UsefulDoses, by = c(SpeciesList$Sample, UsefulDoses$Sample))

merged_df_base <- merge(x = SpeciesList, y = UsefulDoses, by.x = "Sample", by.y = "Sample", all.x = TRUE)

df_cleaned <- na.omit(merged_df_base)

write.csv(df_cleaned, "Final_72hr_GrowthParemeters.csv")


################################################################################
## Curve fitting and parameter extractionf for data extracted at 72 hours
################################################################################

DoseResp <- read.csv("Corrected_Absorbances_72hr_FINAL.csv", header = T)

out_dir <- file.path(getwd(), "dose_response_plots_svg_perRep_72hr_FINAL")
dir.create(out_dir, showWarnings = FALSE)

all_params <- data.frame(
  curveID = character(),
  b = numeric(),
  c = numeric(),
  d = numeric(),
  e = numeric(),
  var_b = numeric(),
  se_b = numeric(),
  var_e = numeric(),
  se_e = numeric(),
  stringsAsFactors = FALSE
)


for (cid in curves) { 
  message("Fitting curve: ", cid) 
  
  sub <- droplevels(DoseResp_NoZero[DoseResp_NoZero$curveID == cid, ]) 
  
  model <- try(
    drm(
      Absorbance ~ LogConc, 
      data = sub, 
      fct = LL2.5()
    ),
    silent = TRUE
  ) 
  
  if (inherits(model, "try-error")) { 
    warning("Model failed for curve ", cid) 
    next 
  } 
  
  ## ---- Extract parameters ----
  cf <- coef(model)
  
  ## ---- Extract variance–covariance matrix safely ----
  vc <- try(vcov(model), silent = TRUE)
  
  if (inherits(vc, "try-error")) {
    warning("vcov failed for curve ", cid)
    var_b <- var_e <- se_b <- se_e <- NA
  } else {
    
    par_names <- names(coef(model))
    names(par_names) <- par_names  # ensure mapping
    
    idx_b <- par_names["b:(Intercept)"]
    idx_e <- par_names["e:(Intercept)"]
    
    idx_b <- which(par_names == "b:(Intercept)")
    idx_e <- which(par_names == "e:(Intercept)")
    
    var_b <- vc[idx_b, idx_b]
    var_e <- vc[idx_e, idx_e]
    
    se_b  <- sqrt(var_b)
    se_e  <- sqrt(var_e)
  }
  
  ## ---- Store results ----
  all_params <- rbind(
    all_params,
    data.frame(
      curveID = cid,
      b = cf["b:(Intercept)"],      # slope
      c = cf["c:(Intercept)"],
      d = cf["d:(Intercept)"],
      e = cf["e:(Intercept)"],      # EC50
      var_b = var_b,
      se_b  = se_b,
      var_e = var_e,
      se_e  = se_e
    )
  )
  
  ## ---- Generate prediction grid ----
  newdata <- data.frame(
    LogConc = seq(min(sub$LogConc), max(sub$LogConc), length.out = 200)
  )
  newdata$Absorbance <- predict(model, newdata)
  
  ## ---- ggplot ----
  p <- ggplot(sub, aes(LogConc, Absorbance)) + 
    geom_point(size = 2) + 
    geom_line(data = newdata, aes(LogConc, Absorbance), linewidth = 1) +
    labs(
      title = paste("Dose–Response Curve:", cid),
      x = "LogConc",
      y = "Absorbance"
    ) + 
    theme_bw()
  
  ## ---- Save as SVG ----
  ggsave(
    filename = file.path(out_dir, paste0("curve_", cid, ".svg")),
    plot = p,
    width = 6,
    height = 4,
    units = "in"
  )
}


all_params

write.csv(all_params, "dose_response_parameters_perRep_72hr_FINAL.csv", row.names = FALSE)


UsefulDoses <- read.csv("dose_response_parameters_perRep_72hr_FINAL.csv", header = T)

SpeciesList <- read.csv("Codes.csv", header = TRUE)  ## used to convert the JAD code to the species name, provided in the data file


NewData <- left_join(SpeciesList, UsefulDoses, by = c(SpeciesList$Sample, UsefulDoses$Sample))

merged_df_base <- merge(x = SpeciesList, y = UsefulDoses, by.x = "Sample", by.y = "Sample", all.x = TRUE)

df_cleaned <- na.omit(merged_df_base)

write.csv(df_cleaned, "Final_72hr_GrowthParemeters.csv")

################################################################################
## Data Analysis
################################################################################

## Combine the two data sets generated manually, adding in the time points as a new collumn titled "Time",
## this dataset is provided as Final_DoseResponse_Rhododendron.csv

UsefulDoses <- read.csv("Final_DoseResponse_RhododendronData.csv", header = T)

## boxplot analysis species by time and by dose parameter ######################

A <- ggplot(UsefulDoses, aes(x = fct_rev(Species), y = Slope))+
  geom_boxplot(outlier.color = "red", aes(fill = Species))+
  coord_flip()+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  geom_hline(yintercept = c(-1, 1), color = "goldenrod", linetype = "solid")+
  scale_fill_manual(values = pal$hex)+
  facet_grid(.~Time)+
  scale_color_manual(values = pal$hex)+
  theme_classic(base_size = 15)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
A  



B <- ggplot(UsefulDoses, aes(x = fct_rev(Species), y = pEC50)) +
  geom_boxplot(outlier.color = "red", aes(fill = Species))+
  coord_flip()+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  #geom_hline(yintercept = c(-1, 1), color = "goldenrod", linetype = "solid")+
  scale_fill_manual(values = pal$hex)+
  facet_grid(.~Time)+
  scale_color_manual(values = pal$hex)+
  theme_classic(base_size = 15)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

B

ggsave("Slope_Boxplot.svg", A, height = 15, width = 20)
ggsave("pEC50_Boxplot.svg", B, height = 15, width = 20)



grouped_tests_results_pEC50_Time_Species <- UsefulDoses %>%
  group_by(Time) %>%
  wilcox_test(pEC50 ~ Species) %>%
  adjust_pvalue(method = "BH") %>%  # Adjust p-values for multiple comparisons
  add_significance()

write.csv("Wilcoxson_TestResults_pEC50byTime.csv", grouped_tests_results_pEC50_Time_Species, header = T)

grouped_tests_results_Slope_Time_Species <- UsefulDoses %>%
  group_by(Time) %>%
  wilcox_test(Slope ~ Species) %>%
  adjust_pvalue(method = "BH") %>%  # Adjust p-values for multiple comparisons
  add_significance()

write.csv("Wilcoxson_TestResults_SlopebyTime.csv", grouped_tests_results_pEC50_Time_Species, header = T)



## Violin Plots for pEC50 and Isolate Effects ##################################

my_comparisons <- combn(unique(UsefulDoses$Isolate), 2, simplify = FALSE)


pEC50_Isolate <- ggplot(UsefulDoses, aes(y = pEC50, x =  Isolate))+
  geom_violin(aes(fill = Isolate))+
  #scale_y_continuous(breaks = c(-5, -2, -1, -0.5, 0, 0.5, 1, 2, 5))+
  stat_summary(
    aes(group = Isolate),  # separate summaries per species
    fun.data = "mean_cl_boot",
    geom = "pointrange",
    color = "black",
    position = position_dodge(width = 0.9),
    linewidth = 0.6
  ) +
  stat_compare_means(
    comparisons = my_comparisons, 
    method = "wilcox.test", 
    label = "p.signif",
  ) +
  theme_classic(base_size = 15)


pEC50_Isolate_Time <- ggplot(UsefulDoses, aes(y = pEC50, x =  Isolate))+
  geom_violin(aes(fill = Isolate))+
  #scale_y_continuous(breaks = c(-5, -2, -1, -0.5, 0, 0.5, 1, 2, 5))+
  facet_grid(.~Time)+
  stat_summary(
    aes(group = Isolate),  # separate summaries per species
    fun.data = "mean_cl_boot",
    geom = "pointrange",
    color = "black",
    position = position_dodge(width = 0.9),
    linewidth = 0.6
  ) +
  stat_compare_means(
    comparisons = my_comparisons, 
    method = "wilcox.test", 
    label = "p.signif",
  ) +
  theme_classic(base_size = 15)


stat.test <- UsefulDoses %>%
  group_by(Isolate) %>%
  wilcox_test(pEC50 ~ Time) %>%
  add_significance("p") %>%
  ungroup() %>%
  mutate(
    xmin = as.numeric(factor(Isolate)) - 0.2,
    xmax = as.numeric(factor(Isolate)) + 0.2,
    y.position = max(UsefulDoses$pEC50, na.rm = TRUE) * 1.1
  )


pEC50_Time_Paired_Whole <- ggplot(UsefulDoses, aes(y = pEC50, x = Isolate, fill = as.factor(Time)))+
  geom_violin()+
  stat_summary(
    aes(group = Time),  # separate summaries per species
    fun.data = "mean_cl_boot",
    geom = "pointrange",
    color = "black",
    position = position_dodge(width = 0.9),
    linewidth = 0.6
  ) +
  stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    inherit.aes = FALSE,
    tip.length = 0.02
  ) +
  theme_classic(base_size = 15)



paired_data <- UsefulDoses %>%
  select(Isolate, Sample, Time, pEC50) %>%
  pivot_wider(
    names_from = Time,
    values_from = pEC50
  )

time_cols <- setdiff(
  names(paired_data),
  c("Isolate", "Sample")
)

stat.test <- paired_data %>%
  filter(if_all(all_of(time_cols), ~ !is.na(.))) %>%
  group_by(Isolate) %>%
  summarise(
    p = wilcox.test(
      .data[[time_cols[1]]],
      .data[[time_cols[2]]],
      paired = TRUE
    )$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    p.signif = case_when(
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  mutate(
    group1 = time_cols[1],
    group2 = time_cols[2],
    xmin = 1,
    xmax = 2,
    y.position = max(UsefulDoses$pEC50, na.rm = TRUE) * 1.1
  )

stat.test <- stat.test %>%
  left_join(
    UsefulDoses %>%
      group_by(Isolate) %>%
      summarise(
        ymax = max(pEC50, na.rm = TRUE),
        ymin = min(pEC50, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Isolate"
  ) %>%
  mutate(
    y.position = ymax + (ymax - ymin) * 0.15
  )

pEC50_Time_Paired_Spaghetti <- ggplot(UsefulDoses, aes(y = pEC50, x = as.factor(Time))) +
  geom_violin(aes(fill = as.factor(Time)), alpha = 0.85) +
  geom_line(aes(group = Sample), color = "grey50", alpha = 0.6, linewidth = 0.6) +
  geom_point(color = "grey50", size = 2.5, alpha = 0.6) +
  facet_wrap(~ Isolate) +
  stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    inherit.aes = FALSE,
    tip.length = 0.02
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Time",
    y = "pEC50"
  )



## Violin Plots for Slope and Isolate Effects ##################################

Slope_Isolate <- ggplot(UsefulDoses, aes(y = Slope, x =  Isolate))+
  geom_violin(aes(fill = Isolate))+
  #scale_y_continuous(breaks = c(-5, -2, -1, -0.5, 0, 0.5, 1, 2, 5))+
  stat_summary(
    aes(group = Isolate),  # separate summaries per species
    fun.data = "mean_cl_boot",
    geom = "pointrange",
    color = "black",
    position = position_dodge(width = 0.9),
    linewidth = 0.6
  ) +
  stat_compare_means(
    comparisons = my_comparisons, 
    method = "wilcox.test", 
    label = "p.signif",
  ) +
  theme_classic(base_size = 15)


Slope_Isolate_Time <- ggplot(UsefulDoses, aes(y = Slope, x =  Isolate))+
  geom_violin(aes(fill = Isolate))+
  facet_grid(.~Time)+
  stat_summary(
    aes(group = Isolate),  # separate summaries per species
    fun.data = "mean_cl_boot",
    geom = "pointrange",
    color = "black",
    position = position_dodge(width = 0.9),
    linewidth = 0.6
  ) +
  stat_compare_means(
    comparisons = my_comparisons, 
    method = "wilcox.test", 
    label = "p.signif",
  ) +
  theme_classic(base_size = 15)


stat.test <- UsefulDoses %>%
  group_by(Isolate) %>%
  wilcox_test(Slope ~ Time) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup() %>%
  mutate(
    xmin = as.numeric(factor(Isolate)) - 0.2,
    xmax = as.numeric(factor(Isolate)) + 0.2,
    y.position = max(UsefulDoses$Slope, na.rm = TRUE) * 1.1
  )
  

Slope_Time_Paired_Whole <- ggplot(UsefulDoses, aes(y = Slope, x = Isolate, fill = as.factor(Time)))+
  geom_violin()+
  stat_summary(
    aes(group = Time),  # separate summaries per species
    fun.data = "mean_cl_boot",
    geom = "pointrange",
    color = "black",
    position = position_dodge(width = 0.9),
    linewidth = 0.6
  ) +
  stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    inherit.aes = FALSE,
    tip.length = 0.02
  ) +
  theme_classic(base_size = 15)



paired_data <- UsefulDoses %>%
  select(Isolate, Sample, Time, Slope) %>%
  pivot_wider(
    names_from = Time,
    values_from = Slope
  )

time_cols <- setdiff(
  names(paired_data),
  c("Isolate", "Sample")
)

stat.test <- paired_data %>%
  filter(if_all(all_of(time_cols), ~ !is.na(.))) %>%
  group_by(Isolate) %>%
  summarise(
    p = wilcox.test(
      .data[[time_cols[1]]],
      .data[[time_cols[2]]],
      paired = TRUE
    )$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    p.signif = case_when(
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  mutate(
    group1 = time_cols[1],
    group2 = time_cols[2],
    xmin = 1,
    xmax = 2,
    y.position = max(UsefulDoses$Slope, na.rm = TRUE) * 1.1
  )

stat.test <- stat.test %>%
  left_join(
    UsefulDoses %>%
      group_by(Isolate) %>%
      summarise(
        ymax = max(Slope, na.rm = TRUE),
        ymin = min(Slope, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Isolate"
  ) %>%
  mutate(
    y.position = ymax + (ymax - ymin) * 0.15
  )

Slope_Time_Paired_Spaghetti <- ggplot(UsefulDoses, aes(y = Slope, x = as.factor(Time))) +
  geom_violin(aes(fill = as.factor(Time)), alpha = 0.85) +
  geom_line(aes(group = Sample), color = "grey50", alpha = 0.6, linewidth = 0.6) +
  geom_point(color = "grey50", size = 2.5, alpha = 0.6) +
  facet_wrap(~ Isolate) +
  stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    inherit.aes = FALSE,
    tip.length = 0.02
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Time",
    y = "Slope"
  )






## Chromophore Correlation Anaylsis within Rhododendron ########################

genus_name <- "Rhododendron"

tree_path <- file.path(getwd(), paste0(genus_name, "_GraftedTree_RhodoTox.tre"))
peaktable_path <- file.path(getwd(), paste0(genus_name, "_AveragePeaks_NewMethod_ForRhodoTox.csv"))
chromophore_path <- file.path(getwd(), paste0(genus_name, "_AllMetabolomicData_byChromophore_ForRhodoTox.csv"))
ecophys_path <- file.path(getwd(),  paste0(genus_name, "_Ecophys_RhodoTox.csv"))


## General Summary Statistics for Chromophore Abundances
Chromophores <- read.csv(chromophore_path)


df <- Chromophores[, -1]



summarize_columns <- function(df) {
  # Keep only numeric columns
  num_df <- df[sapply(df, is.numeric)]
  
  # Apply summary statistics across columns
  res <- sapply(num_df, function(x) {
    c(
      Mean  = mean(x, na.rm = TRUE),
      SD    = sd(x, na.rm = TRUE),
      SE    = sd(x, na.rm = TRUE)/sqrt(as.numeric(length(x))),
      Min   = min(x, na.rm = TRUE),
      Max   = max(x, na.rm = TRUE),
      Range = diff(range(x, na.rm = TRUE))
    )
  })
  
  # Transpose and convert to a data frame
  as.data.frame(t(res))
}

# Example usage with the built-in mtcars dataset:
summary_table <- summarize_columns(Chromophores[, -1])
print(summary_table)




## Pie charts for Chromophore abundances
norm_path <- "C:/Users/justi/OneDrive - Louisiana State University/Putnam_Reimagined/normalized_data/Finished"

genus <- "Rhododendron"

genus_file <- paste0(norm_path, "/",genus, "_FinalChemData_ForRhodoTox.csv")

Genus_Data <- read.csv(genus_file, header = T)

chromophore <- Genus_Data[-1,]$Chromophore
class_label <- Genus_Data[-1, ]$Class
peak_label <- Genus_Data[-1, ]$Sample

names_to_remove <- c("Chromophore", "Sample", "Class")

Genus_Use <- Genus_Data[, !colnames(Genus_Data) %in% names_to_remove]

# Row 1 contains the grouping labels
group_labels <- as.character(Genus_Use[1, ])

# Remove row 1
numeric_data <- Genus_Use[-1, , drop = FALSE]

# Convert all remaining data to numeric
numeric_data <- numeric_data %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

# Find unique groups
groups <- unique(group_labels)

# Average columns belonging to each group
result <- sapply(groups, function(group) {
  
  # Identify columns belonging to this group
  cols <- group_labels == group
  
  # Calculate row-wise mean
  rowMeans(
    numeric_data[, cols, drop = FALSE],
    na.rm = TRUE
  )
})

# Convert to data frame
result <- as.data.frame(result)

# Preserve original row names
rownames(result) <- rownames(numeric_data)

# Give the resulting columns the group names
colnames(result) <- groups

Final_Averages <- cbind(peak_label, class_label, chromophore, result)

Class_Sums <- Final_Averages %>%
  group_by(class_label) %>%
  summarise(
    across(
      where(is.numeric),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

Class_Pivot <- Class_Sums %>%
  pivot_longer(
    cols = -class_label,
    names_to = "Species",
    values_to = "value"
  )


Class_Pie <- Class_Pivot %>%
  group_by(class_label) %>%
  summarise(
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  )


Pie1 <- ggplot(Class_Pie, aes(x = "", y = value, fill = class_label)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y", start = 0) +
  theme_void() + # Removes background lines and axes grids
  labs(fill = "Chemical Class") +
  scale_fill_brewer(palette = "Set2") # Clean color scheme


Class_Data <- as.data.frame(t(Class_Sums))

colnames(Class_Data) <- Class_Sums$class_label
Class_Data_Final <- Class_Data[-1, ]

Class_Data_Final <- Class_Data_Final[order(rownames(Class_Data_Final)), ]


Chromophore_Sums <- Final_Averages %>%
  group_by(chromophore) %>%
  summarise(
    across(
      where(is.numeric),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )


Chromophore_Pivot <- Chromophore_Sums %>%
  pivot_longer(
    cols = -chromophore,
    names_to = "Species",
    values_to = "value"
  )


counts <- as.data.frame(table(chromophore))

Chromophore_Pie <- Chromophore_Pivot %>%
  group_by(chromophore) %>%
  summarise(
    total_amount = sum(value, na.rm = TRUE),
    n_individuals = n(),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = (total_amount / sum(total_amount)) * 100,
    round_percent = round(percentage, 1),
    csum = rev(cumsum(rev(round_percent))), 
    pos = round_percent/2 + lead(csum, 1),
    pos = if_else(is.na(pos), round_percent/2, pos)
  )

Chromophore_Pie <- merge(Chromophore_Pie, counts, 1)

Chromophore_Pie2 <- Chromophore_Pie%>%
  mutate(
    csum_freq = rev(cumsum(rev(Freq))), 
    pos_freq = Freq / 2 + lead(csum_freq, 1),
    pos_freq = if_else(is.na(pos_freq), Freq / 2, pos_freq)
  )


cols <- qualpal(length((unique(Chromophore_Sums$chromophore))))

Pie2 <- ggplot(Chromophore_Pie, aes(x = "", y = round_percent, fill = chromophore)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_label_repel(data = Chromophore_Pie,
                   aes(y = pos, label = paste0(as.numeric(round(total_amount, 1)), ", ", round_percent, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  coord_polar(theta = "y", start = 0) +
  theme_void() + # Removes background lines and axes grids
  labs(fill = "Chromophore Class") +
  scale_fill_manual(values = cols$hex) # Clean color scheme


Chromophore_Pie2 <- Chromophore_Pie%>%
  mutate(
    percentage_compounds = (Freq / sum(Freq)) * 100,
    round_percent_freq = round(percentage_compounds, 1),
    csum_freq = rev(cumsum(rev(round_percent_freq))), 
    pos_freq = round_percent_freq / 2 + lead(csum_freq, 1),
    pos_freq = if_else(is.na(pos_freq), Freq / 2, pos_freq)
  )

Pie3 <- ggplot(Chromophore_Pie2, aes(x = "", y = round_percent_freq, fill = chromophore)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_label_repel(data = Chromophore_Pie2,
                   aes(y = pos_freq, label = paste0(Freq, ", ", round_percent_freq, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  coord_polar(theta = "y", start = 0) +
  theme_void() + # Removes background lines and axes grids
  labs(fill = "Chromophore Class") +
  scale_fill_manual(values = cols$hex) # Clean color scheme

horiz_pies <- ggarrange(Pie3, Pie2, ncol = 2, nrow = 1, labels = c("A)", "B)"), common.legend = TRUE, legend = "bottom")
ggsave("Horiz_RhodoTox_Pie.svg", horiz_pies, height = 10, width = 15)
vertic_pies <- ggarrange(Pie3, Pie2, ncol = 1, nrow = 2, labels = c("A)", "B)"), common.legend = TRUE, legend = "right")
ggsave("Vertical_RhodoTox_Pie.svg", vertic_pies, height = 15, width = 10)








macro_chemicals <- read.csv(chromophore_path)
colnames(macro_chemicals)[1] <- "TreeTipName"

macro_chemicals_noNA <- na.omit(macro_chemicals)

chem_tree<-drop.tip(tree,setdiff(tree$tip.label,macro_chemicals_noNA$TreeTipName))

names_to_remove <-setdiff(macro_chemicals_noNA$TreeTipName, tree$tip.label)

if(length(names_to_remove) > 0) {
  macro_chemicals_noNA <- macro_chemicals_noNA[!(macro_chemicals_noNA$TreeTipName %in% names_to_remove), ]
} else{
  macro_chemicals_noNA <- macro_chemicals_noNA
}


# ---- Identify numeric trait columns ----
trait_cols <- names(macro_chemicals_noNA)[sapply(macro_chemicals, is.numeric)]

# ---- Identify numeric trait columns ----
trait_cols <- names(macro_chemicals_noNA)[sapply(macro_chemicals, is.numeric)]

# ---- Initialize results storage ----
results_list <- list()

data <- macro_chemicals_noNA
pmat <- r2mat <- rmat <- matrix(NA,ncol(data)-1,ncol(data)-1,dimnames = list(colnames(data)[-1],colnames(data[-1])))

for(i in 2:ncol(data)){
  for(j in 2:ncol(data))
  {
    if(i==j) next
    cat(paste(i,j,ncol(data)," "))
    dat <- data.frame(species=data[,1],a=data[,i],b=data[,j])
    dat <- dat[which(complete.cases(dat[,2:3])),]
    dat[,1] <- as.character(dat[,1])
    mod <- phylopars(trait_data = dat,tree = tree,phylo_correlated = FALSE,pheno_correlated = FALSE)
    ab <- mod$anc_recon[1:length(tree$tip.label),]
    mod <- summary(lm(pic(ab[,1],tree)~pic(ab[,2],tree)-1))
    
    p <- mod$coefficients[4]
    r2 <- mod$r.squared
    r <- sqrt(r2)*sign(coef(mod)[1])
    
    
    pmat[i-1,j-1] <- p
    r2mat[i-1,j-1] <- r2
    rmat[i-1,j-1] <- r
  }
}

#Generate Correlation Figure
Cor.Matrix <- rmat
fdrs <- p.adjust(pmat, method = "BH")
resAdj <- matrix(fdrs, ncol = dim(pmat))
rownames(resAdj) <- colnames(data[, -1])
colnames(resAdj) <- colnames(data[, -1])

# generate phylo-corrected correlation matrix and plot
phylo_corplot_path <- file.path(paste0(genus_name, "_MacroEvo_Chromophore.svg"))

svg(phylo_corplot_path, height = 8, width = 8)
corrplot(Cor.Matrix, method = "color", col = magma(200), diag = FALSE, insig = "blank",
         p.mat = resAdj, sig.level = 0.05, tl.cex = 1, tl.col = "black", order = "hclust", addgrid.col = NA,  outline = FALSE)
dev.off()



################################################################################
# Total Effects of Chromophores on Dose Response Parameters
################################################################################

## Setting up the 24hr subset
DoseData <- read.csv("Final_DoseResponse_RhododendronData.csv", header = T)

DoseData$Isolate <- factor(
  DoseData$Isolate,
  levels = unique(DoseData$Isolate)
)

DoseData$Species <- factor(
  DoseData$Species,
  levels = unique(DoseData$Species)
)


DoseData24<- DoseData %>%
  filter(Time == 24) %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))


Dose24_Slope <- DoseData24$Slope
names(Dose24_Slope) <- DoseData24$Species

Dose24_pEC50 <- DoseData24$pEC50
names(Dose24_pEC50) <- DoseData24$Species

Dose24_EC50 <- DoseData24$EC50
names(Dose24_EC50) <- DoseData24$Species

# Phylogenetic Signal at the dose-response parameter level
phylosig(tree, Dose24_Slope, method = "K", test = TRUE)
phylosig(tree, Dose24_pEC50, method = "K", test = TRUE)
phylosig(tree, Dose24_EC50, method = "K", test = TRUE)

phylosig(tree, Dose24_Slope, method = "lambda", test = TRUE)
phylosig(tree, Dose24_pEC50, method = "lambda", test = TRUE)
phylosig(tree, Dose24_EC50, method = "lambda", test = TRUE)






## PGLS Effect Modeling at 24 hours ############################################
DoseData24_NoMean <- DoseData %>%
  filter(Time == 24)

# Slope predicting pEC50
Model1 <- gls(pEC50 ~ Slope, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- summary(Model1)$tTable
write.csv(stats, "Model1_pEC50vsSlope_24.csv")

Model2 <- gls(pEC50 ~ Slope + 0 + Isolate, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model2)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model2_pEC50vsSlope_24.csv")

Model3 <- gls(pEC50 ~ Slope + 0 + Species, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model3)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model3_pEC50vsSlope_24.csv")

Model4 <- gls(pEC50 ~ Slope + 0 + Isolate + Species, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model4)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model4_pEC50vsSlope_24.csv")


## pEC50 predicting Slope
Model1 <- gls(Slope ~ pEC50, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- summary(Model1)$tTable
write.csv(stats, "Model1_SlopevspEC50_24.csv")

Model2 <- gls(Slope ~ pEC50 + 0 + Isolate, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model2)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model2_SlopevspEC50_24.csv")

Model3 <- gls(Slope ~ pEC50 + 0 + Species, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model3)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model3_SlopevspEC50_24.csv")

Model4 <- gls(Slope ~ pEC50 + 0 + Isolate + Species, data = DoseData24_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model4)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model4_SlopevspEC50_24.csv")



## Chromophore Effects on Dose Response Parameters at 24 hr

macro_chemicals <- read.csv(chromophore_path)

DoseData24$Species <- gsub("\\. ", ".", DoseData24$Species)

ChemRed <- setdiff(macro_chemicals$X, DoseData24$Species)

ChemData <- macro_chemicals[!macro_chemicals$X %in% ChemRed, ]


AllData_24 <- merge(ChemData, DoseData24[, -2], by = 1)

colnames(AllData_24)[1] <- "TreeTipName"


# ---- Identify numeric trait columns ----
trait_cols <- names(AllData_24)[sapply(AllData_24, is.numeric)]

# ---- Initialize results storage ----
results_list <- list()

# ---- Loop through all pairwise combinations ----
counter <- 1


PCA_Parameters <- colnames(ChemData[-1])

Final_data <- list()

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- AllData_24[, c(Tox_Var, PCA_Var, "TreeTipName")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric(scale((df[[Tox_Var]])))
    df[[PCA_Var]] <- as.numeric(scale((df[[PCA_Var]])))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model using the standardized data
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~TreeTipName
      ),
      method = "REML"
    )
    
    # Extract coefficient table
    coef_df <- as.data.frame(
      summary(Model)$tTable[-1, c(
        "Value",
        "Std.Error",
        "p-value"
      ), drop = FALSE]
    )
    
    coef_df$Variable <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add identifying information
    coef_df$Tox_Parameter <- Tox_Var
    coef_df$PCA_Parameter <- PCA_Var
    coef_df$Model <- deparse(Formula)
    
    # Store results
    Final_data[[length(Final_data) + 1]] <- coef_df
  }
}

# ---- Combine into a single dataframe ----
Final_data <- do.call(rbind, Final_data)

rownames(Final_data) <- NULL


HeatmapData <- data.frame(
  X.label = Final_data$Tox_Parameter,
  y.label = Final_data$PCA_Parameter,
  Beta = Final_data$Value,
  Pvalue = Final_data$`p-value`,
  Model = Final_data$Model
)


HeatmapData <- HeatmapData %>%
  adjust_pvalue(
    p.col = "Pvalue",
    output.col = "Padj",
    method = "bonferroni"
  )


HeatmapData$signif <- HeatmapData$Padj <= 0.05

HeatmapData$signif[HeatmapData$Signif == FALSE] <- NA_real_

HeatmapData <- HeatmapData %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

write.csv(HeatmapData, "Rhododendron_PGLS_DoseXChemDiv_24hr.csv")

HeatmapData$y.label <- factor(HeatmapData$y.label, levels = unique(HeatmapData$y.label))

# plot for all parameters
megaPlot_24 <- ggplot(data = HeatmapData, aes(x = X.label, y = fct_rev(y.label), fill = Beta_plot))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )



# plot for only Slope and pEC50
DoseBetas24 <- HeatmapData %>%
  filter(X.label %in% c("Slope", "pEC50"))


plots_24 <- ggplot(DoseBetas24,aes(x = X.label, y = fct_rev(y.label), fill = Beta_plot)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_fixed()

plots_24

#only significant plots 
DoseBetas24_sigonly <- DoseBetas24 %>%
  filter(Padj <= 0.05)


plots_24_sigonly <- ggplot(DoseBetas24_sigonly,aes(x = X.label, y = fct_rev(y.label), fill = Beta)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue"
                       #na.value = "ivory4"
  )+
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_fixed()






















## Setting up the 72 hour set
DoseData72 <- DoseData %>%
  filter(Time == 72) %>%
  group_by(Species) %>% 
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))


DoseData72_NoMean <- DoseData %>%
  filter(Time == 72)


Dose72_Slope <- DoseData72$Slope
names(Dose72_Slope) <- DoseData72$Species

Dose72_pEC50 <- DoseData72$pEC50
names(Dose72_pEC50) <- DoseData72$Species

Dose72_EC50 <- DoseData72$EC50
names(Dose72_EC50) <- DoseData72$Species


# Phylogenetic Signal at the dose-response parameter level
phylosig(tree, Dose72_Slope, method = "K", test = TRUE)
phylosig(tree, Dose72_pEC50, method = "K", test = TRUE)
phylosig(tree, Dose72_EC50, method = "K", test = TRUE)

phylosig(tree, Dose72_Slope, method = "lambda", test = TRUE)
phylosig(tree, Dose72_pEC50, method = "lambda", test = TRUE)
phylosig(tree, Dose72_EC50, method = "lambda", test = TRUE)


## PGLS Effect Modeling at 72 hours ############################################

# Slope predicting pEC50
Model1 <- gls(pEC50 ~ Slope, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- summary(Model1)$tTable
write.csv(stats, "Model1_pEC50vsSlope_72.csv")

Model2 <- gls(pEC50 ~ Slope + 0 + Isolate, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model2)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model2_pEC50vsSlope_72.csv")

Model3 <- gls(pEC50 ~ Slope + 0 + Species, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model3)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model3_pEC50vsSlope_72.csv")

Model4 <- gls(pEC50 ~ Slope + 0 + Isolate + Species, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model4)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model4_pEC50vsSlope_72.csv")





# pEC50 predicting slope
Model1 <- gls(Slope ~ pEC50, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- summary(Model1)$tTable
write.csv(stats, "Model1_SlopevspEC50_72.csv")

Model2 <- gls(Slope ~ pEC50 + 0 + Isolate, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model2)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model2_SlopevspEC50_72.csv")

Model3 <- gls(Slope ~ pEC50 + 0 + Species, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model3)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model3_SlopevspEC50_72.csv")

Model4 <- gls(Slope ~ pEC50 + 0 + Isolate + Species, data = DoseData72_NoMean, correlation = corBrownian(1, phy = chem_tree, form = ~Species), 
              method = "REML")
stats <- as.data.frame(summary(Model4)$tTable)
stats$P.adj <- p.adjust(stats$`p-value`, method = "bonferroni")
write.csv(stats, "Model4_SlopevspEC50_72.csv")


DoseData72$Species <- gsub("\\. ", ".", DoseData72$Species)
AllData_72 <- merge(ChemData, DoseData72[, -2], by = 1)


colnames(AllData_72)[1] <- "TreeTipName"


# ---- Identify numeric trait columns ----
trait_cols <- names(AllData_72)[sapply(AllData_72, is.numeric)]

# ---- Initialize results storage ----
results_list <- list()

# ---- Loop through all pairwise combinations ----
counter <- 1


PCA_Parameters <- colnames(ChemData[-1])

Final_data <- list()

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- AllData_72[, c(Tox_Var, PCA_Var, "TreeTipName")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric(scale((df[[Tox_Var]])))
    df[[PCA_Var]] <- as.numeric(scale((df[[PCA_Var]])))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model using the standardized data
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~TreeTipName
      ),
      method = "REML"
    )
    
    # Extract coefficient table
    coef_df <- as.data.frame(
      summary(Model)$tTable[-1, c(
        "Value",
        "Std.Error",
        "p-value"
      ), drop = FALSE]
    )
    
    coef_df$Variable <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add identifying information
    coef_df$Tox_Parameter <- Tox_Var
    coef_df$PCA_Parameter <- PCA_Var
    coef_df$Model <- deparse(Formula)
    
    # Store results
    Final_data[[length(Final_data) + 1]] <- coef_df
  }
}

# ---- Combine into a single dataframe ----
Final_data <- do.call(rbind, Final_data)

rownames(Final_data) <- NULL


HeatmapData <- data.frame(
  X.label = Final_data$Tox_Parameter,
  y.label = Final_data$PCA_Parameter,
  Beta = Final_data$Value,
  Pvalue = Final_data$`p-value`,
  Model = Final_data$Model
)


HeatmapData <- HeatmapData %>%
  adjust_pvalue(
    p.col = "Pvalue",
    output.col = "Padj",
    method = "bonferroni"
  )


HeatmapData$signif <- HeatmapData$Padj <= 0.05

HeatmapData$signif[HeatmapData$Signif == FALSE] <- NA_real_

HeatmapData <- HeatmapData %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

write.csv(HeatmapData, "Rhododendron_PGLS_DoseXChemDiv_72hr.csv")



HeatmapData$y.label <- factor(HeatmapData$y.label, levels = unique(HeatmapData$y.label))


megaPlot_72 <- ggplot(data = HeatmapData, aes(x = X.label, y = fct_rev(y.label), fill = Beta_plot))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )

dev.off()

ggsave("Final_RhodoTox_AllData_ChromophoreClass_72.svg", megaPlot_72, height = 15, width = 10)



DoseBetas72 <- HeatmapData %>%
  filter(X.label %in% c("Slope", "pEC50"))



plots_72 <- ggplot(DoseBetas72,aes(x = X.label, y = fct_rev(y.label), fill = Beta_plot)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_fixed()

plots_72

DoseBetas72_sigonly <- DoseBetas72 %>%
  filter(Padj <= 0.05)


plots_72_sigonly <- ggplot(DoseBetas72_sigonly,aes(x = X.label, y = fct_rev(y.label), fill = Beta)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue"
                       #na.value = "ivory4"
  )+
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_fixed()


# Final plot for main text
Dose_Plots <- ggarrange(plots_24, plots_72, ncol = 2, nrow = 1)

ggsave("RhodoTox_DoseVar_ByChromophore.svg", Dose_Plots, height = 10, width = 15)



################################################################################
## PCA Combining Ecophysiological Data and Chromophore Data
## and observing effects on dose-response parameters
################################################################################


EcoData <- read.csv(ecophys_path)

EcoMeans <- EcoData[, -c(1:3)]

EcoMeans <- EcoMeans %>%
  group_by(Species) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

qual <- qualpal(as.numeric(length(EcoMeans$Species)))
EcoPCA <- prcomp(EcoMeans[, -1], center = TRUE, scale. = TRUE)

Eco_PCAPlot <- fviz_pca_biplot(EcoPCA)

Eco_PCAPlot <- fviz_pca_var(
  EcoPCA,
  axes = c(1, 2),
  label = "var",
  alpha.var = "cos2",
  col.var = "cos2",
  # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE     # Avoid text overlapping
)+
  labs(title = NULL)



#produce plot of individuals
Eco_PCA.IND <- fviz_pca_ind(
  EcoPCA,
  geom.ind = "point",
  pointshape = 19,
  pointsize = 3,
  label = EcoMeans$Species,
  col.ind = EcoMeans$Species,
  mean.point = F
)+
  scale_color_manual(values = qual$hex)+
  labs(title = NULL)

EcoPlots <- ggarrange(Eco_PCAPlot, Eco_PCA.IND, nrow = 2, ncol = 1, labels = c("A)", "B)"))

ggsave("PCA_EcoPlots.svg", EcoPlots, height = 15, width = 10)




EcoChem <- merge(macro_chemicals, EcoMeans, 1)

EcoChemPCA <- prcomp(EcoChem[, -1], center = TRUE, scale. = TRUE)

EcoChem_PCAPlot <- fviz_pca_biplot(EcoChemPCA)



EcoChem_PCAPlot <- fviz_pca_var(
  EcoChemPCA,
  axes = c(1, 2),
  label = "var",
  alpha.var = "cos2",
  col.var = "cos2",
  # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE     # Avoid text overlapping
)+
  labs(title = NULL)



#produce plot of individuals
EcoChem_PCA.IND <- fviz_pca_ind(
  EcoChemPCA,
  geom.ind = "point",
  pointshape = 19,
  pointsize = 3,
  label = EcoMeans$Species,
  col.ind = EcoMeans$Species,
  mean.point = F
)+
  scale_color_manual(values = qual$hex)+
  labs(title = NULL)

EcoChemPlots <- ggarrange(EcoChem_PCAPlot, EcoChem_PCA.IND, nrow = 2, ncol = 1, labels = c("A)", "B)"))

ggsave("PCA_EcoChemPlots.svg", EcoChemPlots, height = 15, width = 10)

ChemPCA <- prcomp(macro_chemicals[, -1], center = TRUE, scale. = TRUE)

Chem_PCAPLot <- fviz_pca_biplot(ChemPCA)

Chem_PCAPLot<- fviz_pca_var(
  ChemPCA,
  axes = c(1, 2),
  label = "var",
  alpha.var = "cos2",
  col.var = "cos2",
  # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE     # Avoid text overlapping
)+
  labs(title = NULL)



#produce plot of individuals
Chem_PCA.IND <- fviz_pca_ind(
  ChemPCA,
  geom.ind = "point",
  pointshape = 19,
  pointsize = 3,
  label = EcoMeans$Species,
  col.ind = EcoMeans$Species,
  mean.point = F
)+
  scale_color_manual(values = qual$hex)+
  labs(title = NULL)

ChemPlots <- ggarrange(Chem_PCAPLot, Chem_PCA.IND, nrow = 2, ncol = 1, labels = c("A)", "B)"))

ggsave("PCA_ChemPlots.svg", ChemPlots, height = 15, width = 10)


PCAS <- ggarrange(Chem_PCAPLot, EcoChem_PCAPlot, Eco_PCAPlot, ncol = 2, nrow = 2)
ggsave("RhodoTox_PCAs_Alternative.svg", PCAS, height = 15, width = 15)

plot(ChemPCA$x[, 1], EcoPCA$x[, 1])
abline(line)

EcoChemData <- data.frame(
  Species = EcoMeans$Species,
  ChemPC1 = ChemPCA$x[, 1],
  EcoPC1 = EcoPCA$x[, 1],
  ChemEcoPC1 = EcoChemPCA$x[, 1],
  ChemPC2 = ChemPCA$x[, 2],
  EcoPC2 = EcoPCA$x[, 2],
  ChemEcoPC2 = EcoChemPCA$x[, 2]
)



trait_cols <- colnames(EcoChemData)[-1]

for (i in 1:(length(trait_cols))) {
  for (j in 1:length(trait_cols)) {
    
    trait1 <- trait_cols[i]
    trait2 <- trait_cols[j]
    
    # Subset data safely
    if (!all(c(trait1, trait2, "Species") %in% colnames(EcoChemData))) {
      message("Skipping missing columns: ", trait1, ", ", trait2)
      next
    }
    
    df <- EcoChemData[, c(trait1, trait2, "Species")]
    df <- na.omit(df)
    
    # ---- Scale traits ----
    df[[trait1]] <- df[[trait1]]
    df[[trait2]] <- df[[trait2]]
    
    if (nrow(df) < 3) next
    
    # ---- Build formula ----
    pgls_formula <- as.formula(paste(trait1, "~", trait2))
    
    # ---- Fit model ----
    fit <- tryCatch({
      gls(
        pgls_formula,
        data = df,
        correlation = corBrownian(1, phy = tree, form = ~Species),
        method = "ML"
      )
    }, error = function(e) {
      message("Failed: ", trait1, " vs ", trait2)
      return(NULL)
    })
    
    # ---- Extract results ----
    if (!is.null(fit)) {
      summ <- summary(fit)$tTable
      
      # Extract slope (trait2 row)
      if (trait2 %in% rownames(summ)) {
        beta <- summ[trait2, "Value"]
        pval <- summ[trait2, "p-value"]
        
        results_list[[counter]] <- data.frame(
          Trait1 = trait1,
          Trait2 = trait2,
          Beta = beta,
          P_value = pval,
          stringsAsFactors = FALSE
        )
        
        counter <- counter + 1
      }
    }
  }
}

# Combine into a single dataframe
results_df <- do.call(rbind, results_list)

results_df$adjusted_p <- p.adjust(results_df$P_value, method = "BH")

results_df$signif <- results_df$adjusted_p <= 0.05

write.csv(results_df, "Rhododendron_PGLS_PCAChemvChem.csv")



## PCA and Dose Data for 24 hours
DoseData24 <- DoseData %>%
  filter(Time == 24) %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

DoseData24_NoMean <- DoseData %>%
  filter(Time == 24)

DoseData24$Species <- gsub("\\. ", ".", DoseData24$Species)

Tox_Chem <- merge(DoseData24, EcoChemData, 1)

Tox_Parameters <- c("Slope", "c", "d", "pEC50", "var_Slope", "var_pEC50")

PCA_Parameters <- c("ChemPC1", "EcoPC1", "ChemPC2", "EcoPC2", "ChemEcoPC1", "ChemEcoPC2")

phylo_tree <- read.tree(tree_path)


Final_data <- list()

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- Tox_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric((df[[Tox_Var]]))
    df[[PCA_Var]] <- as.numeric((df[[PCA_Var]]))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model using the standardized data
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~Species
      ),
      method = "REML"
    )
    # Extract coefficient table
    coef_df <- as.data.frame(
      summary(Model)$tTable[-1, c(
        "Value",
        "Std.Error",
        "p-value"
      ), drop = FALSE]
    )
    
    coef_df$Variable <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add identifying information
    coef_df$Tox_Parameter <- Tox_Var
    coef_df$PCA_Parameter <- PCA_Var
    coef_df$Model <- deparse(Formula)
    
    # Store results
    Final_data[[length(Final_data) + 1]] <- coef_df
  }
}

# Combine all results
Final_data <- do.call(rbind, Final_data)
rownames(Final_data) <- NULL


HeatmapData <- data.frame(
  X.label = Final_data$Tox_Parameter,
  y.label = Final_data$PCA_Parameter,
  Beta = Final_data$Value,
  Pvalue = Final_data$`p-value`,
  Model = Final_data$Model
)
library(rstatix)


HeatmapData <- HeatmapData %>%
  adjust_pvalue(
    p.col = "Pvalue",
    output.col = "Padj",
    method = "bonferroni"
  )


HeatmapData$signif <- HeatmapData$Padj <= 0.05

HeatmapData$signif[HeatmapData$Signif == FALSE] <- NA_real_

HeatmapData <- HeatmapData %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

write.csv(HeatmapData, "Rhododendron_DosevPCA_24hr.csv")


megaPlot <- ggplot(data = HeatmapData, aes(x = X.label, y = y.label, fill = Beta_plot))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )
ggsave("AllDoseTraits_Unscaled_24hrs.svg", megaPlot, height = 10, width =10)

Heatmap_Subset <- HeatmapData %>%
  filter(X.label %in% c("Slope", "pEC50")) %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

heatmap <- ggplot(
  data = Heatmap_Subset,
  aes(x = y.label, y = X.label, fill = Beta_plot)
) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "darkred",
    mid = "white",
    high = "dodgerblue",
    na.value = "ivory4"
  ) +
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )

ggsave("Unscaled_onlyEC50andSlope_72hrs.svg", heatmap, height = 5, width = 10)


## PGLS regressions of all pc axes against all dose response parameters

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- Tox_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric((df[[Tox_Var]]))
    df[[PCA_Var]] <- as.numeric((df[[PCA_Var]]))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model using the standardized data
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~Species
      ),
      method = "REML"
    )
    
    
    svg(paste0(Tox_Var, "_predictedBy_", PCA_Var, "24hr.svg"), height = 5, width = 5)  
    plot(df[[PCA_Var]], df[[Tox_Var]], 
         xlab = PCA_Var, ylab = Tox_Var)
    abline(Model, col = "blue")
    dev.off()  
    
  }
}


# Open a single PDF
pdf("Tox_PCA_plots_24.pdf", width = 12, height = 16)

# Arrange plots in a 4 x 6 grid
par(mfrow = c(6, 4),
    mar = c(4, 4, 2, 1))

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- Tox_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Convert variables to numeric
    df[[Tox_Var]] <- as.numeric(df[[Tox_Var]])
    df[[PCA_Var]] <- as.numeric(df[[PCA_Var]])
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit GLS model
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~Species
      ),
      method = "REML"
    )
    
    # Plot
    plot(
      df[[PCA_Var]],
      df[[Tox_Var]],
      xlab = PCA_Var,
      ylab = Tox_Var,
      main = paste(Tox_Var, "vs", PCA_Var)
    )
    
    # Add regression line
    abline(Model, col = "blue", lwd = 2)
  }
}

# Close the PDF
dev.off()






## PCA and Dose Data for 72 hours ##############################################

DoseData72 <- DoseData %>%
  filter(Time == 72) %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

DoseData72_NoMean <- DoseData %>%
  filter(Time == 72)

DoseData72$Species <- gsub("\\. ", ".", DoseData24$Species)

Tox_Chem <- merge(DoseData72, EcoChemData, 1)

Tox_Parameters <- c("Slope", "c", "d", "pEC50", "var_Slope", "var_pEC50")

PCA_Parameters <- c("ChemPC1", "EcoPC1", "ChemPC2", "EcoPC2", "ChemEcoPC1", "ChemEcoPC2")

phylo_tree <- read.tree(tree_path)


Final_data <- list()

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- Tox_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric((df[[Tox_Var]]))
    df[[PCA_Var]] <- as.numeric((df[[PCA_Var]]))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model using the standardized data
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~Species
      ),
      method = "REML"
    )
    # Extract coefficient table
    coef_df <- as.data.frame(
      summary(Model)$tTable[-1, c(
        "Value",
        "Std.Error",
        "p-value"
      ), drop = FALSE]
    )
    
    coef_df$Variable <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add identifying information
    coef_df$Tox_Parameter <- Tox_Var
    coef_df$PCA_Parameter <- PCA_Var
    coef_df$Model <- deparse(Formula)
    
    # Store results
    Final_data[[length(Final_data) + 1]] <- coef_df
  }
}

# Combine all results
Final_data <- do.call(rbind, Final_data)
rownames(Final_data) <- NULL


HeatmapData <- data.frame(
  X.label = Final_data$Tox_Parameter,
  y.label = Final_data$PCA_Parameter,
  Beta = Final_data$Value,
  Pvalue = Final_data$`p-value`,
  Model = Final_data$Model
)
library(rstatix)


HeatmapData <- HeatmapData %>%
  adjust_pvalue(
    p.col = "Pvalue",
    output.col = "Padj",
    method = "bonferroni"
  )


HeatmapData$signif <- HeatmapData$Padj <= 0.05

HeatmapData$signif[HeatmapData$Signif == FALSE] <- NA_real_

HeatmapData <- HeatmapData %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

write.csv(HeatmapData, "Rhododendron_DosevPCA_72hr.csv")


megaPlot <- ggplot(data = HeatmapData, aes(x = X.label, y = y.label, fill = Beta_plot))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )
ggsave("AllDoseTraits_Unscaled_72hrs.svg", megaPlot, height = 10, width =10)


Heatmap_Subset <- HeatmapData %>%
  filter(X.label %in% c("Slope", "pEC50")) %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

heatmap <- ggplot(
  data = Heatmap_Subset,
  aes(x = y.label, y = X.label, fill = Beta_plot)
) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "darkred",
    mid = "white",
    high = "dodgerblue",
    na.value = "ivory4"
  ) +
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )

ggsave("Unscaled_onlyEC50andSlope_72hrs.svg", heatmap, height = 5, width = 10)





## PGLS regressions of all pc axes against all dose response parameters

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- Tox_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric((df[[Tox_Var]]))
    df[[PCA_Var]] <- as.numeric((df[[PCA_Var]]))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model using the standardized data
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~Species
      ),
      method = "REML"
    )
    
    
    svg(paste0(Tox_Var, "_predictedBy_", PCA_Var, "72hr.svg"), height = 5, width = 5)  
    plot(df[[PCA_Var]], df[[Tox_Var]], 
         xlab = PCA_Var, ylab = Tox_Var)
    abline(Model, col = "blue")
    dev.off()  
    
  }
}


# Open a single PDF
pdf("Tox_PCA_plots_72.pdf", width = 12, height = 16)

# Arrange plots in a 4 x 6 grid
par(mfrow = c(6, 4),
    mar = c(4, 4, 2, 1))

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(PCA_Parameters)) {
    
    PCA_Var <- PCA_Parameters[j]
    
    # Create data frame containing the two variables
    df <- Tox_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Convert variables to numeric
    df[[Tox_Var]] <- as.numeric(df[[Tox_Var]])
    df[[PCA_Var]] <- as.numeric(df[[PCA_Var]])
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit GLS model
    Model <- gls(
      Formula,
      data = df,
      correlation = corBrownian(
        1,
        phylo_tree,
        ~Species
      ),
      method = "REML"
    )
    
    # Plot
    plot(
      df[[PCA_Var]],
      df[[Tox_Var]],
      xlab = PCA_Var,
      ylab = Tox_Var,
      main = paste(Tox_Var, "vs", PCA_Var)
    )
    
    # Add regression line
    abline(Model, col = "blue", lwd = 2)
  }
}

# Close the PDF
dev.off()



################################################################################
## Effects of Individual Signals on Dose Response Parameters
################################################################################

Peak_Parameters <- colnames(Metabolites[-1])

Tox_Chem <- merge(DoseData24, EcoChemData, 1)

ToxPeak_Chem <- merge(Tox_Chem, Metabolites, by = 1)

Final_data <- list()

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(Peak_Parameters)) {
    
    PCA_Var <- Peak_Parameters[j]
    
    # Create data frame containing the two variables
    df <- ToxPeak_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric(scale(df[[Tox_Var]]))
    df[[PCA_Var]] <- as.numeric(scale(df[[PCA_Var]]))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model
    Model <- tryCatch({
      gls(
        Formula,
        data = df,
        correlation = corBrownian(1, phy = tree, form = ~Species),
        method = "REML"
      )
    }, error = function(e) {
      message("Failed: ", Tox_Var, " vs ", PCA_Var,
              " | ", e$message)
      return(NULL)
    })
    
    # Skip failed models
    if (is.null(Model)) next
    
    # Extract coefficient table
    coef_df <- as.data.frame(
      summary(Model)$tTable[-1, c(
        "Value",
        "Std.Error",
        "p-value"
      ), drop = FALSE]
    )
    
    coef_df$Variable <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add identifying information
    coef_df$Tox_Parameter <- Tox_Var
    coef_df$PCA_Parameter <- PCA_Var
    coef_df$Model <- deparse(Formula)
    
    # Store results
    Final_data[[length(Final_data) + 1]] <- coef_df
  }
}

# Combine results
if (length(Final_data) > 0) {
  Final_data <- do.call(rbind, Final_data)
} else {
  Final_data <- data.frame()
}


rownames(Final_data) <- NULL

df_split <- Final_data %>% 
  separate_wider_delim(
    cols = PCA_Parameter, 
    delim = "_", 
    names = c("RT", "Chromophore")
  )



HeatmapData <- data.frame(
  X.label = df_split$Tox_Parameter,
  y.label = df_split$RT,
  chromophore = df_split$Chromophore,
  Beta = df_split$Value,
  Pvalue = df_split$`p-value`,
  Model = df_split$Model
)



HeatmapData <- HeatmapData %>%
  adjust_pvalue(
    p.col = "Pvalue",
    output.col = "Padj",
    method = "bonferroni"
  )


HeatmapData$signif <- HeatmapData$Padj <= 0.05

HeatmapData$signif[HeatmapData$Signif == FALSE] <- NA_real_

HeatmapData <- HeatmapData %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))


megaPlot <- ggplot(data = HeatmapData, aes(x = X.label, y = y.label, fill = Beta_plot))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )

Heatmap_Subset <- HeatmapData %>%
  filter(X.label %in% c("Slope", "pEC50")) %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

write.csv(Heatmap_Subset, "Rhododendron_PeakXDose_24hr.csv")


# Name of your grouping variable
group_var <- "chromophore"

# Split the data into separate data frames
heatmap_data <- split(Heatmap_Subset, Heatmap_Subset[[group_var]])

# Generate a heatmap for each group
heatmap_list <- map(
  names(heatmap_data),
  function(g) {
    
    df <- heatmap_data[[g]]
    
    ggplot(
      data = df,
      aes(x = X.label, y = y.label, fill = Beta_plot)
    ) +
      geom_tile(color = "black") +
      scale_fill_gradient2(
        low = "darkred",
        mid = "white",
        high = "dodgerblue",
        na.value = "ivory4"
      ) +
      theme_classic2(base_size = 15) +
      theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      ggtitle(g)
  }
)

# Combine all heatmaps into a single multipanel figure
heatmap <- wrap_plots(
  heatmap_list,
  ncol = 2
)




ggsave("MultiPanel_ChromophoreHeatmap_24hr.svg", heatmap, height = 20, width = 10)


Sig_only_peaks <- Heatmap_Subset %>%
  filter(Padj <= 0.05)

qual <- qualpal(length(unique(Sig_only_peaks$chromophore)))

# One chromophore value per y.label
chromophore_points <- Sig_only_peaks %>%
  distinct(y.label, chromophore)

# Make sure y-axis ordering is consistent
y_levels <- unique(Sig_only_peaks$y.label)

chromophore_points <- chromophore_points %>%
  mutate(
    y.label = factor(y.label, levels = y_levels)
  )

Sig_only_peaks <- Sig_only_peaks %>%
  mutate(
    y.label = factor(y.label, levels = y_levels)
  )

megaPlot <- ggplot(
  Sig_only_peaks,
  aes(x = X.label, y = y.label)
) +
  
  # Heatmap
  geom_tile(
    aes(fill = Beta_plot),
    color = "black"
  ) +
  
  scale_fill_gradient2(
    low = "darkred",
    mid = "white",
    high = "dodgerblue",
    midpoint = 0,
    na.value = "ivory4",
    name = "Beta"
  ) +
  
  # Chromophore circles
  geom_point(
    data = chromophore_points,
    aes(
      x = -0.65,
      y = y.label,
      color = chromophore
    ),
    size = 4,
    shape = 16,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = qual$hex,
    name = "Chromophore"
  ) +
  
  # Give the plot space to the left of the heatmap
  scale_x_discrete(
    expand = expansion(mult = c(0.15, 0.05))
  ) +
  
  # Allow points outside normal panel boundaries
  coord_cartesian(
    xlim = c(-1, length(unique(Sig_only_peaks$X.label)) + 0.5),
    clip = "off"
  ) +
  
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    
    # Extra space on left for chromophore circles
    plot.margin = margin(
      t = 5.5,
      r = 1,
      b = 5.5,
      l = 5
    )
  )

ggsave("Significant_only_Plot_bypeak24hr.svg", megaPlot, height = 10, width = 10)



## Signals described at 72 hours

Tox_Chem <- merge(DoseData72, EcoChemData, 1)

ToxPeak_Chem <- merge(Tox_Chem, Metabolites, by = 1)

Final_data <- list()

for (i in seq_along(Tox_Parameters)) {
  
  Tox_Var <- Tox_Parameters[i]
  
  for (j in seq_along(Peak_Parameters)) {
    
    PCA_Var <- Peak_Parameters[j]
    
    # Create data frame containing the two variables
    df <- ToxPeak_Chem[, c(Tox_Var, PCA_Var, "Species")]
    
    # Standardize both variables
    df[[Tox_Var]] <- as.numeric(scale(df[[Tox_Var]]))
    df[[PCA_Var]] <- as.numeric(scale(df[[PCA_Var]]))
    
    # Construct formula
    Formula <- reformulate(
      PCA_Var,
      response = Tox_Var
    )
    
    # Fit model
    Model <- tryCatch({
      gls(
        Formula,
        data = df,
        correlation = corBrownian(1, phy = tree, form = ~Species),
        method = "REML"
      )
    }, error = function(e) {
      message("Failed: ", Tox_Var, " vs ", PCA_Var,
              " | ", e$message)
      return(NULL)
    })
    
    # Skip failed models
    if (is.null(Model)) next
    
    # Extract coefficient table
    coef_df <- as.data.frame(
      summary(Model)$tTable[-1, c(
        "Value",
        "Std.Error",
        "p-value"
      ), drop = FALSE]
    )
    
    coef_df$Variable <- rownames(coef_df)
    rownames(coef_df) <- NULL
    
    # Add identifying information
    coef_df$Tox_Parameter <- Tox_Var
    coef_df$PCA_Parameter <- PCA_Var
    coef_df$Model <- deparse(Formula)
    
    # Store results
    Final_data[[length(Final_data) + 1]] <- coef_df
  }
}

# Combine results
if (length(Final_data) > 0) {
  Final_data <- do.call(rbind, Final_data)
} else {
  Final_data <- data.frame()
}


rownames(Final_data) <- NULL

df_split <- Final_data %>% 
  separate_wider_delim(
    cols = PCA_Parameter, 
    delim = "_", 
    names = c("RT", "Chromophore")
  )



HeatmapData <- data.frame(
  X.label = df_split$Tox_Parameter,
  y.label = df_split$RT,
  chromophore = df_split$Chromophore,
  Beta = df_split$Value,
  Pvalue = df_split$`p-value`,
  Model = df_split$Model
)



HeatmapData <- HeatmapData %>%
  adjust_pvalue(
    p.col = "Pvalue",
    output.col = "Padj",
    method = "bonferroni"
  )


HeatmapData$signif <- HeatmapData$Padj <= 0.05

HeatmapData$signif[HeatmapData$Signif == FALSE] <- NA_real_

HeatmapData <- HeatmapData %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))


megaPlot <- ggplot(data = HeatmapData, aes(x = X.label, y = y.label, fill = Beta_plot))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "dodgerblue",
                       na.value = "ivory4")+
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )

Heatmap_Subset <- HeatmapData %>%
  filter(X.label %in% c("Slope", "pEC50")) %>%
  mutate(Beta_plot = ifelse(signif == 0, NA, Beta))

write.csv(Heatmap_Subset, "Rhododendron_PeakXDose_72hr.csv")


# Name of your grouping variable
group_var <- "chromophore"

# Split the data into separate data frames
heatmap_data <- split(Heatmap_Subset, Heatmap_Subset[[group_var]])

# Generate a heatmap for each group
heatmap_list <- map(
  names(heatmap_data),
  function(g) {
    
    df <- heatmap_data[[g]]
    
    ggplot(
      data = df,
      aes(x = X.label, y = y.label, fill = Beta_plot)
    ) +
      geom_tile(color = "black") +
      scale_fill_gradient2(
        low = "darkred",
        mid = "white",
        high = "dodgerblue",
        na.value = "ivory4"
      ) +
      theme_classic2(base_size = 15) +
      theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      ggtitle(g)
  }
)

# Combine all heatmaps into a single multipanel figure
heatmap <- wrap_plots(
  heatmap_list,
  ncol = 2
)

heatmap



ggsave("MultiPanel_ChromophoreHeatmap_72hr.svg", heatmap, height = 20, width = 10)



Sig_only_peaks <- Heatmap_Subset %>%
  filter(Padj <= 0.05)

qual <- qualpal(length(unique(Sig_only_peaks$chromophore)))

# One chromophore value per y.label
chromophore_points <- Sig_only_peaks %>%
  distinct(y.label, chromophore)

# Make sure y-axis ordering is consistent
y_levels <- unique(Sig_only_peaks$y.label)

chromophore_points <- chromophore_points %>%
  mutate(
    y.label = factor(y.label, levels = y_levels)
  )

Sig_only_peaks <- Sig_only_peaks %>%
  mutate(
    y.label = factor(y.label, levels = y_levels)
  )

megaPlot <- ggplot(
  Sig_only_peaks,
  aes(x = X.label, y = y.label)
) +
  
  # Heatmap
  geom_tile(
    aes(fill = Beta_plot),
    color = "black"
  ) +
  
  scale_fill_gradient2(
    low = "darkred",
    mid = "white",
    high = "dodgerblue",
    midpoint = 0,
    na.value = "ivory4",
    name = "Beta"
  ) +
  
  # Chromophore circles
  geom_point(
    data = chromophore_points,
    aes(
      x = -0.65,
      y = y.label,
      color = chromophore
    ),
    size = 4,
    shape = 16,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = qual$hex,
    name = "Chromophore"
  ) +
  
  # Give the plot space to the left of the heatmap
  scale_x_discrete(
    expand = expansion(mult = c(0.15, 0.05))
  ) +
  
  # Allow points outside normal panel boundaries
  coord_cartesian(
    xlim = c(-1, length(unique(Sig_only_peaks$X.label)) + 0.5),
    clip = "off"
  ) +
  
  theme_classic2(base_size = 15) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    
    # Extra space on left for chromophore circles
    plot.margin = margin(
      t = 5.5,
      r = 1,
      b = 5.5,
      l = 5
    )
  )

ggsave("Significant_only_Plot_bypeak72hr.svg", megaPlot, height = 10, width = 10)
















