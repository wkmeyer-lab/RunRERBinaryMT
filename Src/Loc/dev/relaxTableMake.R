####################################################################################################################
####    Initial script written by Alexander Seaver                                                              ####
####    Select sections have been revised utilizing Google Gemini 3.6, manually reviewed by Alexander Seaver    ####
####################################################################################################################

# Load necessary libraries
library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)

# Define base path and target directories
base_dir <- "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output"

analysis_dirs <- c(
  "ComplexDietCentralAnalysisSimplify2",
  "ComplexDietCentralAnalysisSimplifyExclusivePred"
)

rate_folders <- c("Rate_2", "Rate_3", "Rate_4")

# Helper function to access list elements case-insensitively
get_case_insensitive_key <- function(lst, target_key) {
  if (is.null(lst)) return(NULL)
  nms <- names(lst)
  if (is.null(nms)) return(NULL)
  match_idx <- which(tolower(nms) == tolower(target_key))
  if (length(match_idx) > 0) {
    return(lst[[match_idx[1]]])
  }
  return(NULL)
}

# Helper function to safely extract single values that might be missing
safe_extract <- function(x) {
  if (!is.null(x) && length(x) > 0) return(x[[1]]) else return(NA)
}

# Helper function to extract up to 4-category omegas and site proportions from RELAX alternative
extract_rate_distribution <- function(json_data, branch_set) {
  res <- list(
    omega_0 = NA_real_, omega_1 = NA_real_, omega_2 = NA_real_, omega_3 = NA_real_,
    prop_0  = NA_real_, prop_1  = NA_real_, prop_2  = NA_real_, prop_3  = NA_real_
  )

  fits <- get_case_insensitive_key(json_data, "fits")
  alt_fit <- get_case_insensitive_key(fits, "RELAX alternative")
  rate_dists <- get_case_insensitive_key(alt_fit, "Rate Distributions")
  branch_dist <- get_case_insensitive_key(rate_dists, branch_set)

  if (is.null(branch_dist)) return(res)

  omegas <- NULL
  props  <- NULL

  # Case 1: Data frame (standard jsonlite default)
  if (is.data.frame(branch_dist)) {
    cols <- tolower(colnames(branch_dist))

    if ("omega" %in% cols) {
      omegas <- branch_dist[[which(cols == "omega")[1]]]
    } else if (ncol(branch_dist) >= 1) {
      omegas <- branch_dist[[1]]
    }

    if ("proportion" %in% cols) {
      props <- branch_dist[[which(cols == "proportion")[1]]]
    } else if ("p" %in% cols) {
      props <- branch_dist[[which(cols == "p")[1]]]
    } else if ("weight" %in% cols) {
      props <- branch_dist[[which(cols == "weight")[1]]]
    } else if (ncol(branch_dist) >= 2) {
      props <- branch_dist[[2]]
    }

  # Case 2: Matrix
  } else if (is.matrix(branch_dist)) {
    cols <- tolower(colnames(branch_dist))
    if (!is.null(cols) && "omega" %in% cols) {
      omegas <- branch_dist[, which(cols == "omega")[1]]
    } else if (ncol(branch_dist) >= 1) {
      omegas <- branch_dist[, 1]
    }

    if (!is.null(cols) && "proportion" %in% cols) {
      props <- branch_dist[, which(cols == "proportion")[1]]
    } else if (!is.null(cols) && "p" %in% cols) {
      props <- branch_dist[, which(cols == "p")[1]]
    } else if (!is.null(cols) && "weight" %in% cols) {
      props <- branch_dist[, which(cols == "weight")[1]]
    } else if (ncol(branch_dist) >= 2) {
      props <- branch_dist[, 2]
    }

  # Case 3: List / nested structure
  } else if (is.list(branch_dist)) {
    extracted_om <- c()
    extracted_pr <- c()
    for (i in seq_along(branch_dist)) {
      item <- branch_dist[[i]]
      val_om <- NULL
      val_pr <- NULL
      if (is.list(item) || is.data.frame(item)) {
        val_om <- safe_extract(get_case_insensitive_key(item, "omega"))
        val_pr <- safe_extract(get_case_insensitive_key(item, "proportion"))
        if (is.null(val_pr)) val_pr <- safe_extract(get_case_insensitive_key(item, "p"))
        if (is.null(val_pr)) val_pr <- safe_extract(get_case_insensitive_key(item, "weight"))
      } else if (is.numeric(item)) {
        if (length(item) >= 1) val_om <- item[1]
        if (length(item) >= 2) val_pr <- item[2]
      }
      if (!is.null(val_om)) extracted_om <- c(extracted_om, as.numeric(val_om))
      if (!is.null(val_pr)) extracted_pr <- c(extracted_pr, as.numeric(val_pr))
    }
    if (length(extracted_om) > 0) omegas <- extracted_om
    if (length(extracted_pr) > 0) props  <- extracted_pr
  }

  # Assign extracted Omegas
  if (!is.null(omegas) && length(omegas) > 0) {
    if (length(omegas) >= 1) res$omega_0 <- as.numeric(omegas[1])
    if (length(omegas) >= 2) res$omega_1 <- as.numeric(omegas[2])
    if (length(omegas) >= 3) res$omega_2 <- as.numeric(omegas[3])
    if (length(omegas) >= 4) res$omega_3 <- as.numeric(omegas[4])
  }

  # Assign extracted Proportions
  if (!is.null(props) && length(props) > 0) {
    if (length(props) >= 1) res$prop_0 <- as.numeric(props[1])
    if (length(props) >= 2) res$prop_1 <- as.numeric(props[2])
    if (length(props) >= 3) res$prop_2 <- as.numeric(props[3])
    if (length(props) >= 4) res$prop_3 <- as.numeric(props[4])
  }

  return(res)
}

# Initialize an empty list to collect data rows efficiently
all_rows <- list()

# Nested loops across outer directories and rate subdirectories
for (dir_name in analysis_dirs) {
  for (folder in rate_folders) {

    folder_path <- file.path(base_dir, dir_name, "Hyphy", folder)

    # Skip missing directories gracefully
    if (!dir.exists(folder_path)) {
      warning(paste("Directory not found, skipping:", folder_path))
      next
    }

    json_files <- list.files(path = folder_path, pattern = "\\.json$", full.names = TRUE)

    for (file_path in json_files) {
      file_name <- basename(file_path)

      # Extract Gene and Foreground number from filename
      gene_name     <- str_match(file_name, "relax-(.*?)-Foreground")[2]
      foreground_id <- str_match(file_name, "Foreground_(\\d+)")[2]

      # Parse JSON with safety check
      json_data <- tryCatch({
        fromJSON(file_path)
      }, error = function(e) {
        warning(paste("Failed to read or parse JSON:", file_name))
        return(NULL)
      })

      if (is.null(json_data)) next

      # Extract statistical parameters
      test_results <- get_case_insensitive_key(json_data, "test results")
      k_param      <- safe_extract(get_case_insensitive_key(test_results, "relaxation or intensification parameter"))
      p_val        <- safe_extract(get_case_insensitive_key(test_results, "p-value"))
      lrt          <- safe_extract(get_case_insensitive_key(test_results, "LRT"))

      # Extract rate distributions (Omegas for Ref & Test, shared Proportions from Ref)
      dist_ref  <- extract_rate_distribution(json_data, "Reference")
      dist_test <- extract_rate_distribution(json_data, "Test")

      # Extract AICc values
      fits     <- get_case_insensitive_key(json_data, "fits")
      gen_fit  <- get_case_insensitive_key(fits, "General descriptive")
      alt_fit  <- get_case_insensitive_key(fits, "RELAX alternative")
      nul_fit  <- get_case_insensitive_key(fits, "RELAX null")
      par_fit  <- get_case_insensitive_key(fits, "RELAX partitioned descriptive")

      aicc_gen <- safe_extract(get_case_insensitive_key(gen_fit, "AIC-c"))
      aicc_alt <- safe_extract(get_case_insensitive_key(alt_fit, "AIC-c"))
      aicc_nul <- safe_extract(get_case_insensitive_key(nul_fit, "AIC-c"))
      aicc_par <- safe_extract(get_case_insensitive_key(par_fit, "AIC-c"))

      # Build individual row dataframe
      row_data <- data.frame(
        Analysis_Directory       = dir_name,
        Gene                      = gene_name,
        Rate_Folder               = folder,
        Foreground                = foreground_id,
        K_Parameter               = as.numeric(k_param),
        P_Value                   = as.numeric(p_val),
        Likelihood_Ratio          = as.numeric(lrt),

        # Omegas
        Omega_Ref_0               = dist_ref$omega_0,
        Omega_Ref_1               = dist_ref$omega_1,
        Omega_Ref_2               = dist_ref$omega_2,
        Omega_Ref_3               = dist_ref$omega_3,
        Omega_Test_0              = dist_test$omega_0,
        Omega_Test_1              = dist_test$omega_1,
        Omega_Test_2              = dist_test$omega_2,
        Omega_Test_3              = dist_test$omega_3,

        # Global Rate Class Site Proportions
        Prop_0                    = dist_ref$prop_0,
        Prop_1                    = dist_ref$prop_1,
        Prop_2                    = dist_ref$prop_2,
        Prop_3                    = dist_ref$prop_3,

        AICC_General_Descriptive = as.numeric(aicc_gen),
        AICC_RELAX_Alternative   = as.numeric(aicc_alt),
        AICC_RELAX_Null          = as.numeric(aicc_nul),
        AICC_RELAX_Partitioned   = as.numeric(aicc_par),
        stringsAsFactors          = FALSE
      )

      all_rows[[length(all_rows) + 1]] <- row_data
    }
  }
}

# Combine all rows into a single dataframe
all_results <- bind_rows(all_rows)

# Replace raw p-values equal to 0 with 1e-16 prior to correction
all_results <- all_results %>%
  mutate(P_Value = if_else(P_Value == 0, 1e-16, P_Value))

# Flag the row with the lowest AICc for RELAX Alternative per Analysis, Gene, and Foreground group
all_results <- all_results %>%
  group_by(Analysis_Directory, Gene, Foreground) %>%
  mutate(
    Is_Lowest_AICc_Alt = !is.na(AICC_RELAX_Alternative) & 
                          AICC_RELAX_Alternative == min(AICC_RELAX_Alternative, na.rm = TRUE)
  ) %>%
  ungroup()

# Initialize P_Value_Adj_BH as NA for all rows
all_results$P_Value_Adj_BH <- NA_real_

# Calculate Benjamini-Hochberg FDR correction strictly on the lowest AICc subset
lowest_idx <- which(all_results$Is_Lowest_AICc_Alt)
all_results$P_Value_Adj_BH[lowest_idx] <- p.adjust(all_results$P_Value[lowest_idx], method = "BH")

# Reorder columns logically
all_results <- all_results %>%
  select(
    Analysis_Directory, Gene, Rate_Folder, Foreground,
    K_Parameter, P_Value, P_Value_Adj_BH, Likelihood_Ratio,
    Is_Lowest_AICc_Alt,
    Omega_Ref_0, Omega_Ref_1, Omega_Ref_2, Omega_Ref_3,
    Omega_Test_0, Omega_Test_1, Omega_Test_2, Omega_Test_3,
    Prop_0, Prop_1, Prop_2, Prop_3,
    everything()
  )

# Sort the final dataframe logically
all_results <- all_results %>%
  arrange(Analysis_Directory, Gene, as.numeric(Foreground), Rate_Folder)

# Filter dataframe for best rate class per Analysis, Gene, and Foreground combination
lowest_aicc_results <- all_results %>%
  filter(Is_Lowest_AICc_Alt)

# Save combined results to CSV (na = "" ensures non-lowest AICc rows are written as blank cells)
output_csv_combined <- file.path(base_dir, "RELAX_Extracted_Results_Combined.csv")
write.csv(all_results, file = output_csv_combined, row.names = FALSE, na = "")

# Save filtered best-fitting rate class results to separate CSV
output_csv_best <- file.path(base_dir, "RELAX_Extracted_Results_Lowest_AICc.csv")
write.csv(lowest_aicc_results, file = output_csv_best, row.names = FALSE, na = "")

message("Extraction complete!")
message("  - Combined results saved (", nrow(all_results), " rows) to: ", output_csv_combined)
message("  - Lowest AICc results saved (", nrow(lowest_aicc_results), " rows) to: ", output_csv_best)
