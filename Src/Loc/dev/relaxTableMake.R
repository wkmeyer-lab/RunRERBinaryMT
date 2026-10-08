# Load necessary libraries
library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)

# Define base path and target directories
base_dir <- "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output"

analysis_dirs <- c(
  "ComplexDietCentralAnalysisSimplify2",
  "ComplexDietCentralAnalysisSimplifyEqualDrop",
  "ComplexDietCentralAnalysisSimplifyStrictPred"
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

# Helper function to extract up to 4-category omega distributions for Reference or Test
extract_omega_values <- function(json_data, branch_set) {
  res <- list(
    omega_0 = NA_real_,
    omega_1 = NA_real_,
    omega_2 = NA_real_,
    omega_3 = NA_real_
  )
  
  fits <- get_case_insensitive_key(json_data, "fits")
  alt_fit <- get_case_insensitive_key(fits, "RELAX alternative")
  rate_dists <- get_case_insensitive_key(alt_fit, "Rate Distributions")
  branch_dist <- get_case_insensitive_key(rate_dists, branch_set)
  
  if (is.null(branch_dist)) return(res)
  
  omegas <- NULL
  
  # Case 1: Data frame (standard jsonlite default)
  if (is.data.frame(branch_dist)) {
    if ("omega" %in% colnames(branch_dist)) {
      omegas <- branch_dist[["omega"]]
    } else if (ncol(branch_dist) >= 1) {
      omegas <- branch_dist[[1]]
    }
  # Case 2: Matrix
  } else if (is.matrix(branch_dist)) {
    if ("omega" %in% colnames(branch_dist)) {
      omegas <- branch_dist[, "omega"]
    } else if (ncol(branch_dist) >= 1) {
      omegas <- branch_dist[, 1]
    }
  # Case 3: List / nested structure
  } else if (is.list(branch_dist)) {
    extracted <- c()
    for (i in seq_along(branch_dist)) {
      item <- branch_dist[[i]]
      val <- NULL
      if (is.list(item) || is.data.frame(item)) {
        val <- safe_extract(item[["omega"]])
      } else if (is.numeric(item)) {
        val <- item
      }
      if (!is.null(val)) extracted <- c(extracted, as.numeric(val))
    }
    if (length(extracted) > 0) omegas <- extracted
  }
  
  if (!is.null(omegas) && length(omegas) > 0) {
    if (length(omegas) >= 1) res$omega_0 <- as.numeric(omegas[1])
    if (length(omegas) >= 2) res$omega_1 <- as.numeric(omegas[2])
    if (length(omegas) >= 3) res$omega_2 <- as.numeric(omegas[3])
    if (length(omegas) >= 4) res$omega_3 <- as.numeric(omegas[4])
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

      # Extract rate distributions (omega values for Reference and Test)
      omega_ref  <- extract_omega_values(json_data, "Reference")
      omega_test <- extract_omega_values(json_data, "Test")

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
        Gene                     = gene_name,
        Rate_Folder              = folder,
        Foreground               = foreground_id,
        K_Parameter              = as.numeric(k_param),
        P_Value                  = as.numeric(p_val),
        Likelihood_Ratio         = as.numeric(lrt),
        Omega_Ref_0              = omega_ref$omega_0,
        Omega_Ref_1              = omega_ref$omega_1,
        Omega_Ref_2              = omega_ref$omega_2,
        Omega_Ref_3              = omega_ref$omega_3,
        Omega_Test_0             = omega_test$omega_0,
        Omega_Test_1             = omega_test$omega_1,
        Omega_Test_2             = omega_test$omega_2,
        Omega_Test_3             = omega_test$omega_3,
        AICC_General_Descriptive = as.numeric(aicc_gen),
        AICC_RELAX_Alternative   = as.numeric(aicc_alt),
        AICC_RELAX_Null          = as.numeric(aicc_nul),
        AICC_RELAX_Partitioned   = as.numeric(aicc_par),
        stringsAsFactors         = FALSE
      )

      all_rows[[length(all_rows) + 1]] <- row_data
    }
  }
}

# Combine all rows into a single dataframe
all_results <- bind_rows(all_rows)

# Multiple hypothesis testing correction (Benjamini-Hochberg FDR)
all_results <- all_results %>%
  mutate(P_Value_Adj_BH = p.adjust(P_Value, method = "BH")) %>%
  select(
    Analysis_Directory, Gene, Rate_Folder, Foreground,
    K_Parameter, P_Value, P_Value_Adj_BH, Likelihood_Ratio,
    Omega_Ref_0, Omega_Ref_1, Omega_Ref_2, Omega_Ref_3,
    Omega_Test_0, Omega_Test_1, Omega_Test_2, Omega_Test_3,
    everything()
  )

# Sort the final dataframe logically
all_results <- all_results %>%
  arrange(Analysis_Directory, Gene, as.numeric(Foreground), Rate_Folder)

# Save results to CSV
output_csv <- file.path(base_dir, "RELAX_Extracted_Results_Combined.csv")
write.csv(all_results, file = output_csv, row.names = FALSE)

message("Extraction complete! Saved ", nrow(all_results), " rows to: ", output_csv)
