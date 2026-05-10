get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

print_usage <- function(script_name = "run_abundance_prediction.R") {
  cat(
    paste(
      "Usage:",
      sprintf("  Rscript %s --species \"Anas crecca\"", script_name),
      sprintf("  Rscript %s --species \"Anas crecca,Ardea alba\" --n-iter 5", script_name),
      sprintf("  Rscript %s --species-file species.txt --output-folder .", script_name),
      sprintf("  Rscript %s --all-species --checklist-folder \"D:/path/to/checklists\"", script_name),
      sprintf("  Rscript %s --list-species", script_name),
      "",
      "Options:",
      "  --species                   Scientific name(s), comma-separated or repeated.",
      "  --species-file              TXT/CSV file containing scientific names.",
      "  --all-species               Run all species in the checklist folder sequentially.",
      "  --list-species              List available species and exit.",
      "  --eu-shp-path               Path to the EU grid shapefile.",
      "  --env-folder                Folder containing yearly environmental CSV files.",
      "  --lc-path                   Path to the land cover CSV file.",
      "  --checklist-folder          Folder containing species checklist CSV files.",
      "  --output-folder             Root output folder. Default: script folder.",
      "  --start-year                Start year for filtering. Default: 2021.",
      "  --end-year                  End year for filtering. Default: 2022.",
      "  --protocol                  Protocol(s), comma-separated. Default: Traveling.",
      "  --obs-quantile-cutoff       Observation count cutoff. Default: 0.99.",
      "  --observer-quantile-cutoff  Observer cutoff. Default: 0.7.",
      "  --grid-quantile-cutoff      Grid cutoff. Default: 0.5.",
      "  --validation-ratio          Validation split ratio. Default: 0.1.",
      "  --n-iter                    Number of model iterations. Default: 3.",
      "  --seed                      Random seed. Default: 123.",
      "  --help                      Show this help message.",
      sep = "\n"
    )
  )
}

trim_ws <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  script_dir <- get_script_dir()
  opts <- list(
    species = character(),
    species_file = NULL,
    all_species = FALSE,
    list_species = FALSE,
    help = FALSE,
    eu_shp_path = file.path(script_dir, "EU_100km_fishnet_simple_by_distance", "EU_100km_fishnet_simple_by_distance.shp"),
    env_folder = file.path(script_dir, "gee_data", "era5_2016_2022"),
    lc_path = file.path(script_dir, "gee_data", "land_cover_2016_2022.csv"),
    checklist_folder = file.path(script_dir, "ebird_filtered_checklist"),
    output_folder = script_dir,
    start_year = 2021,
    end_year = 2022,
    protocol = "Traveling",
    obs_quantile_cutoff = 0.99,
    observer_quantile_cutoff = 0.7,
    grid_quantile_cutoff = 0.5,
    validation_ratio = 0.1,
    n_iter = 3,
    seed = 123
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]

    if (arg %in% c("--help", "-h")) {
      opts$help <- TRUE
      i <- i + 1
      next
    }
    if (arg == "--all-species") {
      opts$all_species <- TRUE
      i <- i + 1
      next
    }
    if (arg == "--list-species") {
      opts$list_species <- TRUE
      i <- i + 1
      next
    }

    if (!startsWith(arg, "--")) {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
    if (i == length(args) || startsWith(args[i + 1], "--")) {
      stop("Missing value for argument: ", arg, call. = FALSE)
    }

    value <- args[i + 1]
    key <- sub("^--", "", arg)

    if (key == "species") {
      split_species <- unlist(strsplit(value, ",", fixed = TRUE))
      split_species <- trim_ws(split_species)
      split_species <- split_species[nzchar(split_species)]
      opts$species <- c(opts$species, split_species)
    } else if (key == "species-file") {
      opts$species_file <- value
    } else if (key == "eu-shp-path") {
      opts$eu_shp_path <- value
    } else if (key == "env-folder") {
      opts$env_folder <- value
    } else if (key == "lc-path") {
      opts$lc_path <- value
    } else if (key == "checklist-folder") {
      opts$checklist_folder <- value
    } else if (key == "output-folder") {
      opts$output_folder <- value
    } else if (key == "start-year") {
      opts$start_year <- as.integer(value)
    } else if (key == "end-year") {
      opts$end_year <- as.integer(value)
    } else if (key == "protocol") {
      opts$protocol <- value
    } else if (key == "obs-quantile-cutoff") {
      opts$obs_quantile_cutoff <- as.numeric(value)
    } else if (key == "observer-quantile-cutoff") {
      opts$observer_quantile_cutoff <- as.numeric(value)
    } else if (key == "grid-quantile-cutoff") {
      opts$grid_quantile_cutoff <- as.numeric(value)
    } else if (key == "validation-ratio") {
      opts$validation_ratio <- as.numeric(value)
    } else if (key == "n-iter") {
      opts$n_iter <- as.integer(value)
    } else if (key == "seed") {
      opts$seed <- as.integer(value)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }

    i <- i + 2
  }

  opts$checklist_folder <- normalizePath(opts$checklist_folder, winslash = "/", mustWork = FALSE)
  opts$output_folder <- normalizePath(opts$output_folder, winslash = "/", mustWork = FALSE)
  opts$eu_shp_path <- normalizePath(opts$eu_shp_path, winslash = "/", mustWork = FALSE)
  opts$env_folder <- normalizePath(opts$env_folder, winslash = "/", mustWork = FALSE)
  opts$lc_path <- normalizePath(opts$lc_path, winslash = "/", mustWork = FALSE)
  opts$protocols <- trim_ws(unlist(strsplit(opts$protocol, ",", fixed = TRUE)))
  opts$protocols <- opts$protocols[nzchar(opts$protocols)]

  opts
}

extract_species_name_from_file <- function(filename) {
  name <- sub("\\.csv$", "", filename, ignore.case = TRUE)
  name <- sub("^\\^", "", name)
  name <- sub("_filtered_2019to2022$", "", name)
  name <- sub("_filtered2019to2022$", "", name)
  trim_ws(name)
}

build_species_catalog <- function(checklist_folder) {
  if (!dir.exists(checklist_folder)) {
    stop("Checklist folder not found: ", checklist_folder, call. = FALSE)
  }

  csv_paths <- list.files(
    path = checklist_folder,
    pattern = "\\.csv$",
    full.names = TRUE
  )

  if (length(csv_paths) == 0) {
    stop("No CSV files found in checklist folder: ", checklist_folder, call. = FALSE)
  }

  catalog <- data.frame(
    scientific_name = vapply(basename(csv_paths), extract_species_name_from_file, character(1)),
    source_file = basename(csv_paths),
    source_path = normalizePath(csv_paths, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  catalog[order(catalog$scientific_name), , drop = FALSE]
}

read_species_file <- function(path) {
  if (!file.exists(path)) {
    stop("Species file not found: ", path, call. = FALSE)
  }

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if ("scientific_name" %in% names(df)) {
      species <- df$scientific_name
    } else {
      species <- df[[1]]
    }
  } else {
    species <- readLines(path, warn = FALSE)
  }

  species <- trim_ws(species)
  unique(species[nzchar(species)])
}

resolve_species_paths <- function(opts, species_catalog) {
  if (opts$all_species) {
    return(species_catalog)
  }

  requested <- opts$species
  if (!is.null(opts$species_file)) {
    requested <- c(requested, read_species_file(opts$species_file))
  }
  requested <- unique(trim_ws(requested))
  requested <- requested[nzchar(requested)]

  if (length(requested) == 0) {
    stop(
      "Please provide at least one of --species, --species-file, --all-species, or --list-species.",
      call. = FALSE
    )
  }

  matched <- species_catalog[species_catalog$scientific_name %in% requested, , drop = FALSE]
  missing_species <- setdiff(requested, matched$scientific_name)
  if (length(missing_species) > 0) {
    stop(
      "Species not found in checklist folder: ",
      paste(missing_species, collapse = ", "),
      call. = FALSE
    )
  }

  matched[match(requested, matched$scientific_name), , drop = FALSE]
}

validate_options <- function(opts) {
  if (opts$start_year > opts$end_year) {
    stop("--start-year must be less than or equal to --end-year.", call. = FALSE)
  }
  if (length(opts$protocols) == 0) {
    stop("At least one protocol is required.", call. = FALSE)
  }
  if (!dir.exists(opts$checklist_folder)) {
    stop("Checklist folder not found: ", opts$checklist_folder, call. = FALSE)
  }
  if (!dir.exists(opts$env_folder)) {
    stop("Environment folder not found: ", opts$env_folder, call. = FALSE)
  }
  if (!file.exists(opts$lc_path)) {
    stop("Land cover file not found: ", opts$lc_path, call. = FALSE)
  }
  if (!file.exists(opts$eu_shp_path)) {
    stop("EU shapefile not found: ", opts$eu_shp_path, call. = FALSE)
  }
  if (opts$validation_ratio <= 0 || opts$validation_ratio >= 1) {
    stop("--validation-ratio must be between 0 and 1.", call. = FALSE)
  }
  for (cutoff_name in c("obs_quantile_cutoff", "observer_quantile_cutoff", "grid_quantile_cutoff")) {
    cutoff <- opts[[cutoff_name]]
    if (is.na(cutoff) || cutoff <= 0 || cutoff >= 1) {
      stop(sprintf("--%s must be between 0 and 1.", gsub("_", "-", cutoff_name)), call. = FALSE)
    }
  }
  if (is.na(opts$n_iter) || opts$n_iter < 1) {
    stop("--n-iter must be at least 1.", call. = FALSE)
  }
}

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_name <- if (length(script_file_arg) > 0) {
  basename(sub("^--file=", "", script_file_arg[1]))
} else {
  "run_abundance_prediction.R"
}

opts <- parse_args()

if (opts$help) {
  print_usage(script_name)
  quit(save = "no", status = 0)
}

species_catalog <- build_species_catalog(opts$checklist_folder)

if (opts$list_species) {
  write.table(
    species_catalog["scientific_name"],
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  quit(save = "no", status = 0)
}

validate_options(opts)
selected_species <- resolve_species_paths(opts, species_catalog)

suppressPackageStartupMessages({
  library(mgcv)
  library(parallel)
  library(pbapply)
  library(pdp)
  library(viridis)
  library(rnaturalearth)
  library(sf)
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(xgboost)
  library(caret)
  library(lubridate)
  library(tidyr)
  library(zoo)
  library(MASS)
  library(glmmTMB)
  library(car)
  library(DHARMa)
  library(data.table)
  library(tibble)
  library(grid)
  library(gridExtra)
  library(png)
})

year_range_num <- seq(opts$start_year, opts$end_year)
year_range_chr <- as.character(year_range_num)

CCI_output_folder <- file.path(opts$output_folder, "glmm_performance")
validtion_data_folder <- file.path(opts$output_folder, "validation_data")
abundance_model_output_folder <- file.path(opts$output_folder, "abundance_spatiotemporal_sampling_method")

dir.create(opts$output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(CCI_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(validtion_data_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(abundance_model_output_folder, recursive = TRUE, showWarnings = FALSE)

get_EU_map_and_neighbors <- function(shp_path) {
  EU.map <- st_read(shp_path, quiet = TRUE)
  EU_map <- as.data.frame(EU.map)
  neighbors <- st_touches(EU.map)

  matrix_data <- matrix(nrow = length(EU.map$Id), ncol = 9)
  colnames(matrix_data) <- c("id", paste0("neighbor", 1:8))
  matrix_data[, 1] <- EU.map$Id

  for (i in seq_along(EU.map$Id)) {
    neighbors_indices <- neighbors[[i]]
    if (length(neighbors_indices) > 0) {
      neighbors_id <- EU.map[neighbors_indices, ]
      for (j in seq_len(min(length(neighbors_id$Id), 8))) {
        matrix_data[i, j + 1] <- neighbors_id$Id[j]
      }
    }
  }

  neighbors_df <- as.data.frame(matrix_data)
  list(EU.map = EU.map, EU_map = EU_map, neighbors_df = neighbors_df)
}

load_environmental_data <- function(start_year, end_year, base_path) {
  en_df <- NULL

  for (year_num in start_year:end_year) {
    path <- file.path(base_path, paste0(year_num, "_median_combined_result.csv"))
    part_df <- read.csv(path)
    if (is.null(en_df)) {
      en_df <- part_df
    } else {
      en_df <- rbind(en_df, part_df)
    }
  }

  en_df <- en_df %>%
    separate(Month, into = c("year_number", "month_number"), sep = "-", convert = TRUE)

  en_df$month_number <- as.factor(en_df$month_number)
  en_df$year_number <- as.factor(en_df$year_number)
  en_df
}

load_land_cover_data <- function(file_path, year_range) {
  land_cover_df <- read.csv(file_path)

  land_cover_df <- land_cover_df %>%
    arrange(Id, year_number, month_number)

  land_cover_df$month_number <- as.factor(land_cover_df$month_number)
  land_cover_df$year_number <- as.factor(land_cover_df$year_number)

  land_cover_df %>%
    filter(year_number %in% year_range)
}

process_observation_data <- function(data_Id, quantile_threshold = 0.99, year_range, protocol_range) {
  checklist_feature_col <- c(
    "observation_count", "locality_type", "observation_date", "observer_id",
    "protocol_type", "duration_mins", "effort_km", "Id", "month_number"
  )
  data_Id <- data_Id[, checklist_feature_col]

  data_Id$effort_km[is.na(data_Id$effort_km)] <- 0
  data_Id$duration_mins[is.na(data_Id$duration_mins)] <- 0

  data_Id <- data_Id %>%
    mutate(year_number = year(as.Date(observation_date))) %>%
    filter(year_number %in% year_range)

  X_count <- nrow(data_Id[which(data_Id$observation_count == "X"), ])
  X_proportion <- X_count / nrow(data_Id)

  data_Id$observation_count <- ifelse(data_Id$observation_count == "X", "0", data_Id$observation_count)
  data_Id$observation_count <- as.numeric(data_Id$observation_count)

  data_Id <- data_Id %>%
    filter(protocol_type %in% protocol_range)

  data_Id <- merge(data_Id, en_df, by = c("Id", "year_number", "month_number"), all = FALSE)
  data_Id <- merge(data_Id, lc_df, by = c("Id", "year_number", "month_number"), all = FALSE)

  ids_with_na <- data_Id %>%
    group_by(Id) %>%
    filter(any(is.na(across(everything())))) %>%
    pull(Id) %>%
    unique()
  data_Id <- data_Id[!data_Id$Id %in% ids_with_na, ]

  obs_summary <- summary(data_Id$observation_count)
  num_row_data <- nrow(data_Id)

  ob_count_threshold <- quantile(data_Id$observation_count, probs = quantile_threshold, na.rm = TRUE)
  print(paste0(quantile_threshold * 100, "% quantile = ", ob_count_threshold))
  data_Id <- data_Id %>%
    filter(observation_count <= ob_count_threshold)

  list(
    X_proportion = X_proportion,
    observation_summary = obs_summary,
    ob_count_threshold = ob_count_threshold,
    num_row_data = num_row_data,
    filtered_data = data_Id
  )
}

filter_observers_by_quantile <- function(data_Id, X = 0.7) {
  observer_freq <- data_Id %>%
    count(observer_id, name = "observation_times")

  probs <- seq(0.5, 0.95, by = 0.05)

  retained_stats <- lapply(probs, function(p) {
    threshold <- quantile(observer_freq$observation_times, probs = p)
    top_observers <- observer_freq %>%
      filter(observation_times >= threshold) %>%
      pull(observer_id)
    retained_n <- data_Id %>%
      filter(observer_id %in% top_observers) %>%
      nrow()
    data.frame(
      Percentile = paste0(p * 100, "%"),
      Threshold = as.numeric(threshold),
      Retained_N = retained_n,
      Retained_Ratio = retained_n / nrow(data_Id)
    )
  }) %>% bind_rows()

  retained_plot <- ggplot(retained_stats, aes(x = Percentile, y = Retained_Ratio)) +
    geom_line(group = 1, color = "steelblue", linewidth = 1.2) +
    geom_point(size = 3, color = "darkred") +
    geom_text(aes(label = paste0(round(Retained_Ratio * 100, 1), "%")), vjust = -0.5, size = 4.5) +
    labs(
      title = "sample retention ratio vs observer observation number quantile threshold",
      x = "quantile threshold of observer observation times",
      y = "sample retention ratio"
    ) +
    ylim(0, 1) +
    theme_minimal()

  threshold_X <- quantile(observer_freq$observation_times, probs = X)
  top_observers_X <- observer_freq %>%
    filter(observation_times >= threshold_X) %>%
    pull(observer_id)

  data_top_X <- data_Id %>%
    filter(observer_id %in% top_observers_X)

  retained_ratio <- nrow(data_top_X) / nrow(data_Id)

  list(
    plot = retained_plot,
    threshold = threshold_X,
    retained_ratio = retained_ratio,
    data_top_X = data_top_X
  )
}

filter_grids_with_plot <- function(data_Id, quantile_threshold = 0.5) {
  id_freq <- data_Id %>%
    count(Id, name = "observation_times")

  probs <- seq(0.3, 0.95, by = 0.05)

  retained_stats <- lapply(probs, function(p) {
    threshold <- quantile(id_freq$observation_times, probs = p)
    retained_ids <- id_freq$Id[id_freq$observation_times >= threshold]
    retained_n <- data_Id %>%
      filter(Id %in% retained_ids) %>%
      nrow()
    data.frame(
      Percentile = paste0(p * 100, "%"),
      Threshold = as.numeric(threshold),
      Retained_N = retained_n,
      Retained_Ratio = retained_n / nrow(data_Id)
    )
  }) %>% bind_rows()

  retained_plot <- ggplot(retained_stats, aes(x = Percentile, y = Retained_Ratio)) +
    geom_line(group = 1, color = "steelblue", linewidth = 1.2) +
    geom_point(size = 3, color = "darkred") +
    geom_text(aes(label = paste0(round(Retained_Ratio * 100, 1), "%")), vjust = -0.5, size = 4.5) +
    labs(
      title = "Sample retention ratio vs quantile threshold of grid observation times",
      x = "quantile threshold of grid observation times",
      y = "sample retention ratio"
    ) +
    ylim(0, 1) +
    theme_minimal()

  threshold <- quantile(id_freq$observation_times, probs = quantile_threshold)
  retained_ids <- id_freq$Id[id_freq$observation_times >= threshold]

  filtered_data <- data_Id %>%
    filter(Id %in% retained_ids) %>%
    left_join(id_freq, by = "Id") %>%
    rename(id_observation_times = observation_times)

  num_grids_above <- length(retained_ids)
  total_grids <- nrow(id_freq)
  proportion_above <- num_grids_above / total_grids

  cat("Quantile threshold =", quantile_threshold, "\n")
  cat("Observation count threshold =", threshold, "\n")
  cat("The number of grids that meet the threshold =", num_grids_above, "/", total_grids, "\n")
  cat("The ratio =", round(proportion_above, 4), "(", round(proportion_above * 100, 2), "%)\n")

  list(
    threshold = threshold,
    valid_grid_count = num_grids_above,
    total_grids = total_grids,
    proportion = proportion_above,
    filtered_data = filtered_data,
    plot = retained_plot
  )
}

calculate_observer_cci <- function(data_Id, birdname, cci_output_folder) {
  cci_randomeffect_output_folder <- file.path(cci_output_folder, "cci_random_effect")
  glmm_performance_output_folder <- file.path(cci_output_folder, "glmm_performance")

  if (!dir.exists(cci_randomeffect_output_folder)) {
    dir.create(cci_randomeffect_output_folder, recursive = TRUE)
  }
  if (!dir.exists(glmm_performance_output_folder)) {
    dir.create(glmm_performance_output_folder, recursive = TRUE)
  }

  data_Id <- data_Id %>%
    mutate(
      duration_mins_log = log(duration_mins + 1),
      effort_km_log = log(effort_km + 1),
      observer_id = as.factor(observer_id)
    )

  if (all(c(1, 0) %in% unique(data_Id$protocol_type))) {
    model_poisson <- glmmTMB(
      observation_count ~ duration_mins_log + (1 | observer_id),
      family = poisson,
      data = data_Id
    )

    model_nbinom <- glmmTMB(
      observation_count ~ duration_mins_log + (1 | observer_id),
      family = nbinom2,
      data = data_Id
    )
  } else if (0 %in% unique(data_Id$protocol_type)) {
    model_poisson <- glmmTMB(
      observation_count ~ duration_mins_log + (1 | observer_id),
      family = poisson,
      data = data_Id
    )

    model_nbinom <- glmmTMB(
      observation_count ~ duration_mins_log + (1 | observer_id),
      family = nbinom2,
      data = data_Id
    )
  } else {
    model_poisson <- glmmTMB(
      observation_count ~ duration_mins_log + effort_km_log + (1 | observer_id),
      family = poisson,
      data = data_Id
    )

    model_nbinom <- glmmTMB(
      observation_count ~ duration_mins_log + effort_km_log + (1 | observer_id),
      family = nbinom2,
      data = data_Id
    )
  }

  cat("poisson model:\n")
  print(summary(model_poisson))
  cat("nbinom model:\n")
  print(summary(model_nbinom))

  cat("AIC:\n")
  aic_values <- AIC(model_poisson, model_nbinom)
  cat("BIC:\n")
  bic_values <- BIC(model_poisson, model_nbinom)

  cat("=== AIC performance ===\n")
  print(aic_values)
  cat("=== BIC performance===\n")
  print(bic_values)

  glmm_performance_df <- data.frame(
    birdname = birdname,
    AIC_poisson = aic_values$AIC[1],
    AIC_nbinom = aic_values$AIC[2],
    BIC_poisson = bic_values$BIC[1],
    BIC_nbinom = bic_values$BIC[2]
  )

  write.csv(
    glmm_performance_df,
    file.path(glmm_performance_output_folder, paste0(birdname, ".csv")),
    row.names = FALSE
  )

  random_effects <- ranef(model_nbinom)
  observer_effects <- random_effects$cond$observer_id

  observer_effects_df <- data.frame(
    observer_id = rownames(observer_effects),
    cci = scale(observer_effects[, "(Intercept)"]),
    row.names = NULL
  )

  write.csv(
    observer_effects_df,
    file = file.path(cci_randomeffect_output_folder, paste0(birdname, ".csv")),
    row.names = FALSE
  )

  cat("=== Observer CCI summary ===\n")
  cat("Median:", median(observer_effects_df$cci), "\n")
  cat("Mean:", mean(observer_effects_df$cci), "\n")
  cat("SD:", sd(observer_effects_df$cci), "\n")

  data_Id$observer_id <- as.character(data_Id$observer_id)
  merge(data_Id, observer_effects_df, by = "observer_id", all.x = TRUE)
}

calculate_poisson_deviance_explained <- function(y_true, y_pred) {
  y_pred_safe <- pmax(y_pred, 1e-15)
  y_mean <- mean(y_true)
  y_mean_safe <- pmax(y_mean, 1e-15)

  D_model <- 2 * sum(
    ifelse(
      y_true == 0,
      y_pred_safe,
      y_true * log(y_true / y_pred_safe) - (y_true - y_pred_safe)
    ),
    na.rm = TRUE
  )

  D_null <- 2 * sum(
    ifelse(
      y_true == 0,
      y_mean_safe,
      y_true * log(y_true / y_mean_safe) - (y_true - y_mean_safe)
    ),
    na.rm = TRUE
  )

  1 - (D_model / D_null)
}

spatiotemporal_sampling <- function(
  data_Id, validation_data, unique_data_Id, en_df, lc_df,
  n_iter, seed = 123, bird_name, output_folder, prediction_years
) {
  set.seed(seed)

  train_rmse_list <- numeric(n_iter)
  oob_rmse_list <- numeric(n_iter)
  train_size_list <- numeric(n_iter)
  oob_size_list <- numeric(n_iter)
  train_spearman_list <- numeric(n_iter)
  oob_spearman_list <- numeric(n_iter)
  train_pde_list <- numeric(n_iter)
  oob_pde_list <- numeric(n_iter)

  abundance_output_folder <- file.path(output_folder, "abundance_prediction", bird_name)
  performance_output_folder <- file.path(output_folder, "model_performance")
  validation_output_folder <- file.path(output_folder, "validation_prediction", bird_name)
  importance_output_folder <- file.path(output_folder, "feature_importance", bird_name)

  dir.create(abundance_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(performance_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(validation_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(importance_output_folder, recursive = TRUE, showWarnings = FALSE)

  for (j in 1:12) {
    validation_data[[paste0("month_number", j)]] <- ifelse(validation_data$month_number == j, 1, 0)
  }
  for (j in prediction_years) {
    validation_data[[paste0("year_number", j)]] <- ifelse(as.integer(as.character(validation_data$year_number)) == j, 1, 0)
  }

  for (i in seq_len(n_iter)) {
    cat("\nIteration", i, "\n")

    zero_proportion <- sum(data_Id$observation_count == 0) / nrow(data_Id)
    sample_size <- ceiling(mean(table(data_Id$year_number)))
    sample_Id_size <- ceiling(sample_size / length(table(data_Id$Id)))

    Id_list <- table(data_Id$Id)
    sampled_data <- data_Id %>%
      filter(Id %in% names(Id_list)) %>%
      group_by(Id, year_number) %>%
      group_modify(~ {
        if (nrow(.x) == 0) {
          return(tibble())
        }
        slice_sample(.x, n = sample_Id_size, replace = TRUE)
      }) %>%
      ungroup()

    sample_zero_count <- ceiling(length(sampled_data$observation_count) * zero_proportion)
    sample_nonzero_count <- length(sampled_data$observation_count) - sample_zero_count

    zero_sample <- sampled_data %>%
      filter(observation_count == 0) %>%
      slice_sample(n = sample_zero_count, replace = TRUE)

    nonzero_sample <- sampled_data %>%
      filter(observation_count > 0) %>%
      slice_sample(n = sample_nonzero_count, replace = TRUE)

    sampled_data <- rbind(zero_sample, nonzero_sample)

    for (j in 1:12) {
      sampled_data[[paste0("month_number", j)]] <- ifelse(sampled_data$month_number == j, 1, 0)
    }
    for (j in prediction_years) {
      sampled_data[[paste0("year_number", j)]] <- ifelse(as.integer(as.character(sampled_data$year_number)) == j, 1, 0)
    }

    if (all(c(0, 1) %in% unique(data_Id$protocol_type))) {
      exclude_cols <- c(
        "observation_date", "observation_count", "Id", "observer_id", "binary_count",
        "year_number", "month_number", "duration_mins_log", "effort_km_log", "sample_key"
      )
    } else {
      exclude_cols <- c(
        "observation_date", "observation_count", "Id", "observer_id", "binary_count",
        "year_number", "month_number", "duration_mins_log", "effort_km_log", "sample_key", "protocol_type"
      )
    }

    feature_cols <- setdiff(colnames(data_Id), exclude_cols)
    filtered_data <- as.data.frame(sampled_data)[, feature_cols]
    X_train <- model.matrix(~ . - 1, data = filtered_data)
    y_train <- sampled_data$observation_count
    dtrain <- xgb.DMatrix(data = X_train, label = y_train)

    model <- xgboost(
      data = dtrain, objective = "count:poisson", eta = 0.3, nrounds = 500,
      eval_metric = "rmse", colsample_bytree = 0.8, nthread = 12, verbose = 0
    )

    importance_matrix <- xgb.importance(model = model)
    write.csv(
      importance_matrix,
      file.path(importance_output_folder, paste0("iteration", i, ".csv")),
      row.names = FALSE
    )

    full_data <- data_Id
    for (j in 1:12) {
      full_data[[paste0("month_number", j)]] <- ifelse(full_data$month_number == j, 1, 0)
    }
    for (j in prediction_years) {
      full_data[[paste0("year_number", j)]] <- ifelse(as.integer(as.character(full_data$year_number)) == j, 1, 0)
    }

    full_X <- model.matrix(~ . - 1, data = as.data.frame(full_data)[, feature_cols])
    pred <- predict(model, full_X)

    dvalid <- model.matrix(~ . - 1, data = as.data.frame(validation_data)[, feature_cols])
    valid_pred <- predict(model, dvalid)
    valid_residual <- valid_pred - validation_data$observation_count

    valid_pred_df <- setNames(
      data.frame(
        Id = validation_data$Id,
        year_number = validation_data$year_number,
        month_number = validation_data$month_number,
        sample_key = validation_data$sample_key,
        valid_pred
      ),
      c("Id", "year_number", "month_number", "sample_key", paste0("abundance", i))
    )
    write.csv(
      valid_pred_df,
      file.path(validation_output_folder, sprintf("abundance_iteration%d.csv", i)),
      row.names = FALSE
    )

    residual <- pred - full_data$observation_count
    full_data$set <- ifelse(full_data$sample_key %in% sampled_data$sample_key, "train", "oob")

    train_rmse <- sqrt(mean(residual[full_data$set == "train"]^2, na.rm = TRUE))
    oob_rmse <- sqrt(mean(residual[full_data$set == "oob"]^2, na.rm = TRUE))
    valid_rmse <- sqrt(mean(valid_residual^2, na.rm = TRUE))

    train_spearman <- cor(
      pred[full_data$set == "train"],
      full_data$observation_count[full_data$set == "train"],
      method = "spearman"
    )
    oob_spearman <- cor(
      pred[full_data$set == "oob"],
      full_data$observation_count[full_data$set == "oob"],
      method = "spearman"
    )
    valid_spearman <- cor(
      valid_pred,
      validation_data$observation_count,
      method = "spearman"
    )

    train_pde_list[i] <- calculate_poisson_deviance_explained(
      full_data$observation_count[full_data$set == "train"],
      pred[full_data$set == "train"]
    )

    oob_pde_list[i] <- calculate_poisson_deviance_explained(
      full_data$observation_count[full_data$set == "oob"],
      pred[full_data$set == "oob"]
    )
    valid_pde <- calculate_poisson_deviance_explained(
      validation_data$observation_count,
      valid_pred
    )

    cat(
      "Train RMSE:", round(train_rmse, 4), " | OOB RMSE:", round(oob_rmse, 4),
      " | validation RMSE:", round(valid_rmse, 4), "\n"
    )
    cat(
      "Train size:", sum(full_data$set == "train"), " | OOB size:", sum(full_data$set == "oob"),
      " | validation size:", nrow(validation_data), "\n"
    )
    cat(
      "Train Spearman:", round(train_spearman, 4), " | OOB Spearman:", round(oob_spearman, 4),
      " | validation Spearman:", round(valid_spearman, 4), "\n"
    )
    cat(
      "Train P-DE:", round(train_pde_list[i], 4), " | OOB P-DE:", round(oob_pde_list[i], 4),
      " | validation P-DE:", round(valid_pde, 4), "\n"
    )

    train_rmse_list[i] <- train_rmse
    oob_rmse_list[i] <- oob_rmse
    train_size_list[i] <- sum(full_data$set == "train")
    oob_size_list[i] <- sum(full_data$set == "oob")
    train_spearman_list[i] <- train_spearman
    oob_spearman_list[i] <- oob_spearman

    Id_list1 <- unique(unique_data_Id$Id)
    abundance_sublist <- vector("list", length(Id_list1))

    for (k in seq_along(Id_list1)) {
      gid <- Id_list1[k]
      gid1 <- unique_data_Id$id_observation_times[unique_data_Id$Id == gid][1]

      predict_data <- expand.grid(
        Id = gid,
        year_number = prediction_years,
        month_number = 1:12,
        id_observation_times = gid1
      )
      predict_data$locality_type <- 1
      predict_data$protocol_type <- 1
      predict_data$duration_mins <- 60
      predict_data$effort_km <- 2
      predict_data$cci <- 0.5

      for (j in 1:12) {
        predict_data[[paste0("month_number", j)]] <- ifelse(predict_data$month_number == j, 1, 0)
      }
      for (j in prediction_years) {
        predict_data[[paste0("year_number", j)]] <- ifelse(predict_data$year_number == j, 1, 0)
      }

      predict_data$year_number <- as.factor(predict_data$year_number)
      predict_data$month_number <- as.factor(predict_data$month_number)

      predict_data <- merge(predict_data, en_df, by = c("Id", "year_number", "month_number"), all = FALSE)
      predict_data <- merge(predict_data, lc_df, by = c("Id", "year_number", "month_number"), all = FALSE)

      predict_data_X <- predict_data[, intersect(feature_cols, colnames(predict_data))]
      pred_features <- model.matrix(~ . - 1, data = predict_data_X)
      dpred <- xgb.DMatrix(data = pred_features)
      predictions_poisson <- predict(model, dpred)

      result_df <- data.frame(
        Id = predict_data$Id,
        year_number = predict_data$year_number,
        month_number = predict_data$month_number,
        abundance = predictions_poisson
      )
      abundance_sublist[[k]] <- result_df
    }

    write.csv(
      bind_rows(abundance_sublist),
      file.path(abundance_output_folder, sprintf("abundance_iteration%d.csv", i)),
      row.names = FALSE
    )
  }

  result_df <- data.frame(
    iteration = seq_along(train_rmse_list),
    train_rmse = train_rmse_list,
    oob_rmse = oob_rmse_list,
    train_size = train_size_list,
    oob_size = oob_size_list,
    train_spearman_list = train_spearman_list,
    oob_spearman_list = oob_spearman_list,
    train_pde_list = train_pde_list,
    oob_pde_list = oob_pde_list
  )

  write.csv(
    result_df,
    file.path(performance_output_folder, paste0(bird_name, ".csv")),
    row.names = FALSE
  )
  result_df
}

EU_map_info_list <- get_EU_map_and_neighbors(opts$eu_shp_path)
EU.map <- EU_map_info_list$EU.map
EU_map <- EU_map_info_list$EU_map
neighbors_df <- EU_map_info_list$neighbors_df

en_df_raw <- load_environmental_data(
  start_year = opts$start_year,
  end_year = opts$end_year,
  base_path = opts$env_folder
)

lc_df_raw <- load_land_cover_data(
  file_path = opts$lc_path,
  year_range = year_range_chr
)

for (row_idx in seq_len(nrow(selected_species))) {
  path_csv <- selected_species$source_path[row_idx]

  print("Load species checklist...")
  oringinal_data_Id <- read.csv(path_csv)

  bird_name <- selected_species$scientific_name[row_idx]
  print(paste0("species name : ", bird_name))

  en_df <- en_df_raw
  lc_df <- lc_df_raw

  env_vars <- colnames(en_df)[grepl("_MEDIAN$", colnames(en_df))]
  lc_vars <- colnames(lc_df)[grepl("_MEAN$", colnames(lc_df))]
  en_df[env_vars] <- scale(en_df[env_vars])
  lc_df[lc_vars] <- scale(lc_df[lc_vars])

  print("Data filtering step1 ...")
  filtered_result <- process_observation_data(
    oringinal_data_Id,
    quantile_threshold = opts$obs_quantile_cutoff,
    year_range = year_range_num,
    protocol_range = opts$protocols
  )
  num_row_data <- filtered_result$num_row_data
  data_Id <- filtered_result$filtered_data
  print(table(data_Id$protocol_type))

  print("Data filtering step2 ...")
  filtered_result <- filter_observers_by_quantile(
    data_Id,
    X = opts$observer_quantile_cutoff
  )
  data_Id <- filtered_result$data_top_X
  print(filtered_result$plot)

  print("Data filtering step3 ...")
  filtered_result <- filter_grids_with_plot(
    data_Id,
    quantile_threshold = opts$grid_quantile_cutoff
  )
  data_Id <- filtered_result$filtered_data
  print(filtered_result$plot)

  print(
    paste0(
      "The ratio of the number of filtered samples to the number of original samples = ",
      nrow(data_Id) / num_row_data
    )
  )
  print("Data filtering completed")

  if (nrow(data_Id) < 10) {
    warning("Skipping species due to insufficient records after filtering: ", bird_name)
    next
  }

  data_Id$locality_type <- ifelse(data_Id$locality_type == "H", 1, 0)
  data_Id$binary_count <- ifelse(data_Id$observation_count > 0, 1, 0)
  data_Id$protocol_type <- ifelse(data_Id$protocol_type == "Traveling", 1, 0)
  cat("protocol_type:\n")
  print(table(data_Id$protocol_type))
  cat("locality_type:\n")
  print(table(data_Id$locality_type))
  cat("observation count binary format:\n")
  print(table(data_Id$binary_count))

  zero_proportion <- sum(data_Id$observation_count == 0) / length(data_Id$observation_count)
  print(paste0("sample zero proportion = ", zero_proportion))

  print("Calculate the observation bias correction index...")
  data_Id <- calculate_observer_cci(
    data_Id,
    birdname = bird_name,
    cci_output_folder = CCI_output_folder
  )

  data_Id_only_effort <- data_Id[, c(
    "observer_id", "Id", "year_number", "month_number", "observation_count", "locality_type",
    "protocol_type", "duration_mins", "effort_km", "id_observation_times", "cci"
  )]

  max_cci_id <- data_Id_only_effort %>%
    filter(cci == max(cci, na.rm = TRUE)) %>%
    pull(observer_id) %>%
    unique()

  min_cci_id <- data_Id_only_effort %>%
    filter(cci == min(cci, na.rm = TRUE)) %>%
    pull(observer_id) %>%
    unique()

  cat("observer_id with highest CCI:", max_cci_id, "\n")
  cat("observer_id with lowest CCI:", min_cci_id, "\n")

  print("Split training set and testing set...")
  data_Id$sample_key <- seq_len(nrow(data_Id))
  unique_data_Id <- data_Id[, c("Id", "year_number", "month_number", "id_observation_times")]
  unique_data_Id <- unique_data_Id %>%
    distinct()
  print(paste0("number of pair (grid, time) : ", nrow(unique_data_Id)))

  set.seed(opts$seed)
  n <- nrow(data_Id)
  val_size <- floor(opts$validation_ratio * n)
  if (val_size < 1) {
    stop("Validation split is too small. Increase data size or validation ratio.", call. = FALSE)
  }

  val_indices <- sample(seq_len(n), size = val_size)
  validation_data <- data_Id[val_indices, ]
  data_Id <- data_Id[-val_indices, ]
  validation_observation_count_data <- validation_data[, c("sample_key", "observation_count")]

  write.csv(
    validation_observation_count_data,
    file.path(validtion_data_folder, paste0(bird_name, ".csv")),
    row.names = FALSE
  )

  print("start training the model and predicting results...")
  result <- spatiotemporal_sampling(
    data_Id = data_Id,
    validation_data = validation_data,
    unique_data_Id = unique_data_Id,
    en_df = en_df,
    lc_df = lc_df,
    n_iter = opts$n_iter,
    seed = opts$seed,
    bird_name = bird_name,
    output_folder = abundance_model_output_folder,
    prediction_years = year_range_num
  )
}
