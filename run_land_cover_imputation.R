required_packages <- c("sf", "dplyr", "zoo", "xgboost", "parallel", "pbapply", "data.table")

load_required_packages <- function() {
  for (pkg in required_packages) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

all_land_cover_types <- c(
  "bare", "built", "crops", "flooded_vegetation", "grass",
  "shrub_and_scrub", "snow_and_ice", "trees", "water"
)

print_help <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript run_land_cover_imputation.R [options]",
      "",
      "Options:",
      "  --mode <sampling|imputation|aggregate|all>   Execution mode (default: all)",
      "  --eu-shp-path <path>                         European fishnet shapefile path",
      "  --input-csv <path>                          Land cover + climate CSV with missing values",
      "  --output-folder <path>                      Output root folder (default: project folder)",
      "  --seed-start <int>                          First seed number (default: 123)",
      "  --seed-end <int>                            Last seed number (default: 222)",
      "  --n-cores <int>                             Number of cores for parallel steps (default: 2)",
      "  --land-cover-types <csv>                    Comma-separated land cover types (default: all 9 types)",
      "  --help                                      Show this message",
      "",
      "Examples:",
      "  Rscript run_land_cover_imputation.R --mode all --seed-start 123 --seed-end 130 --n-cores 2",
      "  Rscript run_land_cover_imputation.R --mode aggregate",
      sep = "\n"
    )
  )
}

stop_with_message <- function(message_text) {
  stop(message_text, call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || is.na(x) || identical(x, "")) y else x
}

get_script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  normalizePath(dirname(sub("^--file=", "", script_arg[1])), winslash = "/", mustWork = TRUE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    mode = "all",
    eu_shp_path = file.path(
      get_script_dir(),
      "EU_100km_fishnet_simple_by_distance",
      "EU_100km_fishnet_simple_by_distance.shp"
    ),
    input_csv = file.path(
      get_script_dir(),
      "gee_data",
      "EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv"
    ),
    output_folder = get_script_dir(),
    seed_start = 123L,
    seed_end = 222L,
    n_cores = 2L,
    land_cover_types = all_land_cover_types
  )

  if ("--help" %in% args) {
    print_help()
    quit(save = "no", status = 0)
  }

  parsed <- defaults
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--mode") {
      i <- i + 1L
      parsed$mode <- args[[i]]
    } else if (arg == "--eu-shp-path") {
      i <- i + 1L
      parsed$eu_shp_path <- args[[i]]
    } else if (arg == "--input-csv") {
      i <- i + 1L
      parsed$input_csv <- args[[i]]
    } else if (arg == "--output-folder") {
      i <- i + 1L
      parsed$output_folder <- args[[i]]
    } else if (arg == "--seed-start") {
      i <- i + 1L
      parsed$seed_start <- as.integer(args[[i]])
    } else if (arg == "--seed-end") {
      i <- i + 1L
      parsed$seed_end <- as.integer(args[[i]])
    } else if (arg == "--n-cores") {
      i <- i + 1L
      parsed$n_cores <- as.integer(args[[i]])
    } else if (arg == "--land-cover-types") {
      i <- i + 1L
      parsed$land_cover_types <- trimws(strsplit(args[[i]], ",", fixed = TRUE)[[1]])
    } else {
      stop_with_message(paste("Unknown argument:", arg))
    }
    i <- i + 1L
  }

  valid_modes <- c("sampling", "imputation", "aggregate", "all")
  if (!parsed$mode %in% valid_modes) {
    stop_with_message(paste("Invalid --mode. Choose one of:", paste(valid_modes, collapse = ", ")))
  }
  if (is.na(parsed$seed_start) || is.na(parsed$seed_end) || parsed$seed_start > parsed$seed_end) {
    stop_with_message("--seed-start and --seed-end must be valid integers with seed_start <= seed_end.")
  }
  if (is.na(parsed$n_cores) || parsed$n_cores < 1L) {
    stop_with_message("--n-cores must be a positive integer.")
  }
  if (length(parsed$land_cover_types) == 0 || any(parsed$land_cover_types == "")) {
    stop_with_message("--land-cover-types must contain at least one valid type.")
  }
  invalid_land_cover_types <- setdiff(parsed$land_cover_types, all_land_cover_types)
  if (length(invalid_land_cover_types) > 0) {
    stop_with_message(
      paste(
        "Invalid --land-cover-types:",
        paste(invalid_land_cover_types, collapse = ", "),
        "| valid values are:",
        paste(all_land_cover_types, collapse = ", ")
      )
    )
  }

  parsed$eu_shp_path <- normalizePath(parsed$eu_shp_path, winslash = "/", mustWork = FALSE)
  parsed$input_csv <- normalizePath(parsed$input_csv, winslash = "/", mustWork = FALSE)
  parsed$output_folder <- normalizePath(parsed$output_folder, winslash = "/", mustWork = FALSE)
  parsed$seed_list <- seq.int(parsed$seed_start, parsed$seed_end)
  parsed$land_cover_types <- unique(parsed$land_cover_types)
  parsed
}

build_paths <- function(output_folder) {
  land_cover_root <- file.path(output_folder, "land_cover_imputation")
  list(
    root = land_cover_root,
    sampling_out_folder = file.path(land_cover_root, "sampling_data"),
    ml_imputation_output_folder = file.path(land_cover_root, "ml_prediction_output"),
    two_method_performance_output_folder = file.path(land_cover_root, "two_method_performance"),
    two_method_test_imputation_output_folder = file.path(land_cover_root, "two_method_test_output"),
    final_output_folder = file.path(land_cover_root, "final_output")
  )
}

get_final_output_path <- function(paths, selected_land_cover_types) {
  if (length(selected_land_cover_types) == 1) {
    lc_type <- selected_land_cover_types[[1]]
    lc_output_folder <- file.path(paths$final_output_folder, lc_type)
    if (!dir.exists(lc_output_folder)) {
      dir.create(lc_output_folder, recursive = TRUE, showWarnings = FALSE)
    }
    return(file.path(
      lc_output_folder,
      paste0("EU_2016_2022_land_cover_imputation_by_xgboost_", lc_type, ".csv")
    ))
  }
  file.path(
    paths$final_output_folder,
    "EU_2016_2022_land_cover_imputation_by_xgboost.csv"
  )
}

ensure_output_dirs <- function(paths) {
  dirs <- c(
    paths$root,
    paths$sampling_out_folder,
    paths$ml_imputation_output_folder,
    paths$two_method_performance_output_folder,
    paths$two_method_test_imputation_output_folder,
    paths$final_output_folder
  )
  for (dir_path in dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  for (lc_type in all_land_cover_types) {
    lc_dir <- file.path(paths$ml_imputation_output_folder, lc_type)
    if (!dir.exists(lc_dir)) {
      dir.create(lc_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

validate_input_paths <- function(config) {
  if (!file.exists(config$eu_shp_path)) {
    stop_with_message(paste("Shapefile not found:", config$eu_shp_path))
  }
  if (!file.exists(config$input_csv)) {
    stop_with_message(paste("Input CSV not found:", config$input_csv))
  }
}

read_inputs <- function(config) {
  EU.map <- st_read(config$eu_shp_path, quiet = TRUE)
  EU_centroids <- suppressWarnings(st_centroid(EU.map))
  EU_center_coordinate_df <- data.frame(
    Id = EU_centroids$Id,
    lat = st_coordinates(EU_centroids)[, "Y"],
    lon = st_coordinates(EU_centroids)[, "X"]
  )

  land_cover_data <- read.csv(config$input_csv)
  land_cover_data <- merge(land_cover_data, EU_center_coordinate_df, by = "Id")
  land_cover_data <- land_cover_data[order(
    land_cover_data$Id,
    land_cover_data$year_number,
    land_cover_data$month_number
  ), ]
  list(
    land_cover_data = land_cover_data,
    grid_Id_table = table(land_cover_data$Id)
  )
}

run_sampling_interpolation <- function(seed_num, land_cover_data, grid_Id_table, sampling_out_folder) {
  land_cover_interpolation <- data.frame()

  for (j in seq_along(grid_Id_table)) {
    set.seed(seed_num)
    Id_num <- names(grid_Id_table)[j]
    single_grid_data <- land_cover_data[land_cover_data$Id == Id_num, ]
    grass_MEAN <- single_grid_data$grass_MEAN
    non_na <- !is.na(grass_MEAN)

    rle_result <- rle(non_na)
    segments <- which(rle_result$values & rle_result$lengths >= 3)

    continuous_segments <- lapply(segments, function(idx) {
      start <- if (idx == 1) 1 else sum(rle_result$lengths[1:(idx - 1)]) + 1
      end <- start + rle_result$lengths[idx] - 1
      list(start = start, end = end)
    })

    sample_positions <- lapply(continuous_segments, function(segment) {
      start <- segment$start
      end <- segment$end
      length <- end - start + 1

      if (length == 3) {
        sampled <- start + 1
      } else if (length >= 4 && length <= 6) {
        sampled <- sample((start + 1):(end - 1), 1)
      } else if (length >= 7 && length <= 12) {
        sampled <- sample((start + 1):(end - 1), 2)
      } else if (length > 12) {
        candidates <- (start + 1):(end - 1)
        sampled <- sample(candidates, 2)
        sampled <- sort(sampled)
        if (sampled[2] - sampled[1] != 1) {
          sampled[2] <- sampled[1] + 1
        }
        other_candidates <- setdiff(candidates, sampled)
        sampled <- c(sampled, sample(other_candidates, 1))
      } else {
        sampled <- NA
      }
      sampled
    })

    sample_positions_df <- do.call(rbind, lapply(seq_along(sample_positions), function(i) {
      data.frame(segment_id = i, position = sample_positions[[i]])
    }))

    if (is.null(sample_positions_df) || nrow(sample_positions_df) == 0) {
      next
    }

    land_cover_values <- list(
      grass = single_grid_data$grass_MEAN,
      shrub_and_scrub = single_grid_data$shrub_and_scrub_MEAN,
      snow_and_ice = single_grid_data$snow_and_ice_MEAN,
      trees = single_grid_data$trees_MEAN,
      water = single_grid_data$water_MEAN,
      bare = single_grid_data$bare_MEAN,
      built = single_grid_data$built_MEAN,
      crops = single_grid_data$crops_MEAN,
      flooded_vegetation = single_grid_data$flooded_vegetation_MEAN
    )

    final_interpolated_results <- do.call(rbind, lapply(unique(sample_positions_df$segment_id), function(seg_id) {
      segment_positions <- sample_positions_df[sample_positions_df$segment_id == seg_id, "position"]
      start_idx <- continuous_segments[[seg_id]]$start
      end_idx <- continuous_segments[[seg_id]]$end

      result_rows <- data.frame(segment_id = seg_id, position = segment_positions)
      for (lc_type in all_land_cover_types) {
        segment_values <- land_cover_values[[lc_type]]
        current_segment_values <- segment_values[start_idx:end_idx]
        current_segment_values[segment_positions - start_idx + 1] <- NA
        interpolated_segment <- na.approx(current_segment_values, na.rm = FALSE)
        interpolated_values <- interpolated_segment[segment_positions - start_idx + 1]
        original_values <- segment_values[segment_positions]
        residuals <- abs(interpolated_values - original_values)

        result_rows[[paste0(lc_type, "_original")]] <- original_values
        result_rows[[paste0(lc_type, "_interpolated")]] <- interpolated_values
        result_rows[[paste0(lc_type, "_residual")]] <- residuals
      }
      result_rows
    }))

    single_grid_result <- single_grid_data[, 1:3]
    index <- 3
    for (lc_type in all_land_cover_types) {
      index <- index + 1
      single_grid_result[[paste0(lc_type, "_original")]] <- single_grid_data[, index]
      single_grid_result[[paste0(lc_type, "_interpolated")]] <- NA
      single_grid_result[[paste0(lc_type, "_residual")]] <- NA
    }

    for (i in seq_len(nrow(final_interpolated_results))) {
      index <- final_interpolated_results$position[i]
      for (lc_type in all_land_cover_types) {
        single_grid_result[[paste0(lc_type, "_interpolated")]][index] <- final_interpolated_results[[paste0(lc_type, "_interpolated")]][i]
        single_grid_result[[paste0(lc_type, "_residual")]][index] <- final_interpolated_results[[paste0(lc_type, "_residual")]][i]
      }
    }

    land_cover_interpolation <- rbind(land_cover_interpolation, single_grid_result)
  }

  out_path <- file.path(sampling_out_folder, paste0("land_cover_seed", seed_num, ".csv"))
  write.csv(land_cover_interpolation, out_path, row.names = FALSE)
  paste("Seed", seed_num, "sampling finished")
}

run_imputation <- function(seed_num, sampling_out_folder, land_cover_data, paths, selected_land_cover_types) {
  lc_path <- file.path(sampling_out_folder, paste0("land_cover_seed", seed_num, ".csv"))
  if (!file.exists(lc_path)) {
    stop_with_message(paste("Sampling file not found for seed", seed_num, ":", lc_path))
  }

  lc_imputation <- read.csv(lc_path)
  imputation_cols_to_check <- 4:30
  lc_imputation_train_test_data <- lc_imputation[!apply(
    lc_imputation[, imputation_cols_to_check],
    1,
    function(row) all(is.na(row))
  ), ]

  lc_imputation_predict_data <- lc_imputation[apply(
    lc_imputation[, imputation_cols_to_check],
    1,
    function(row) all(is.na(row))
  ), ]

  lc_imputation_test_data <- lc_imputation_train_test_data[!is.na(lc_imputation_train_test_data$grass_interpolated), ]
  lc_imputation_train_data <- lc_imputation_train_test_data[is.na(lc_imputation_train_test_data$grass_interpolated), ]

  lc_test_data <- lc_imputation_test_data[, c(1:3)]
  lc_train_data <- lc_imputation_train_data[, c(1:3)]
  lc_predict_data <- lc_imputation_predict_data[, c(1:3)]

  lc_train_data <- merge(
    lc_train_data,
    land_cover_data,
    by = c("Id", "year_number", "month_number"),
    all.x = FALSE
  )
  lc_test_data <- merge(
    lc_test_data,
    land_cover_data,
    by = c("Id", "year_number", "month_number"),
    all.x = FALSE
  )
  lc_predict_data <- merge(
    lc_predict_data,
    land_cover_data,
    by = c("Id", "year_number", "month_number"),
    all.x = FALSE
  )

  lc_train_data$Id <- as.factor(lc_train_data$Id)
  lc_train_data$year_number <- as.factor(lc_train_data$year_number)
  lc_train_data$month_number <- as.factor(lc_train_data$month_number)

  lc_test_data$Id <- as.factor(lc_test_data$Id)
  lc_test_data$year_number <- as.factor(lc_test_data$year_number)
  lc_test_data$month_number <- as.factor(lc_test_data$month_number)

  lc_predict_data$Id <- as.factor(lc_predict_data$Id)
  lc_predict_data$year_number <- as.factor(lc_predict_data$year_number)
  lc_predict_data$month_number <- as.factor(lc_predict_data$month_number)

  exclude_columns <- c(
    "Id", "year_number", "month_number", "grass_MEAN",
    "shrub_and_scrub_MEAN", "snow_and_ice_MEAN", "trees_MEAN",
    "water_MEAN", "bare_MEAN", "built_MEAN", "crops_MEAN",
    "flooded_vegetation_MEAN"
  )
  remaining_columns <- setdiff(names(lc_train_data), exclude_columns)
  train_features <- as.matrix(lc_train_data[, remaining_columns])
  test_features <- as.matrix(lc_test_data[, remaining_columns])
  predict_features <- as.matrix(lc_predict_data[, remaining_columns])

  output_two_method_df <- lc_imputation_test_data[, 1:3]
  output_ML_predict_df <- lc_predict_data[, 1:3]
  output_validation_df <- data.frame()

  for (lc_type in selected_land_cover_types) {
    set.seed(123)
    train_target <- log(
      lc_train_data[[paste0(lc_type, "_MEAN")]] /
        (1 - lc_train_data[[paste0(lc_type, "_MEAN")]])
    )
    train_target <- as.numeric(train_target)

    dtrain <- xgb.DMatrix(data = train_features, label = train_target)
    dtest <- xgb.DMatrix(data = test_features)
    dpred <- xgb.DMatrix(data = predict_features)

    xgb_model <- xgboost(
      data = dtrain,
      objective = "reg:squarederror",
      max_depth = 15,
      eta = 0.1,
      nrounds = 300,
      subsample = 1,
      colsample_bytree = 1,
      nthread = 5,
      verbose = 0
    )

    predictions_test <- 1 / (1 + exp(-predict(xgb_model, dtest)))
    na_predictions <- 1 / (1 + exp(-predict(xgb_model, dpred)))

    ML_rmse_test <- sqrt(mean((predictions_test - lc_test_data[[paste0(lc_type, "_MEAN")]])^2))
    ML_mae_test <- mean(abs(predictions_test - lc_test_data[[paste0(lc_type, "_MEAN")]]))
    linear_rmse_test <- sqrt(mean((lc_imputation_test_data[[paste0(lc_type, "_residual")]])^2))
    linear_mae_test <- mean(lc_imputation_test_data[[paste0(lc_type, "_residual")]])

    message(paste("Seed", seed_num, "-", lc_type, "ML_RMSE:", ML_rmse_test))
    message(paste("Seed", seed_num, "-", lc_type, "ML_mae:", ML_mae_test))

    ML_column_name <- paste0(lc_type, "_ML_interpolation")
    linear_column_name <- paste0(lc_type, "_linear_interpolation")
    original_column_name <- paste0(lc_type, "_original")
    predict_ML_column_name <- paste0(lc_type, "_MEAN")

    ML_comparison_df <- cbind(
      lc_test_data[, 1:3],
      setNames(data.frame(predictions_test), ML_column_name)
    )
    linear_comparison_df <- cbind(
      lc_imputation_test_data[, 1:3],
      setNames(data.frame(lc_imputation_test_data[[paste0(lc_type, "_original")]]), original_column_name),
      setNames(data.frame(lc_imputation_test_data[[paste0(lc_type, "_interpolated")]]), linear_column_name)
    )

    two_method_df <- merge(
      linear_comparison_df,
      ML_comparison_df,
      by = c("Id", "year_number", "month_number"),
      all.x = FALSE
    )
    two_method_df[[paste0(lc_type, "_linear_residual")]] <- abs(
      two_method_df[[paste0(lc_type, "_original")]] -
        two_method_df[[paste0(lc_type, "_linear_interpolation")]]
    )
    two_method_df[[paste0(lc_type, "_ML_residual")]] <- abs(
      two_method_df[[paste0(lc_type, "_original")]] -
        two_method_df[[paste0(lc_type, "_ML_interpolation")]]
    )

    ML_predict_df <- cbind(
      lc_predict_data[, 1:3],
      setNames(data.frame(na_predictions), predict_ML_column_name)
    )
    validation_df <- data.frame(
      linear_rmse = linear_rmse_test,
      ML_rmse = ML_rmse_test,
      linear_mae = linear_mae_test,
      ML_mae = ML_mae_test,
      row.names = lc_type
    )

    output_two_method_df <- merge(
      output_two_method_df,
      two_method_df,
      by = c("Id", "year_number", "month_number"),
      all.x = FALSE
    )
    output_validation_df <- rbind(output_validation_df, validation_df)
    output_ML_predict_df <- merge(
      output_ML_predict_df,
      ML_predict_df,
      by = c("Id", "year_number", "month_number"),
      all.x = FALSE
    )

    write.csv(
      ML_predict_df,
      file.path(paths$ml_imputation_output_folder, lc_type, paste0("seed", seed_num, ".csv")),
      row.names = FALSE
    )
  }

  write.csv(
    output_validation_df,
    file.path(paths$two_method_performance_output_folder, paste0("seed", seed_num, ".csv"))
  )
  write.csv(
    output_two_method_df,
    file.path(paths$two_method_test_imputation_output_folder, paste0("seed", seed_num, ".csv")),
    row.names = FALSE
  )
  write.csv(
    output_ML_predict_df,
    file.path(paths$ml_imputation_output_folder, paste0("seed", seed_num, ".csv")),
    row.names = FALSE
  )

  paste("Seed", seed_num, "imputation finished")
}

run_seed_process <- function(seed_num, land_cover_data, grid_Id_table, paths, selected_land_cover_types) {
  run_sampling_interpolation(seed_num, land_cover_data, grid_Id_table, paths$sampling_out_folder)
  run_imputation(seed_num, paths$sampling_out_folder, land_cover_data, paths, selected_land_cover_types)
  paste("Seed", seed_num, "all finished")
}

run_sampling_parallel <- function(seed_list, n_cores, land_cover_data, grid_Id_table, paths) {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c(
    "land_cover_data", "grid_Id_table", "run_sampling_interpolation",
    "all_land_cover_types", "paths"
  ), envir = environment())
  clusterEvalQ(cl, {
    library(zoo)
    NULL
  })
  pbapply::pblapply(seed_list, function(seed_num) {
    run_sampling_interpolation(seed_num, land_cover_data, grid_Id_table, paths$sampling_out_folder)
  }, cl = cl)
}

run_imputation_parallel <- function(seed_list, n_cores, land_cover_data, paths, selected_land_cover_types) {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c(
    "run_imputation", "land_cover_data", "selected_land_cover_types",
    "paths", "stop_with_message"
  ), envir = environment())
  clusterEvalQ(cl, {
    library(xgboost)
    NULL
  })
  pbapply::pblapply(seed_list, function(seed_num) {
    run_imputation(seed_num, paths$sampling_out_folder, land_cover_data, paths, selected_land_cover_types)
  }, cl = cl)
}

aggregate_seed_predictions <- function(seed_list, land_cover_data, paths, selected_land_cover_types) {
  all_lc_df <- land_cover_data[, c("Id", "year_number", "month_number", paste0(all_land_cover_types, "_MEAN"))]

  for (lc_type in selected_land_cover_types) {
    single_lc_df <- land_cover_data[, c("Id", "year_number", "month_number", paste0(lc_type, "_MEAN"))]
    temp_lc_df <- single_lc_df[, c("Id", "year_number", "month_number")]

    for (seed_num in seed_list) {
      seed_path <- file.path(paths$ml_imputation_output_folder, lc_type, paste0("seed", seed_num, ".csv"))
      if (!file.exists(seed_path)) {
        fallback_seed_path <- file.path(paths$ml_imputation_output_folder, paste0("seed", seed_num, ".csv"))
        if (file.exists(fallback_seed_path)) {
          seed_path <- fallback_seed_path
        } else {
          stop_with_message(paste("Missing ML prediction file for aggregation:", seed_path))
        }
      }
      lc_ip_df <- read.csv(seed_path)
      lc_ip_df <- lc_ip_df[, c("Id", "year_number", "month_number", paste0(lc_type, "_MEAN"))]
      colnames(lc_ip_df)[4] <- paste0("ML_", lc_type, "_MEAN_", seed_num)
      temp_lc_df <- merge(
        temp_lc_df,
        lc_ip_df,
        by = c("Id", "year_number", "month_number"),
        all = TRUE
      )
    }

    result_df <- temp_lc_df %>%
      group_by(Id, year_number, month_number) %>%
      summarise(
        mean_value = mean(c_across(starts_with(paste0("ML_", lc_type, "_MEAN_"))), na.rm = TRUE),
        sd_value = sd(c_across(starts_with(paste0("ML_", lc_type, "_MEAN_"))), na.rm = TRUE),
        .groups = "drop"
      )

    single_lc_df <- merge(
      single_lc_df,
      result_df[, 1:4],
      by = c("Id", "year_number", "month_number"),
      all = TRUE
    )
    single_lc_df[, 4][is.na(single_lc_df[, 4]) & !is.na(single_lc_df[, 5])] <-
      single_lc_df[, 5][is.na(single_lc_df[, 4]) & !is.na(single_lc_df[, 5])]
    single_lc_df <- single_lc_df[, 1:4]
    single_lc_df <- single_lc_df[!is.na(single_lc_df[, 4]), ]
    single_lc_df <- single_lc_df[order(
      single_lc_df$Id,
      single_lc_df$year_number,
      single_lc_df$month_number
    ), ]

    all_lc_df <- all_lc_df[, setdiff(names(all_lc_df), paste0(lc_type, "_MEAN")), drop = FALSE]
    all_lc_df <- merge(
      all_lc_df,
      single_lc_df,
      by = c("Id", "year_number", "month_number"),
      all = TRUE
    )
  }

  all_lc_df_copy <- all_lc_df
  selected_output_columns <- match(paste0(selected_land_cover_types, "_MEAN"), names(all_lc_df_copy))
  all_lc_df_copy <- all_lc_df_copy[complete.cases(all_lc_df_copy[, selected_output_columns, drop = FALSE]), ]
  final_output_csv <- get_final_output_path(paths, selected_land_cover_types)
  write.csv(all_lc_df_copy, final_output_csv, row.names = FALSE)
  final_output_csv
}

main <- function() {
  config <- parse_args()
  load_required_packages()
  validate_input_paths(config)
  paths <- build_paths(config$output_folder)
  ensure_output_dirs(paths)

  message("Reading shapefile and land cover data...")
  input_data <- read_inputs(config)
  land_cover_data <- input_data$land_cover_data
  grid_Id_table <- input_data$grid_Id_table

  if (config$mode == "sampling") {
    message("Running sampling interpolation validation...")
    results <- run_sampling_parallel(
      config$seed_list, config$n_cores, land_cover_data, grid_Id_table, paths
    )
    cat(unlist(results), sep = "\n")
  } else if (config$mode == "imputation") {
    message("Running ML imputation...")
    results <- run_imputation_parallel(
      config$seed_list, config$n_cores, land_cover_data, paths, config$land_cover_types
    )
    cat(unlist(results), sep = "\n")
  } else if (config$mode == "aggregate") {
    message("Aggregating multi-seed predictions...")
    final_path <- aggregate_seed_predictions(config$seed_list, land_cover_data, paths, config$land_cover_types)
    message("Final imputed CSV written to: ", final_path)
  } else if (config$mode == "all") {
    message("Running sampling interpolation validation...")
    sampling_results <- run_sampling_parallel(
      config$seed_list, config$n_cores, land_cover_data, grid_Id_table, paths
    )
    cat(unlist(sampling_results), sep = "\n")
    message("Running ML imputation...")
    imputation_results <- run_imputation_parallel(
      config$seed_list, config$n_cores, land_cover_data, paths, config$land_cover_types
    )
    cat(unlist(imputation_results), sep = "\n")
    message("Aggregating multi-seed predictions...")
    final_path <- aggregate_seed_predictions(config$seed_list, land_cover_data, paths, config$land_cover_types)
    message("Final imputed CSV written to: ", final_path)
  }
}

main()
