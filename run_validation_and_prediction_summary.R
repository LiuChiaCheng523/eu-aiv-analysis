get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

trim_ws <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

print_usage <- function(script_name = "run_validation_and_prediction_summary.R") {
  cat(
    paste(
      "Usage:",
      sprintf("  Rscript %s --species \"Anas crecca\"", script_name),
      sprintf("  Rscript %s --all-species", script_name),
      sprintf("  Rscript %s --species-file species.txt --n-iter 5 --mad-k 1.5", script_name),
      sprintf("  Rscript %s --list-species", script_name),
      "",
      "Options:",
      "  --species             Scientific name(s), comma-separated or repeated.",
      "  --species-file        TXT/CSV file containing scientific names.",
      "  --all-species         Process all species available under abundance_prediction.",
      "  --list-species        List species currently available for post-processing and exit.",
      "  --abundance-root      Root folder created by run_abundance_prediction.R.",
      "  --validation-folder   Folder containing validation_data CSV files.",
      "  --checklist-folder    Folder containing checklist CSV files.",
      "  --bird-type-csv       Bird type lookup CSV for future plotting support.",
      "  --eu-shp-path         EU grid shapefile path for abundance maps.",
      "  --output-folder       Output folder. Default: project root.",
      "  --n-iter              Number of iterations to combine. Default: auto-detect.",
      "  --mad-k               MAD filter multiplier. Default: 1.5.",
      "  --map-years           Year(s) to plot for species abundance maps, comma-separated. Default: auto-detect all years.",
      "  --map-months          Month(s) to plot for species abundance maps, comma-separated. Default: 1,4,7,10.",
      "  --help                Show this help message.",
      "",
      "Note:",
      "  This CLI only works for species that already have outputs under abundance_prediction.",
      "  If a species is missing there, run run_abundance_prediction.R first.",
      sep = "\n"
    )
  )
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
    abundance_root = file.path(script_dir, "abundance_spatiotemporal_sampling_method"),
    validation_folder = file.path(script_dir, "validation_data"),
    checklist_folder = file.path(script_dir, "ebird_filtered_checklist"),
    bird_type_csv = file.path(script_dir, "bird_type_lookup.csv"),
    eu_shp_path = file.path(script_dir, "EU_100km_fishnet_simple_by_distance", "EU_100km_fishnet_simple_by_distance.shp"),
    output_folder = script_dir,
    n_iter = NA_integer_,
    mad_k = 1.5,
    map_years = character(),
    map_months = c(1L, 4L, 7L, 10L)
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
      split_species <- trim_ws(unlist(strsplit(value, ",", fixed = TRUE)))
      opts$species <- c(opts$species, split_species[nzchar(split_species)])
    } else if (key == "species-file") {
      opts$species_file <- value
    } else if (key == "abundance-root") {
      opts$abundance_root <- value
    } else if (key == "validation-folder") {
      opts$validation_folder <- value
    } else if (key == "checklist-folder") {
      opts$checklist_folder <- value
    } else if (key == "bird-type-csv") {
      opts$bird_type_csv <- value
    } else if (key == "eu-shp-path") {
      opts$eu_shp_path <- value
    } else if (key == "output-folder") {
      opts$output_folder <- value
    } else if (key == "n-iter") {
      opts$n_iter <- as.integer(value)
    } else if (key == "mad-k") {
      opts$mad_k <- as.numeric(value)
    } else if (key == "map-years") {
      opts$map_years <- trim_ws(unlist(strsplit(value, ",", fixed = TRUE)))
    } else if (key == "map-months") {
      opts$map_months <- as.integer(trim_ws(unlist(strsplit(value, ",", fixed = TRUE))))
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }

    i <- i + 2
  }

  opts$abundance_root <- normalizePath(opts$abundance_root, winslash = "/", mustWork = FALSE)
  opts$validation_folder <- normalizePath(opts$validation_folder, winslash = "/", mustWork = FALSE)
  opts$checklist_folder <- normalizePath(opts$checklist_folder, winslash = "/", mustWork = FALSE)
  opts$bird_type_csv <- normalizePath(opts$bird_type_csv, winslash = "/", mustWork = FALSE)
  opts$eu_shp_path <- normalizePath(opts$eu_shp_path, winslash = "/", mustWork = FALSE)
  opts$output_folder <- normalizePath(opts$output_folder, winslash = "/", mustWork = FALSE)
  opts$species <- unique(trim_ws(opts$species))
  opts$species <- opts$species[nzchar(opts$species)]
  opts$map_years <- unique(trim_ws(opts$map_years))
  opts$map_years <- opts$map_years[nzchar(opts$map_years)]
  opts$map_months <- unique(opts$map_months[!is.na(opts$map_months)])
  opts
}

extract_species_name_from_file <- function(filename) {
  name <- sub("\\.csv$", "", filename, ignore.case = TRUE)
  name <- sub("^\\^", "", name)
  name <- sub("_filtered_2019to2022$", "", name)
  name <- sub("_filtered2019to2022$", "", name)
  trim_ws(name)
}

build_checklist_catalog <- function(checklist_folder) {
  csv_paths <- list.files(checklist_folder, pattern = "\\.csv$", full.names = TRUE)
  data.frame(
    scientific_name = vapply(basename(csv_paths), extract_species_name_from_file, character(1)),
    checklist_path = normalizePath(csv_paths, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

build_available_species_catalog <- function(abundance_root) {
  abundance_prediction_root <- file.path(abundance_root, "abundance_prediction")
  if (!dir.exists(abundance_prediction_root)) {
    stop(
      "abundance_prediction folder not found under: ", abundance_root,
      "\nPlease run run_abundance_prediction.R first.",
      call. = FALSE
    )
  }

  dirs <- list.dirs(abundance_prediction_root, recursive = FALSE, full.names = TRUE)
  data.frame(
    scientific_name = basename(dirs),
    abundance_prediction_path = normalizePath(dirs, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

read_species_file <- function(path) {
  if (!file.exists(path)) {
    stop("Species file not found: ", path, call. = FALSE)
  }
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if ("scientific_name" %in% names(df)) {
      species <- df$scientific_name
    } else if ("birdname" %in% names(df)) {
      species <- df$birdname
    } else {
      species <- df[[1]]
    }
  } else {
    species <- readLines(path, warn = FALSE)
  }
  species <- trim_ws(species)
  unique(species[nzchar(species)])
}

resolve_species <- function(opts, available_catalog) {
  if (opts$all_species) {
    return(sort(available_catalog$scientific_name))
  }

  requested <- opts$species
  if (!is.null(opts$species_file)) {
    requested <- c(requested, read_species_file(opts$species_file))
  }
  requested <- unique(trim_ws(requested))
  requested <- requested[nzchar(requested)]

  if (length(requested) == 0) {
    stop("Please provide --species, --species-file, --all-species, or --list-species.", call. = FALSE)
  }

  missing <- setdiff(requested, available_catalog$scientific_name)
  if (length(missing) > 0) {
    stop(
      "These species do not have abundance_prediction outputs yet: ",
      paste(missing, collapse = ", "),
      "\nPlease run run_abundance_prediction.R first.",
      call. = FALSE
    )
  }

  sort(requested)
}

count_iteration_files <- function(folder_path, pattern_prefix) {
  files <- list.files(folder_path, pattern = paste0("^", pattern_prefix, "\\d+\\.csv$"), full.names = FALSE)
  if (length(files) == 0) {
    return(0L)
  }
  as.integer(length(files))
}

prepare_bird_type_lookup <- function(path, checklist_catalog) {
  if (!file.exists(path)) {
    warning("bird type CSV not found: ", path)
    return(NULL)
  }
  bird_type_df <- read.csv(path, stringsAsFactors = FALSE)
  bird_type_df <- bird_type_df[bird_type_df$birdname %in% checklist_catalog$scientific_name, , drop = FALSE]
  bird_type_df
}

valid_spearman_for_ensemble_mad <- function(folder_path, validation_data_df, n_iter, k = 1.5) {
  merged_df <- read.csv(file.path(folder_path, "abundance_iteration1.csv"), stringsAsFactors = FALSE)
  merged_df <- merged_df[, c(4, 5)]

  if (n_iter == 1) {
    merged_df <- merge(merged_df, validation_data_df, by = "sample_key", all = FALSE)
    return(cor(merged_df$abundance1, merged_df$observation_count, method = "spearman"))
  }

  for (i in 2:n_iter) {
    temp_df <- read.csv(file.path(folder_path, sprintf("abundance_iteration%d.csv", i)), stringsAsFactors = FALSE)
    col_name <- colnames(temp_df)[5]
    merged_df <- merge(
      merged_df,
      temp_df[, c(4, 5)],
      by = "sample_key",
      all = FALSE
    )
    colnames(merged_df)[ncol(merged_df)] <- col_name
  }

  abundance_cols <- grep("^abundance\\d+$", colnames(merged_df), value = TRUE)
  merged_df <- merge(
    merged_df,
    validation_data_df[, c("sample_key", "observation_count")],
    by = "sample_key",
    all = FALSE
  )

  abundance_filtered <- merged_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      abundance_median = median(dplyr::c_across(dplyr::all_of(abundance_cols)), na.rm = TRUE),
      abundance_mad = median(abs(dplyr::c_across(dplyr::all_of(abundance_cols)) - abundance_median), na.rm = TRUE),
      abundance_filtered_mean = mean(
        dplyr::c_across(dplyr::all_of(abundance_cols))[
          abs(dplyr::c_across(dplyr::all_of(abundance_cols)) - abundance_median) <= (k * abundance_mad)
        ],
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup()

  cor(abundance_filtered$abundance_filtered_mean, abundance_filtered$observation_count, method = "spearman")
}

combine_abundance_predictions_MAD <- function(folder_path, n_iter, k = 1.5) {
  merged_df <- read.csv(file.path(folder_path, "abundance_iteration1.csv"), stringsAsFactors = FALSE)
  colnames(merged_df)[4] <- "abundance1"

  for (i in 2:n_iter) {
    temp_df <- read.csv(file.path(folder_path, sprintf("abundance_iteration%d.csv", i)), stringsAsFactors = FALSE)
    colnames(temp_df)[4] <- paste0("abundance", i)
    merged_df <- merge(merged_df, temp_df, by = c("Id", "year_number", "month_number"), all = TRUE)
  }

  abundance_cols <- grep("^abundance\\d+$", colnames(merged_df), value = TRUE)

  merged_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      abundance_median = median(dplyr::c_across(dplyr::all_of(abundance_cols)), na.rm = TRUE),
      abundance_mad = median(abs(dplyr::c_across(dplyr::all_of(abundance_cols)) - abundance_median), na.rm = TRUE),
      abundance_filtered_mean = mean(
        dplyr::c_across(dplyr::all_of(abundance_cols))[
          abs(dplyr::c_across(dplyr::all_of(abundance_cols)) - abundance_median) <= (k * abundance_mad)
        ],
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup()
}

validate_options <- function(opts) {
  if (!dir.exists(opts$abundance_root)) {
    stop("abundance root not found: ", opts$abundance_root, call. = FALSE)
  }
  if (!dir.exists(opts$validation_folder)) {
    stop("validation folder not found: ", opts$validation_folder, call. = FALSE)
  }
  if (!dir.exists(opts$checklist_folder)) {
    stop("checklist folder not found: ", opts$checklist_folder, call. = FALSE)
  }
  if (!file.exists(opts$eu_shp_path)) {
    stop("EU shapefile not found: ", opts$eu_shp_path, call. = FALSE)
  }
  if (!is.na(opts$n_iter) && opts$n_iter < 1) {
    stop("--n-iter must be at least 1.", call. = FALSE)
  }
  if (is.na(opts$mad_k) || opts$mad_k <= 0) {
    stop("--mad-k must be greater than 0.", call. = FALSE)
  }
  if (length(opts$map_months) == 0) {
    stop("--map-months must include at least one month.", call. = FALSE)
  }
  if (any(opts$map_months < 1 | opts$map_months > 12)) {
    stop("--map-months must be between 1 and 12.", call. = FALSE)
  }
}

load_eu_map <- function(shp_path) {
  sf::st_read(shp_path, quiet = TRUE)
}

plot_validation_summary <- function(summary_df, output_path) {
  plot_df <- summary_df
  plot_df$label <- plot_df$birdname

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = validation_n, y = spatial_SRC, color = bird_type)) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = label), nudge_y = 0.01, size = 3.5, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30", linewidth = 0.5) +
    ggplot2::scale_color_manual(
      values = c("waterbird" = "#1f78b4", "non_waterbird" = "#33a02c"),
      na.value = "grey50"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      x = "validation sample size",
      y = "Spearman's rank correlation",
      color = "bird type",
      title = "Validation performance summary"
    )

  ggplot2::ggsave(output_path, p, width = 10, height = 7, dpi = 300)
}

plot_species_map <- function(eu_map, abundance_df, bird_name, map_year, map_month, output_path) {
  year_value <- suppressWarnings(as.integer(as.character(abundance_df$year_number)))
  month_value <- suppressWarnings(as.integer(as.character(abundance_df$month_number)))
  map_slice <- abundance_df[year_value == map_year & month_value == map_month, , drop = FALSE]

  if (nrow(map_slice) == 0) {
    warning("No abundance rows found for map slice: ", bird_name, " ", map_year, "-", map_month)
    return(invisible(NULL))
  }

  eu_map_plot <- dplyr::left_join(eu_map, map_slice, by = "Id")
  p <- ggplot2::ggplot(eu_map_plot) +
    ggplot2::geom_sf(ggplot2::aes(fill = abundance_filtered_mean), color = NA) +
    viridis::scale_fill_viridis(option = "C", direction = -1, na.value = "lightgray") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0(bird_name, " relative abundance ", map_year, "-", sprintf("%02d", map_month)),
      fill = "abundance"
    )

  ggplot2::ggsave(output_path, p, width = 10, height = 7, dpi = 300)
}

build_feature_importance_summary <- function(feature_importance_root, species_list, bird_type_df, n_iter_limit = NULL) {
  vegetation_feature_vector <- c(
    "water_MEAN", "trees_MEAN", "crops_MEAN", "bare_MEAN", "built_MEAN",
    "grass_MEAN", "snow_and_ice_MEAN", "flooded_vegetation_MEAN", "shrub_and_scrub_MEAN"
  )

  if (!dir.exists(feature_importance_root)) {
    warning("Feature importance folder not found: ", feature_importance_root)
    return(NULL)
  }

  classified_species <- if (!is.null(bird_type_df)) {
    intersect(species_list, bird_type_df$birdname)
  } else {
    character()
  }

  if (length(classified_species) == 0) {
    warning("No classified species available for FI radar plot.")
    return(NULL)
  }

  get_group_sum <- function(group_species) {
    sum_df <- data.frame(Feature = vegetation_feature_vector, sum_rank_point = 0, stringsAsFactors = FALSE)

    for (birdname in group_species) {
      bird_folder <- file.path(feature_importance_root, birdname)
      if (!dir.exists(bird_folder)) next

      iter_files <- list.files(bird_folder, pattern = "^iteration\\d+\\.csv$", full.names = TRUE)
      if (length(iter_files) == 0) next
      iter_files <- iter_files[order(iter_files)]
      if (!is.null(n_iter_limit) && !is.na(n_iter_limit)) {
        iter_files <- head(iter_files, n_iter_limit)
      }

      for (iter_file in iter_files) {
        fi_df <- tryCatch(read.csv(iter_file, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(fi_df) || !all(c("Feature", "Gain") %in% names(fi_df))) next

        fi_df <- fi_df[fi_df$Feature %in% vegetation_feature_vector, c("Feature", "Gain"), drop = FALSE]
        if (nrow(fi_df) == 0) next

        fi_df$rank_point <- rank(-fi_df$Gain, ties.method = "first")
        fi_df$rank_point <- 10 - fi_df$rank_point
        fi_df <- fi_df[, c("Feature", "rank_point"), drop = FALSE]

        sum_df <- merge(sum_df, fi_df, by = "Feature", all = TRUE)
        sum_df$sum_rank_point <- sum_df$sum_rank_point + ifelse(is.na(sum_df$rank_point), 0, sum_df$rank_point)
        sum_df <- sum_df[, c("Feature", "sum_rank_point")]
      }
    }

    sum_df[order(-sum_df$sum_rank_point), , drop = FALSE]
  }

  water_species <- intersect(
    classified_species,
    bird_type_df$birdname[bird_type_df$bird_type == "waterbird"]
  )
  nonwater_species <- intersect(
    classified_species,
    bird_type_df$birdname[bird_type_df$bird_type == "non_waterbird"]
  )

  all_sum_fi <- get_group_sum(classified_species)
  water_sum_fi <- get_group_sum(water_species)
  nonwater_sum_fi <- get_group_sum(nonwater_species)

  colnames(all_sum_fi)[2] <- "all"
  colnames(water_sum_fi)[2] <- "water"
  colnames(nonwater_sum_fi)[2] <- "non_water"

  sum_fi <- merge(all_sum_fi, water_sum_fi, by = "Feature", all = TRUE)
  sum_fi <- merge(sum_fi, nonwater_sum_fi, by = "Feature", all = TRUE)
  sum_fi[is.na(sum_fi)] <- 0
  sum_fi <- sum_fi[order(-sum_fi$all), , drop = FALSE]
  sum_fi$check <- sum_fi$water + sum_fi$non_water
  sum_fi
}

plot_feature_importance_radar <- function(sum_fi, output_path) {
  if (is.null(sum_fi) || nrow(sum_fi) == 0) {
    return(invisible(NULL))
  }

  radar_raw <- sum_fi[, c("Feature", "water", "non_water", "all")]
  radar_mat <- t(as.matrix(radar_raw[, c("water", "non_water", "all")]))
  colnames(radar_mat) <- radar_raw$Feature
  rownames(radar_mat) <- c("water", "non_water", "all")

  global_max <- max(radar_raw$all, na.rm = TRUE)
  if (!is.finite(global_max) || global_max <= 0) {
    warning("Feature importance summary has no positive values for radar plot.")
    return(invisible(NULL))
  }

  feature_count <- ncol(radar_mat)
  angles <- seq(0, 2 * pi, length.out = feature_count + 1)
  scale_breaks <- round(seq(0, global_max, length.out = 5))

  png(output_path, width = 2000, height = 1600, res = 200)
  par(mar = c(3, 3, 5, 6))
  plot.new()
  plot.window(xlim = c(-1.4, 1.8), ylim = c(-1.3, 1.3), asp = 1)

  for (grid_val in seq(0.2, 1, by = 0.2)) {
    lines(grid_val * cos(angles), grid_val * sin(angles), col = "grey80", lty = 1)
  }
  for (i in seq_len(feature_count)) {
    lines(c(0, cos(angles[i])), c(0, sin(angles[i])), col = "grey70", lty = 1)
    text(1.15 * cos(angles[i]), 1.15 * sin(angles[i]), labels = colnames(radar_mat)[i], cex = 0.8)
  }
  for (i in seq_along(scale_breaks)) {
    y_pos <- (i - 1) / (length(scale_breaks) - 1)
    text(0.03, y_pos, labels = scale_breaks[i], cex = 0.75, col = "grey30")
  }

  draw_series <- function(values, border_col, fill_col) {
    scaled_vals <- if (global_max > 0) values / global_max else values
    coords_x <- c(scaled_vals * cos(angles[-length(angles)]), scaled_vals[1] * cos(angles[1]))
    coords_y <- c(scaled_vals * sin(angles[-length(angles)]), scaled_vals[1] * sin(angles[1]))
    polygon(coords_x, coords_y, border = border_col, col = fill_col, lwd = 2)
    points(scaled_vals * cos(angles[-length(angles)]), scaled_vals * sin(angles[-length(angles)]), pch = 16, col = border_col, cex = 0.8)
  }

  draw_series(radar_mat["water", ], "#1f78b4", grDevices::adjustcolor("#1f78b4", alpha.f = 0.15))
  draw_series(radar_mat["non_water", ], "#33a02c", grDevices::adjustcolor("#33a02c", alpha.f = 0.15))
  draw_series(radar_mat["all", ], "#e31a1c", grDevices::adjustcolor("#e31a1c", alpha.f = 0.15))

  title("radar chart : land cover feature importance ranking")
  legend(
    "topright",
    inset = c(-0.18, 0),
    legend = c("water", "non_water", "all"),
    col = c("#1f78b4", "#33a02c", "#e31a1c"),
    lwd = 3,
    bty = "n",
    xpd = TRUE
  )
  dev.off()
}

get_map_years_to_plot <- function(abundance_df, requested_years) {
  available_years <- sort(unique(suppressWarnings(as.integer(as.character(abundance_df$year_number)))))
  available_years <- available_years[!is.na(available_years)]

  if (length(requested_years) == 0) {
    return(available_years)
  }

  requested_years_int <- suppressWarnings(as.integer(requested_years))
  requested_years_int <- requested_years_int[!is.na(requested_years_int)]
  intersect(requested_years_int, available_years)
}

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_name <- if (length(script_file_arg) > 0) basename(sub("^--file=", "", script_file_arg[1])) else "run_validation_and_prediction_summary.R"
opts <- parse_args()

if (opts$help) {
  print_usage(script_name)
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(viridis)
})

validate_options(opts)
checklist_catalog <- build_checklist_catalog(opts$checklist_folder)
available_catalog <- build_available_species_catalog(opts$abundance_root)
available_catalog <- available_catalog[available_catalog$scientific_name %in% checklist_catalog$scientific_name, , drop = FALSE]
available_catalog <- available_catalog[order(available_catalog$scientific_name), , drop = FALSE]

if (opts$list_species) {
  write.table(available_catalog["scientific_name"], row.names = FALSE, col.names = TRUE, quote = FALSE)
  quit(save = "no", status = 0)
}

selected_species <- resolve_species(opts, available_catalog)
bird_type_df <- prepare_bird_type_lookup(opts$bird_type_csv, checklist_catalog)
eu_map <- load_eu_map(opts$eu_shp_path)

summary_output_dir <- file.path(opts$output_folder, "validation_prediction_summary")
mad_output_dir <- file.path(summary_output_dir, "mad_filter_abundance")
plot_output_dir <- file.path(summary_output_dir, "plots")
species_plot_output_dir <- file.path(plot_output_dir, "species_maps")
feature_importance_root <- file.path(opts$abundance_root, "feature_importance")
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mad_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(species_plot_output_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()

for (bird_name in selected_species) {
  cat("Processing species:", bird_name, "\n")

  validation_prediction_folder <- file.path(opts$abundance_root, "validation_prediction", bird_name)
  abundance_prediction_folder <- file.path(opts$abundance_root, "abundance_prediction", bird_name)
  validation_data_path <- file.path(opts$validation_folder, paste0(bird_name, ".csv"))

  if (!dir.exists(validation_prediction_folder) || !dir.exists(abundance_prediction_folder) || !file.exists(validation_data_path)) {
    stop(
      "Missing required post-processing inputs for ", bird_name,
      ". Please run run_abundance_prediction.R first.",
      call. = FALSE
    )
  }

  available_validation_iter <- count_iteration_files(validation_prediction_folder, "abundance_iteration")
  available_abundance_iter <- count_iteration_files(abundance_prediction_folder, "abundance_iteration")
  available_iter <- min(available_validation_iter, available_abundance_iter)

  if (available_iter < 1) {
    stop(
      "No iteration CSV files found for ", bird_name,
      ". Please run run_abundance_prediction.R first.",
      call. = FALSE
    )
  }

  use_iter <- if (is.na(opts$n_iter)) available_iter else opts$n_iter
  if (use_iter > available_iter) {
    stop(
      "Requested --n-iter=", use_iter, " but only ", available_iter,
      " iteration file(s) are available for ", bird_name,
      ". Please rerun run_abundance_prediction.R with more iterations, or lower --n-iter.",
      call. = FALSE
    )
  }

  validation_data_df <- read.csv(validation_data_path, stringsAsFactors = FALSE)
  spatial_src <- valid_spearman_for_ensemble_mad(
    folder_path = validation_prediction_folder,
    validation_data_df = validation_data_df,
    n_iter = use_iter,
    k = opts$mad_k
  )

  abund_all <- combine_abundance_predictions_MAD(
    folder_path = abundance_prediction_folder,
    n_iter = use_iter,
    k = opts$mad_k
  )
  species_plot_dir <- file.path(species_plot_output_dir, bird_name)
  dir.create(species_plot_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    abund_all,
    file.path(mad_output_dir, paste0(bird_name, ".csv")),
    row.names = FALSE
  )
  map_years_to_plot <- get_map_years_to_plot(abund_all, opts$map_years)
  if (length(map_years_to_plot) == 0) {
    warning("No matching map years available for ", bird_name, ". Skipping species maps.")
  } else {
    safe_bird_name <- gsub("[\\\\/:*?\"<>|]", "_", bird_name)
    for (plot_year in map_years_to_plot) {
      for (plot_month in opts$map_months) {
        plot_species_map(
          eu_map = eu_map,
          abundance_df = abund_all,
          bird_name = bird_name,
          map_year = plot_year,
          map_month = plot_month,
          output_path = file.path(
            species_plot_dir,
            paste0(safe_bird_name, "_", plot_year, "_", sprintf("%02d", plot_month), ".png")
          )
        )
      }
    }
  }

  bird_type_row <- if (!is.null(bird_type_df)) bird_type_df[bird_type_df$birdname == bird_name, , drop = FALSE] else NULL
  bird_type <- if (!is.null(bird_type_row) && nrow(bird_type_row) > 0) bird_type_row$bird_type[1] else NA
  waterbird <- if (!is.null(bird_type_row) && nrow(bird_type_row) > 0) bird_type_row$waterbird[1] else NA

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    birdname = bird_name,
    n_iter = use_iter,
    mad_k = opts$mad_k,
    validation_n = nrow(validation_data_df),
    spatial_SRC = spatial_src,
    bird_type = bird_type,
    waterbird = waterbird,
    stringsAsFactors = FALSE
  )
}

summary_df <- dplyr::bind_rows(summary_rows) %>%
  dplyr::arrange(dplyr::desc(spatial_SRC))

write.csv(summary_df, file.path(summary_output_dir, "validation_summary.csv"), row.names = FALSE)

if (!is.null(bird_type_df)) {
  write.csv(bird_type_df, file.path(summary_output_dir, "bird_type_lookup_intersected.csv"), row.names = FALSE)
}
plot_validation_summary(summary_df, file.path(plot_output_dir, "validation_summary_scatter.png"))

fi_summary_df <- build_feature_importance_summary(
  feature_importance_root = feature_importance_root,
  species_list = selected_species,
  bird_type_df = bird_type_df,
  n_iter_limit = if (is.na(opts$n_iter)) NULL else opts$n_iter
)
if (!is.null(fi_summary_df)) {
  write.csv(fi_summary_df, file.path(summary_output_dir, "feature_importance_land_cover_summary.csv"), row.names = FALSE)
  plot_feature_importance_radar(
    fi_summary_df,
    file.path(plot_output_dir, "feature_importance_land_cover_radar.png")
  )
}

cat("Done.\n")
cat("Summary CSV:", file.path(summary_output_dir, "validation_summary.csv"), "\n")
cat("MAD-filter abundance folder:", mad_output_dir, "\n")
cat("Plot folder:", plot_output_dir, "\n")
