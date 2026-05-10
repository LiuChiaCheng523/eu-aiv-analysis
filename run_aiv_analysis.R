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

print_usage <- function(script_name = "run_aiv_analysis.R") {
  cat(
    paste(
      "Usage:",
      sprintf("  Rscript %s --outbreak-type Wild --species \"Vanellus vanellus\"", script_name),
      sprintf("  Rscript %s --outbreak-type Domestic --species \"Vanellus vanellus,Branta canadensis\"", script_name),
      sprintf("  Rscript %s --outbreak-type Wild --all-species", script_name),
      sprintf("  Rscript %s --list-species", script_name),
      "",
      "Options:",
      "  --outbreak-type         Wild or Domestic.",
      "  --species               Scientific name(s), comma-separated or repeated.",
      "  --species-file          TXT/CSV file containing scientific names.",
      "  --all-species           Process all species available in mad_filter_abundance.",
      "  --list-species          List species currently available for AIV analysis and exit.",
      "  --eu-shp-path           EU grid shapefile path.",
      "  --chicken-density-path  Chicken livestock density CSV path.",
      "  --duck-density-path     Duck livestock density CSV path.",
      "  --aiv-2021-path         AIV fixed data 2021 CSV path.",
      "  --aiv-2022-path         AIV fixed data 2022 CSV path.",
      "  --bird-abundance-folder Folder containing MAD-filter abundance CSVs.",
      "  --output-folder         Output folder. Default: project root.",
      "  --write-date            Date tag for output CSV names. Default: today's date (YYYYMMDD).",
      "  --help                  Show this help message.",
      "",
      "Note:",
      "  This CLI expects species CSVs from run_validation_and_prediction_summary.R under mad_filter_abundance.",
      sep = "\n"
    )
  )
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  script_dir <- get_script_dir()
  opts <- list(
    outbreak_type = NULL,
    species = character(),
    species_file = NULL,
    all_species = FALSE,
    list_species = FALSE,
    help = FALSE,
    eu_shp_path = file.path(script_dir, "EU_100km_fishnet_simple_by_distance", "EU_100km_fishnet_simple_by_distance.shp"),
    chicken_density_path = file.path(script_dir, "livestock_density_10km", "chicken livestock density 10km.csv"),
    duck_density_path = file.path(script_dir, "livestock_density_10km", "duck livestock density_2015_10km.csv"),
    aiv_2021_path = file.path(script_dir, "aiv_fixed_data", "EU aiv fixed data 2021.csv"),
    aiv_2022_path = file.path(script_dir, "aiv_fixed_data", "EU aiv fixed data 2022.csv"),
    bird_abundance_folder = file.path(script_dir, "validation_prediction_summary", "mad_filter_abundance"),
    output_folder = script_dir,
    write_date = format(Sys.Date(), "%Y%m%d")
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

    if (key == "outbreak-type") {
      opts$outbreak_type <- value
    } else if (key == "species") {
      split_species <- trim_ws(unlist(strsplit(value, ",", fixed = TRUE)))
      opts$species <- c(opts$species, split_species[nzchar(split_species)])
    } else if (key == "species-file") {
      opts$species_file <- value
    } else if (key == "eu-shp-path") {
      opts$eu_shp_path <- value
    } else if (key == "chicken-density-path") {
      opts$chicken_density_path <- value
    } else if (key == "duck-density-path") {
      opts$duck_density_path <- value
    } else if (key == "aiv-2021-path") {
      opts$aiv_2021_path <- value
    } else if (key == "aiv-2022-path") {
      opts$aiv_2022_path <- value
    } else if (key == "bird-abundance-folder") {
      opts$bird_abundance_folder <- value
    } else if (key == "output-folder") {
      opts$output_folder <- value
    } else if (key == "write-date") {
      opts$write_date <- value
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }

    i <- i + 2
  }

  opts$species <- unique(trim_ws(opts$species))
  opts$species <- opts$species[nzchar(opts$species)]
  opts$eu_shp_path <- normalizePath(opts$eu_shp_path, winslash = "/", mustWork = FALSE)
  opts$chicken_density_path <- normalizePath(opts$chicken_density_path, winslash = "/", mustWork = FALSE)
  opts$duck_density_path <- normalizePath(opts$duck_density_path, winslash = "/", mustWork = FALSE)
  opts$aiv_2021_path <- normalizePath(opts$aiv_2021_path, winslash = "/", mustWork = FALSE)
  opts$aiv_2022_path <- normalizePath(opts$aiv_2022_path, winslash = "/", mustWork = FALSE)
  opts$bird_abundance_folder <- normalizePath(opts$bird_abundance_folder, winslash = "/", mustWork = FALSE)
  opts$output_folder <- normalizePath(opts$output_folder, winslash = "/", mustWork = FALSE)
  opts
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

build_species_catalog <- function(folder) {
  if (!dir.exists(folder)) {
    stop(
      "Bird abundance folder not found: ", folder,
      "\nPlease run run_validation_and_prediction_summary.R first.",
      call. = FALSE
    )
  }
  csv_paths <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  data.frame(
    scientific_name = sub("\\.csv$", "", basename(csv_paths)),
    abundance_path = normalizePath(csv_paths, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

resolve_species <- function(opts, species_catalog) {
  if (opts$all_species) {
    return(sort(species_catalog$scientific_name))
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

  missing <- setdiff(requested, species_catalog$scientific_name)
  if (length(missing) > 0) {
    stop(
      "These species are not available in mad_filter_abundance: ",
      paste(missing, collapse = ", "),
      "\nPlease run run_validation_and_prediction_summary.R first.",
      call. = FALSE
    )
  }
  sort(requested)
}

validate_options <- function(opts) {
  if (!is.null(opts$outbreak_type)) {
    opts$outbreak_type <- trim_ws(opts$outbreak_type)
  }
  if (!opts$list_species && is.null(opts$outbreak_type)) {
    stop("--outbreak-type is required unless you use --list-species.", call. = FALSE)
  }
  if (!opts$list_species && !opts$outbreak_type %in% c("Wild", "Domestic")) {
    stop("--outbreak-type must be either Wild or Domestic.", call. = FALSE)
  }
  required_files <- c(opts$eu_shp_path, opts$chicken_density_path, opts$duck_density_path, opts$aiv_2021_path, opts$aiv_2022_path)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"), call. = FALSE)
  }
}

get_neighbors_df <- function(eu_map) {
  neighbors <- sf::st_touches(eu_map)
  matrix_data <- matrix(nrow = length(eu_map$Id), ncol = 9)
  colnames(matrix_data) <- c("Id", paste0("neighbor", 1:8))
  matrix_data[, 1] <- eu_map$Id

  for (i in seq_along(eu_map$Id)) {
    neighbors_indices <- neighbors[[i]]
    if (length(neighbors_indices) > 0) {
      neighbors_id <- eu_map[neighbors_indices, ]
      for (j in seq_len(min(length(neighbors_id$Id), 8))) {
        matrix_data[i, j + 1] <- neighbors_id$Id[j]
      }
    }
  }
  as.data.frame(matrix_data)
}

get_eu_center_coordinate_df <- function(eu_map) {
  eu_centroids <- sf::st_centroid(eu_map)
  data.frame(
    Id = eu_centroids$Id,
    latitude = sf::st_coordinates(eu_centroids)[, "Y"],
    longitude = sf::st_coordinates(eu_centroids)[, "X"]
  )
}

get_livestock_sum_density <- function(path) {
  read.csv(path, stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(Id)) %>%
    dplyr::group_by(Id) %>%
    dplyr::summarise(sum_VALUE = sum(VALUE, na.rm = TRUE), .groups = "drop")
}

get_aiv_records <- function(path, outbreak_type, eu_crs) {
  aiv_df <- read.csv(path, stringsAsFactors = FALSE)
  aiv_df <- aiv_df[!is.na(aiv_df$Id), , drop = FALSE]
  aiv_df$Species <- sub(",.*", "", aiv_df$Species)
  aiv_df <- aiv_df[aiv_df$Species == outbreak_type, , drop = FALSE]
  list(
    df = aiv_df,
    sf = sf::st_as_sf(aiv_df, coords = c("longitude", "latitude"), crs = eu_crs)
  )
}

get_aiv_raw_records <- function(path) {
  aiv_df <- read.csv(path, stringsAsFactors = FALSE)
  aiv_df <- aiv_df[!is.na(aiv_df$Id), , drop = FALSE]
  aiv_df$Species <- sub(",.*", "", aiv_df$Species)
  aiv_df
}

build_outbreak_frequency <- function(eu_map_df, eu_aiv_df, neighbors_df, outbreak_type) {
  id_outbreak_freq <- as.data.frame(table(eu_aiv_df$Id))
  colnames(id_outbreak_freq) <- c("Id", "freq21to22")
  id_outbreak_freq$stage1_freq <- 0
  id_outbreak_freq$stage2_freq <- 0
  id_outbreak_freq$stage3_freq <- 0

  id_to_row <- match(eu_aiv_df$Id, id_outbreak_freq$Id)

  for (i in seq_len(nrow(eu_aiv_df))) {
    row_index <- id_to_row[i]
    if (is.na(row_index)) next
    obs_year <- eu_aiv_df$observation_year[i]
    obs_month <- eu_aiv_df$observation_month[i]

    if (obs_year == 2021) {
      if (obs_month >= 5 && obs_month <= 10) {
        id_outbreak_freq$stage1_freq[row_index] <- id_outbreak_freq$stage1_freq[row_index] + 1
      } else if (obs_month >= 11) {
        id_outbreak_freq$stage2_freq[row_index] <- id_outbreak_freq$stage2_freq[row_index] + 1
      }
    } else if (obs_year == 2022) {
      if (obs_month <= 4) {
        id_outbreak_freq$stage2_freq[row_index] <- id_outbreak_freq$stage2_freq[row_index] + 1
      } else if (obs_month >= 5 && obs_month <= 10) {
        id_outbreak_freq$stage3_freq[row_index] <- id_outbreak_freq$stage3_freq[row_index] + 1
      }
    }
  }

  id_df <- data.frame(Id = eu_map_df$Id)
  id_outbreak_freq <- merge(id_df, id_outbreak_freq, by = "Id", all = TRUE)
  id_outbreak_freq$Id <- as.numeric(id_outbreak_freq$Id)
  id_outbreak_freq$stage1_freq[is.na(id_outbreak_freq$stage1_freq)] <- 0
  id_outbreak_freq$stage2_freq[is.na(id_outbreak_freq$stage2_freq)] <- 0
  id_outbreak_freq$stage3_freq[is.na(id_outbreak_freq$stage3_freq)] <- 0
  id_outbreak_freq$stage1_outbreak_nei <- 0
  id_outbreak_freq$stage2_outbreak_nei <- 0
  id_outbreak_freq$stage3_outbreak_nei <- 0

  for (i in seq_len(nrow(id_outbreak_freq))) {
    id_number <- as.numeric(id_outbreak_freq$Id[i])
    id_neighbors <- neighbors_df[neighbors_df$Id == id_number, 2:9]
    id_neighbors <- as.numeric(id_neighbors[!is.na(id_neighbors)])

    if (length(id_neighbors) > 0) {
      for (nei_id_num in id_neighbors) {
        match_row <- which(id_outbreak_freq$Id == nei_id_num)
        if (length(match_row) == 0) next
        if (id_outbreak_freq$stage1_freq[match_row] > 0) id_outbreak_freq$stage1_outbreak_nei[i] <- id_outbreak_freq$stage1_outbreak_nei[i] + 1
        if (id_outbreak_freq$stage2_freq[match_row] > 0) id_outbreak_freq$stage2_outbreak_nei[i] <- id_outbreak_freq$stage2_outbreak_nei[i] + 1
        if (id_outbreak_freq$stage3_freq[match_row] > 0) id_outbreak_freq$stage3_outbreak_nei[i] <- id_outbreak_freq$stage3_outbreak_nei[i] + 1
      }
    }
  }

  id_outbreak_freq
}

get_month_weight_tables <- function(eu_aiv_df, sum_stage1_freq, sum_stage2_freq, sum_stage3_freq) {
  stage1_start <- as.Date("2021-05-01")
  stage1_end <- as.Date("2021-10-31")
  stage2_start <- as.Date("2021-11-01")
  stage2_end <- as.Date("2022-04-30")
  stage3_start <- as.Date("2022-05-01")
  stage3_end <- as.Date("2022-10-31")

  aiv_year_month_df <- eu_aiv_df[, c("Event.ID", "observation.date", "observation_year", "observation_month")]
  aiv_year_month_df$observation.date <- as.Date(aiv_year_month_df$observation.date)

  stage1_aiv_month <- aiv_year_month_df %>% dplyr::filter(observation.date >= stage1_start & observation.date <= stage1_end)
  stage2_aiv_month <- aiv_year_month_df %>% dplyr::filter(observation.date >= stage2_start & observation.date <= stage2_end)
  stage3_aiv_month <- aiv_year_month_df %>% dplyr::filter(observation.date >= stage3_start & observation.date <= stage3_end)

  stage1_aiv_month_weighted <- stage1_aiv_month %>%
    dplyr::group_by(observation_year, observation_month) %>%
    dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")
  stage2_aiv_month_weighted <- stage2_aiv_month %>%
    dplyr::group_by(observation_year, observation_month) %>%
    dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")
  stage3_aiv_month_weighted <- stage3_aiv_month %>%
    dplyr::group_by(observation_year, observation_month) %>%
    dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")

  stage1_aiv_month_weighted$month_weighted <- stage1_aiv_month_weighted$sample_count / max(sum_stage1_freq, 1)
  stage2_aiv_month_weighted$month_weighted <- stage2_aiv_month_weighted$sample_count / max(sum_stage2_freq, 1)
  stage3_aiv_month_weighted$month_weighted <- stage3_aiv_month_weighted$sample_count / max(sum_stage3_freq, 1)

  list(
    stage1_start = stage1_start,
    stage1_end = stage1_end,
    stage2_start = stage2_start,
    stage2_end = stage2_end,
    stage3_start = stage3_start,
    stage3_end = stage3_end,
    stage1_aiv_month_weighted = stage1_aiv_month_weighted,
    stage2_aiv_month_weighted = stage2_aiv_month_weighted,
    stage3_aiv_month_weighted = stage3_aiv_month_weighted
  )
}

get_weighted_abundance <- function(abundance, stage_dates, stage1_aiv_month_weighted, stage2_aiv_month_weighted, stage3_aiv_month_weighted) {
  weighted_abundance <- abundance
  weighted_abundance$unweighted_abundance <- weighted_abundance$abundance_filtered_mean

  weighted_abundance <- weighted_abundance %>%
    dplyr::mutate(date = as.Date(paste(year_number, month_number, "01", sep = "-")))

  weighted_abundance %>%
    dplyr::mutate(
      stage = dplyr::case_when(
        date >= stage_dates$stage1_start & date <= stage_dates$stage1_end ~ "stage1",
        date >= stage_dates$stage2_start & date <= stage_dates$stage2_end ~ "stage2",
        date >= stage_dates$stage3_start & date <= stage_dates$stage3_end ~ "stage3",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(stage)) %>%
    dplyr::mutate(
      observation_year = as.numeric(format(date, "%Y")),
      observation_month = as.numeric(format(date, "%m"))
    ) %>%
    dplyr::left_join(
      dplyr::bind_rows(
        stage1_aiv_month_weighted %>% dplyr::mutate(stage = "stage1"),
        stage2_aiv_month_weighted %>% dplyr::mutate(stage = "stage2"),
        stage3_aiv_month_weighted %>% dplyr::mutate(stage = "stage3")
      ),
      by = c("observation_year", "observation_month", "stage")
    ) %>%
    dplyr::mutate(
      month_weighted = dplyr::coalesce(month_weighted, 0),
      weighted_abundance = unweighted_abundance * month_weighted
    ) %>%
    dplyr::group_by(Id, stage) %>%
    dplyr::summarise(sum_weighted_abundance = sum(weighted_abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(Id, stage = c("stage1", "stage2", "stage3"), fill = list(sum_weighted_abundance = 0))
}

merge_all_stages_data <- function(all_df, weighted_abundance_result, eu_center_coordinate_df) {
  stage1_result <- weighted_abundance_result[, c(1, 3)][weighted_abundance_result$stage == "stage1", ]
  stage2_result <- weighted_abundance_result[, c(1, 3)][weighted_abundance_result$stage == "stage2", ]
  stage3_result <- weighted_abundance_result[, c(1, 3)][weighted_abundance_result$stage == "stage3", ]
  colnames(stage1_result)[2] <- "weighted_abundance_stage1"
  colnames(stage2_result)[2] <- "weighted_abundance_stage2"
  colnames(stage3_result)[2] <- "weighted_abundance_stage3"

  all_info_df <- merge(all_df, stage1_result, by = "Id", all = TRUE)
  all_info_df <- merge(all_info_df, stage2_result, by = "Id", all = TRUE)
  all_info_df <- merge(all_info_df, stage3_result, by = "Id", all = TRUE)
  all_info_df <- merge(all_info_df, eu_center_coordinate_df, by = "Id", all = TRUE)
  all_info_df <- all_info_df[!is.na(all_info_df$weighted_abundance_stage1), ]
  all_info_df$stage1_freq[is.na(all_info_df$stage1_freq)] <- 0
  all_info_df$stage2_freq[is.na(all_info_df$stage2_freq)] <- 0
  all_info_df$stage3_freq[is.na(all_info_df$stage3_freq)] <- 0
  all_info_df$stage1_outbreak_nei[is.na(all_info_df$stage1_outbreak_nei)] <- 0
  all_info_df$stage2_outbreak_nei[is.na(all_info_df$stage2_outbreak_nei)] <- 0
  all_info_df$stage3_outbreak_nei[is.na(all_info_df$stage3_outbreak_nei)] <- 0
  all_info_df$standard_chicken_density <- as.numeric(scale(all_info_df$chicken_sum_VALUE))
  all_info_df$standard_duck_density <- as.numeric(scale(all_info_df$duck_sum_VALUE))
  all_info_df
}

rbind_all_stages_data <- function(all_df, weighted_abundance_result, eu_center_coordinate_df) {
  all_stage1_df <- all_df[, c("Id", "chicken_sum_VALUE", "duck_sum_VALUE", "stage1_freq", "stage1_outbreak_binary", "stage1_outbreak_nei")]
  all_stage2_df <- all_df[, c("Id", "chicken_sum_VALUE", "duck_sum_VALUE", "stage2_freq", "stage2_outbreak_binary", "stage2_outbreak_nei")]
  all_stage3_df <- all_df[, c("Id", "chicken_sum_VALUE", "duck_sum_VALUE", "stage3_freq", "stage3_outbreak_binary", "stage3_outbreak_nei")]
  colnames(all_stage1_df)[4:6] <- c("outbreak_freq", "outbreak_binary", "outbreak_nei")
  colnames(all_stage2_df)[4:6] <- c("outbreak_freq", "outbreak_binary", "outbreak_nei")
  colnames(all_stage3_df)[4:6] <- c("outbreak_freq", "outbreak_binary", "outbreak_nei")
  all_stage1_df$stage <- "stage1"
  all_stage2_df$stage <- "stage2"
  all_stage3_df$stage <- "stage3"

  all_stage1_df <- merge(all_stage1_df, weighted_abundance_result, by = c("Id", "stage"), all = FALSE)
  all_stage2_df <- merge(all_stage2_df, weighted_abundance_result, by = c("Id", "stage"), all = FALSE)
  all_stage3_df <- merge(all_stage3_df, weighted_abundance_result, by = c("Id", "stage"), all = FALSE)
  all_stage1_df <- merge(all_stage1_df, eu_center_coordinate_df, by = "Id", all = FALSE)
  all_stage2_df <- merge(all_stage2_df, eu_center_coordinate_df, by = "Id", all = FALSE)
  all_stage3_df <- merge(all_stage3_df, eu_center_coordinate_df, by = "Id", all = FALSE)

  all_stage_df <- rbind(all_stage1_df, all_stage2_df, all_stage3_df)
  all_stage_df$outbreak_freq[is.na(all_stage_df$outbreak_freq)] <- 0
  all_stage_df$outbreak_binary[is.na(all_stage_df$outbreak_binary)] <- 0
  all_stage_df$stage <- factor(all_stage_df$stage, levels = c("stage1", "stage2", "stage3"))
  all_stage_df$standard_chicken_density <- as.numeric(scale(all_stage_df$chicken_sum_VALUE))
  all_stage_df$standard_duck_density <- as.numeric(scale(all_stage_df$duck_sum_VALUE))
  all_stage_df
}

get_single_stage_gam_model <- function(birdname, stage_merge_df) {
  auto_gam_fit <- function(stage_var, data, kval_init = 30, max_k = 150, step_k = 20) {
    kval <- kval_init
    acceptable <- FALSE
    gam_fit <- NULL

    while (!acceptable && kval <= max_k) {
      formula_expr <- as.formula(
        paste0(
          stage_var, "_freq ~ weighted_abundance_", stage_var,
          " + standard_chicken_density + standard_duck_density + ",
          "s(latitude, longitude, k = ", kval, ", bs = 'tp')"
        )
      )

      gam_fit <- tryCatch(
        mgcv::gam(formula = formula_expr, family = poisson(link = "log"), method = "REML", data = data),
        warning = function(w) {
          message(paste0(stage_var, " model warning for bird ", birdname, " at k = ", kval, ":\n", conditionMessage(w)))
          NULL
        },
        error = function(e) {
          message(paste0(stage_var, " model failed for bird ", birdname, " at k = ", kval, ":\n", conditionMessage(e)))
          NULL
        }
      )
      if (is.null(gam_fit)) return(NULL)

      kcheck <- tryCatch(mgcv:::k.check(gam_fit), error = function(e) NULL)
      if (is.null(kcheck)) return(NULL)
      pvals <- kcheck[, "p-value"]
      if (all(pvals >= 0.05)) acceptable <- TRUE else kval <- kval + step_k
    }
    gam_fit
  }

  list(
    gam_stage1 = auto_gam_fit("stage1", stage_merge_df),
    gam_stage2 = auto_gam_fit("stage2", stage_merge_df),
    gam_stage3 = auto_gam_fit("stage3", stage_merge_df)
  )
}

get_all_stage_model <- function(birdname, stage_rbind_df, kval_init = 50, max_k = 150, step_k = 20) {
  kval <- kval_init
  acceptable <- FALSE
  gam_stage4 <- NULL

  while (!acceptable && kval <= max_k) {
    gam_stage4 <- tryCatch(
      mgcv::gam(
        outbreak_freq ~ stage * sum_weighted_abundance + standard_chicken_density + standard_duck_density +
          s(latitude, longitude, k = kval, bs = "tp"),
        family = poisson(link = "log"),
        method = "REML",
        data = stage_rbind_df
      ),
      warning = function(w) {
        message(paste("stage all stage model warning for bird:", birdname, "\n", conditionMessage(w)))
        NULL
      },
      error = function(e) {
        message(paste("stage all stage model failed for bird:", birdname, "\n", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(gam_stage4)) return(NULL)
    kcheck <- tryCatch(mgcv:::k.check(gam_stage4), error = function(e) NULL)
    if (is.null(kcheck)) return(NULL)
    pvals <- kcheck[, "p-value"]
    if (all(pvals >= 0.05)) acceptable <- TRUE else kval <- kval + step_k
  }
  gam_stage4
}

get_all_stages_stage_table <- function(birdname, gam_stage4, num_grid) {
  p_coeff <- p_pv <- rep(NA, 8)
  if (!is.null(gam_stage4)) {
    p_coeff <- summary(gam_stage4)$p.coeff
    p_pv <- summary(gam_stage4)$p.pv
  }
  data.frame(
    bird_name = birdname,
    coef_stagestage1 = p_coeff[2],
    p_stagestage1 = p_pv[2],
    coef_stagestage3 = p_coeff[3],
    p_stagestage3 = p_pv[3],
    num_grid = num_grid
  )
}

get_all_stages_density_table <- function(birdname, gam_stage4, num_grid) {
  p_coeff <- p_pv <- rep(NA, 8)
  if (!is.null(gam_stage4)) {
    p_coeff <- summary(gam_stage4)$p.coeff
    p_pv <- summary(gam_stage4)$p.pv
  }
  data.frame(
    bird_name = birdname,
    coef_chicken_density = p_coeff[5],
    p_chicken_density = p_pv[5],
    coef_duck_density = p_coeff[6],
    p_duck_density = p_pv[6],
    num_grid = num_grid
  )
}

get_all_stages_stageabundance_table <- function(birdname, gam_stage4, num_grid) {
  p_coeff <- p_pv <- rep(NA, 8)
  if (!is.null(gam_stage4)) {
    p_coeff <- summary(gam_stage4)$p.coeff
    p_pv <- summary(gam_stage4)$p.pv
  }
  data.frame(
    bird_name = birdname,
    coef_abundance = p_coeff[4],
    p_abundance = p_pv[4],
    coef_stage1abundance = p_coeff[7],
    p_stage1abundance = p_pv[7],
    coef_stage3abundance = p_coeff[8],
    p_stage3abundance = p_pv[8],
    num_grid = num_grid
  )
}

get_single_stage_abundance_table <- function(birdname, single_stage_gam_model, num_grid) {
  gam_stage1 <- single_stage_gam_model$gam_stage1
  gam_stage2 <- single_stage_gam_model$gam_stage2
  gam_stage3 <- single_stage_gam_model$gam_stage3
  p_coeff_stage1 <- p_pv_stage1 <- rep(NA, 4)
  p_coeff_stage2 <- p_pv_stage2 <- rep(NA, 4)
  p_coeff_stage3 <- p_pv_stage3 <- rep(NA, 4)
  if (!is.null(gam_stage1)) {
    p_coeff_stage1 <- summary(gam_stage1)$p.coeff
    p_pv_stage1 <- summary(gam_stage1)$p.pv
  }
  if (!is.null(gam_stage2)) {
    p_coeff_stage2 <- summary(gam_stage2)$p.coeff
    p_pv_stage2 <- summary(gam_stage2)$p.pv
  }
  if (!is.null(gam_stage3)) {
    p_coeff_stage3 <- summary(gam_stage3)$p.coeff
    p_pv_stage3 <- summary(gam_stage3)$p.pv
  }
  data.frame(
    bird_name = birdname,
    coef_abundance_stage1 = ifelse(length(p_coeff_stage1) >= 2, p_coeff_stage1[2], NA),
    p_abundance_stage1 = ifelse(length(p_pv_stage1) >= 2, p_pv_stage1[2], NA),
    coef_abundance_stage2 = ifelse(length(p_coeff_stage2) >= 2, p_coeff_stage2[2], NA),
    p_abundance_stage2 = ifelse(length(p_pv_stage2) >= 2, p_pv_stage2[2], NA),
    coef_abundance_stage3 = ifelse(length(p_coeff_stage3) >= 2, p_coeff_stage3[2], NA),
    p_abundance_stage3 = ifelse(length(p_pv_stage3) >= 2, p_pv_stage3[2], NA),
    num_grid = num_grid
  )
}

get_single_stage_density_table <- function(birdname, single_stage_gam_model, num_grid) {
  gam_stage1 <- single_stage_gam_model$gam_stage1
  gam_stage2 <- single_stage_gam_model$gam_stage2
  gam_stage3 <- single_stage_gam_model$gam_stage3
  p_coeff_stage1 <- p_pv_stage1 <- rep(NA, 4)
  p_coeff_stage2 <- p_pv_stage2 <- rep(NA, 4)
  p_coeff_stage3 <- p_pv_stage3 <- rep(NA, 4)
  if (!is.null(gam_stage1)) {
    p_coeff_stage1 <- summary(gam_stage1)$p.coeff
    p_pv_stage1 <- summary(gam_stage1)$p.pv
  }
  if (!is.null(gam_stage2)) {
    p_coeff_stage2 <- summary(gam_stage2)$p.coeff
    p_pv_stage2 <- summary(gam_stage2)$p.pv
  }
  if (!is.null(gam_stage3)) {
    p_coeff_stage3 <- summary(gam_stage3)$p.coeff
    p_pv_stage3 <- summary(gam_stage3)$p.pv
  }
  data.frame(
    bird_name = birdname,
    coef_chicken_stage1 = ifelse(length(p_coeff_stage1) >= 4, p_coeff_stage1[3], NA),
    p_chicken_stage1 = ifelse(length(p_pv_stage1) >= 4, p_pv_stage1[3], NA),
    coef_duck_stage1 = ifelse(length(p_coeff_stage1) >= 4, p_coeff_stage1[4], NA),
    p_duck_stage1 = ifelse(length(p_pv_stage1) >= 4, p_pv_stage1[4], NA),
    coef_chicken_stage2 = ifelse(length(p_coeff_stage2) >= 4, p_coeff_stage2[3], NA),
    p_chicken_stage2 = ifelse(length(p_pv_stage2) >= 4, p_pv_stage2[3], NA),
    coef_duck_stage2 = ifelse(length(p_coeff_stage2) >= 4, p_coeff_stage2[4], NA),
    p_duck_stage2 = ifelse(length(p_pv_stage2) >= 4, p_pv_stage2[4], NA),
    coef_chicken_stage3 = ifelse(length(p_coeff_stage3) >= 4, p_coeff_stage3[3], NA),
    p_chicken_stage3 = ifelse(length(p_pv_stage3) >= 4, p_pv_stage3[3], NA),
    coef_duck_stage3 = ifelse(length(p_coeff_stage3) >= 4, p_coeff_stage3[4], NA),
    p_duck_stage3 = ifelse(length(p_pv_stage3) >= 4, p_pv_stage3[4], NA),
    num_grid = num_grid
  )
}

save_outbreak_and_abundance_image <- function(eu_map, stage_rbind_df, stage, birdname, output_folder) {
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  bird_folder <- file.path(output_folder, birdname)
  if (!dir.exists(bird_folder)) dir.create(bird_folder, recursive = TRUE)

  stage_df <- stage_rbind_df[stage_rbind_df$stage == stage, c("Id", "outbreak_freq", "outbreak_binary", "sum_weighted_abundance", "stage")]
  outbreak_abundance_dat <- merge(eu_map, stage_df, by = "Id", all = TRUE)

  outbreak_quartiles <- quantile(outbreak_abundance_dat$outbreak_freq, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  abundance_quartiles <- quantile(outbreak_abundance_dat$sum_weighted_abundance, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

  outbreak_Q1 <- outbreak_quartiles[1]
  outbreak_Q2 <- outbreak_quartiles[2]
  outbreak_Q3 <- outbreak_quartiles[3]
  abundance_Q1 <- abundance_quartiles[1]
  abundance_Q2 <- abundance_quartiles[2]
  abundance_Q3 <- abundance_quartiles[3]

  outbreak_abundance_dat$outbreak_freq[is.na(outbreak_abundance_dat$outbreak_freq)] <- -1
  outbreak_abundance_dat$sum_weighted_abundance[is.na(outbreak_abundance_dat$sum_weighted_abundance)] <- -1

  outbreak_abundance_dat$outbreak_freq_range <- 0
  outbreak_abundance_dat$abundance_range <- 0

  for (i in seq_len(nrow(outbreak_abundance_dat))) {
    if (outbreak_abundance_dat$outbreak_freq[i] == -1) {
      outbreak_abundance_dat$outbreak_freq_range[i] <- 0
    } else if (0 <= outbreak_abundance_dat$outbreak_freq[i] && outbreak_abundance_dat$outbreak_freq[i] <= outbreak_Q1) {
      outbreak_abundance_dat$outbreak_freq_range[i] <- 1
    } else if (outbreak_Q1 < outbreak_abundance_dat$outbreak_freq[i] && outbreak_abundance_dat$outbreak_freq[i] <= outbreak_Q2) {
      outbreak_abundance_dat$outbreak_freq_range[i] <- 2
    } else if (outbreak_Q2 < outbreak_abundance_dat$outbreak_freq[i] && outbreak_abundance_dat$outbreak_freq[i] <= outbreak_Q3) {
      outbreak_abundance_dat$outbreak_freq_range[i] <- 3
    } else {
      outbreak_abundance_dat$outbreak_freq_range[i] <- 4
    }
  }

  for (i in seq_len(nrow(outbreak_abundance_dat))) {
    if (outbreak_abundance_dat$sum_weighted_abundance[i] == -1) {
      outbreak_abundance_dat$abundance_range[i] <- 0
    } else if (0 <= outbreak_abundance_dat$sum_weighted_abundance[i] && outbreak_abundance_dat$sum_weighted_abundance[i] <= abundance_Q1) {
      outbreak_abundance_dat$abundance_range[i] <- 1
    } else if (abundance_Q1 < outbreak_abundance_dat$sum_weighted_abundance[i] && outbreak_abundance_dat$sum_weighted_abundance[i] <= abundance_Q2) {
      outbreak_abundance_dat$abundance_range[i] <- 2
    } else if (abundance_Q2 < outbreak_abundance_dat$sum_weighted_abundance[i] && outbreak_abundance_dat$sum_weighted_abundance[i] <= abundance_Q3) {
      outbreak_abundance_dat$abundance_range[i] <- 3
    } else {
      outbreak_abundance_dat$abundance_range[i] <- 4
    }
  }

  outbreak_label1 <- paste0("minimum-", round(outbreak_Q1, 4))
  outbreak_label2 <- paste0(round(outbreak_Q1, 4), "-", round(outbreak_Q2, 4))
  outbreak_label3 <- paste0(round(outbreak_Q2, 4), "-", round(outbreak_Q3, 4))
  outbreak_label4 <- paste0(round(outbreak_Q3, 4), "-maximum")

  abundance_label1 <- paste0("minimum-", round(abundance_Q1, 4))
  abundance_label2 <- paste0(round(abundance_Q1, 4), "-", round(abundance_Q2, 4))
  abundance_label3 <- paste0(round(abundance_Q2, 4), "-", round(abundance_Q3, 4))
  abundance_label4 <- paste0(round(abundance_Q3, 4), "-maximum")

  outbreak_abundance_dat$outbreak_freq_range <- factor(
    outbreak_abundance_dat$outbreak_freq_range,
    levels = 0:4,
    labels = c("NA", outbreak_label1, outbreak_label2, outbreak_label3, outbreak_label4)
  )

  outbreak_abundance_dat$abundance_range <- factor(
    outbreak_abundance_dat$abundance_range,
    levels = 0:4,
    labels = c("NA", abundance_label1, abundance_label2, abundance_label3, abundance_label4)
  )

  p <- ggplot2::ggplot(data = outbreak_abundance_dat) +
    ggplot2::geom_sf(ggplot2::aes(fill = abundance_range), color = "black", linewidth = 0.1) +
    ggplot2::scale_fill_manual(
      values = c("#FFEEEE", "#FFCCCC", "#FFAAAA", "#FF6666", "#FF0000"),
      name = "abundance range"
    ) +
    ggplot2::labs(
      title = paste0(birdname, " weighted abundance in ", stage),
      caption = "Data source: weighted abundance from run_validation_and_prediction_summary.R"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  output_path <- file.path(bird_folder, paste0(stage, "_plot.png"))
  ggplot2::ggsave(output_path, plot = p, width = 12, height = 6, dpi = 300, device = "png")
}

save_overview_maps <- function(eu_map, chicken_livestock_dat, duck_livestock_dat, eu_aiv_sf, overview_dir, outbreak_type) {
  dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)

  chicken_livestock_dat$sum_VALUE_quartile <- dplyr::ntile(chicken_livestock_dat$sum_VALUE, 4)
  duck_livestock_dat$sum_VALUE_quartile <- dplyr::ntile(duck_livestock_dat$sum_VALUE, 4)

  chicken_quartile_breaks <- quantile(chicken_livestock_dat$sum_VALUE, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  chicken_quartile_labels <- c(
    paste0(formatC(chicken_quartile_breaks[1], format = "d", big.mark = ","), " - ", formatC(chicken_quartile_breaks[2], format = "d", big.mark = ",")),
    paste0(formatC(chicken_quartile_breaks[2] + 1, format = "d", big.mark = ","), " - ", formatC(chicken_quartile_breaks[3], format = "d", big.mark = ",")),
    paste0(formatC(chicken_quartile_breaks[3] + 1, format = "d", big.mark = ","), " - ", formatC(chicken_quartile_breaks[4], format = "d", big.mark = ",")),
    paste0(formatC(chicken_quartile_breaks[4] + 1, format = "d", big.mark = ","), " - ", formatC(chicken_quartile_breaks[5], format = "d", big.mark = ","))
  )

  duck_quartile_breaks <- quantile(duck_livestock_dat$sum_VALUE, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  duck_quartile_labels <- c(
    paste0(formatC(duck_quartile_breaks[1], format = "d", big.mark = ","), " - ", formatC(duck_quartile_breaks[2], format = "d", big.mark = ",")),
    paste0(formatC(duck_quartile_breaks[2] + 1, format = "d", big.mark = ","), " - ", formatC(duck_quartile_breaks[3], format = "d", big.mark = ",")),
    paste0(formatC(duck_quartile_breaks[3] + 1, format = "d", big.mark = ","), " - ", formatC(duck_quartile_breaks[4], format = "d", big.mark = ",")),
    paste0(formatC(duck_quartile_breaks[4] + 1, format = "d", big.mark = ","), " - ", formatC(duck_quartile_breaks[5], format = "d", big.mark = ","))
  )

  chicken_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = chicken_livestock_dat, ggplot2::aes(fill = as.factor(sum_VALUE_quartile)), color = "gray30", linewidth = 0.1) +
    ggplot2::scale_fill_brewer(palette = "Reds", name = "range", labels = chicken_quartile_labels) +
    ggplot2::coord_sf(xlim = c(-25, 70), ylim = c(33, 79), expand = FALSE) +
    ggplot2::theme_dark(base_size = 12) +
    ggplot2::ggtitle("chicken livestock density")
  duck_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = duck_livestock_dat, ggplot2::aes(fill = as.factor(sum_VALUE_quartile)), color = "gray30", linewidth = 0.1) +
    ggplot2::scale_fill_brewer(palette = "Reds", name = "range", labels = duck_quartile_labels) +
    ggplot2::coord_sf(xlim = c(-25, 70), ylim = c(33, 79), expand = FALSE) +
    ggplot2::theme_dark(base_size = 12) +
    ggplot2::ggtitle("duck livestock density")
  outbreak_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = eu_map, fill = "white", color = "black", linewidth = 0.1) +
    ggplot2::geom_sf(data = eu_aiv_sf, color = "red", size = 0.7, alpha = 0.9) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = paste0("AIV Outbreak (", outbreak_type, ") 2021-2022"))

  ggplot2::ggsave(file.path(overview_dir, "chicken_density.png"), chicken_map, width = 10, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(overview_dir, "duck_density.png"), duck_map, width = 10, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(overview_dir, paste0(outbreak_type, "_outbreak_map.png")), outbreak_map, width = 10, height = 6, dpi = 300)
}

save_outbreak_monthly_distribution_plot <- function(all_aiv_df, outbreak_type_df, outbreak_type, overview_dir) {
  build_monthly_count_df <- function(df, label) {
    df$observation.date <- suppressWarnings(lubridate::parse_date_time(
      df$observation.date,
      orders = c(
        "Y/m/d HMSz", "Y/m/d HMS", "Y/m/d",
        "Y-m-d HMSz", "Y-m-d HMS", "Y-m-d",
        "ymd HMSz", "ymd HMS", "ymd"
      ),
      tz = "UTC"
    ))
    df <- df[!is.na(df$observation.date), , drop = FALSE]
    if (nrow(df) == 0) {
      return(data.frame(
        observation_year_month = character(),
        outbreak_count = numeric(),
        series = character(),
        stringsAsFactors = FALSE
      ))
    }
    df$observation_year_month <- format(as.Date(df$observation.date), "%Y-%m")

    monthly_df <- as.data.frame(table(df$observation_year_month), stringsAsFactors = FALSE)
    colnames(monthly_df) <- c("observation_year_month", "outbreak_count")
    monthly_df$series <- label
    monthly_df
  }

  all_monthly_df <- build_monthly_count_df(all_aiv_df, "Wild+Domestic")
  type_monthly_df <- build_monthly_count_df(outbreak_type_df, outbreak_type)

  plot_monthly_df <- function(monthly_df, title_text, output_path) {
    p <- ggplot2::ggplot(monthly_df, ggplot2::aes(x = observation_year_month, y = outbreak_count)) +
      ggplot2::geom_col(fill = "skyblue") +
      ggplot2::labs(
        title = title_text,
        x = "Year-Month",
        y = "Number of Events"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
      )

    ggplot2::ggsave(output_path, p, width = 10, height = 6, dpi = 300)
  }

  plot_monthly_df(
    all_monthly_df,
    "Monthly Distribution of Avian Influenza Events (Wild+Domestic, 2021-2022)",
    file.path(overview_dir, "Wild_Domestic_monthly_outbreak_distribution.png")
  )
  plot_monthly_df(
    type_monthly_df,
    paste0("Monthly Distribution of Avian Influenza Events (", outbreak_type, ", 2021-2022)"),
    file.path(overview_dir, paste0(outbreak_type, "_monthly_outbreak_distribution.png"))
  )
}

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_name <- if (length(script_file_arg) > 0) basename(sub("^--file=", "", script_file_arg[1])) else "run_aiv_analysis.R"
opts <- parse_args()

if (opts$help) {
  print_usage(script_name)
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(mgcv)
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(lubridate)
  library(tidyr)
})

validate_options(opts)
species_catalog <- build_species_catalog(opts$bird_abundance_folder)

if (opts$list_species) {
  write.table(species_catalog["scientific_name"], row.names = FALSE, col.names = TRUE, quote = FALSE)
  quit(save = "no", status = 0)
}

selected_species <- resolve_species(opts, species_catalog)

dir.create(opts$output_folder, recursive = TRUE, showWarnings = FALSE)
aiv_analysis_root <- file.path(opts$output_folder, "aiv_analysis")
aiv_analysis_type_folder <- file.path(aiv_analysis_root, opts$outbreak_type)
aiv_analysis_output_folder <- file.path(aiv_analysis_type_folder, opts$write_date)
weighted_abundance_output_folder <- file.path(aiv_analysis_output_folder, "weighted_abundance")
overview_maps_output_folder <- file.path(aiv_analysis_output_folder, "overview_maps")
dir.create(aiv_analysis_root, recursive = TRUE, showWarnings = FALSE)
dir.create(aiv_analysis_type_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(aiv_analysis_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(weighted_abundance_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(overview_maps_output_folder, recursive = TRUE, showWarnings = FALSE)

eu_map <- sf::st_read(opts$eu_shp_path, quiet = TRUE)
eu_map_df <- as.data.frame(eu_map)
neighbors_df <- get_neighbors_df(eu_map)
eu_center_coordinate_df <- get_eu_center_coordinate_df(eu_map)

chicken_livestock_sum_density_df <- get_livestock_sum_density(opts$chicken_density_path)
duck_livestock_sum_density_df <- get_livestock_sum_density(opts$duck_density_path)
chicken_livestock_dat <- merge(eu_map, chicken_livestock_sum_density_df, by = "Id", all = TRUE)
duck_livestock_dat <- merge(eu_map, duck_livestock_sum_density_df, by = "Id", all = TRUE)
chicken_livestock_dat$sum_VALUE[is.na(chicken_livestock_dat$sum_VALUE)] <- 0
duck_livestock_dat$sum_VALUE[is.na(duck_livestock_dat$sum_VALUE)] <- 0

aiv_2021 <- get_aiv_records(opts$aiv_2021_path, opts$outbreak_type, sf::st_crs(eu_map))
aiv_2022 <- get_aiv_records(opts$aiv_2022_path, opts$outbreak_type, sf::st_crs(eu_map))
all_aiv_2021_raw <- get_aiv_raw_records(opts$aiv_2021_path)
all_aiv_2022_raw <- get_aiv_raw_records(opts$aiv_2022_path)
eu_aiv_2021_df <- aiv_2021$df
eu_aiv_2022_df <- aiv_2022$df
eu_aiv_2021to2022_sf <- rbind(aiv_2021$sf, aiv_2022$sf)
eu_aiv_df <- rbind(eu_aiv_2022_df, eu_aiv_2021_df)
all_aiv_df <- rbind(all_aiv_2021_raw, all_aiv_2022_raw)
all_aiv_df <- all_aiv_df[all_aiv_df$Species %in% c("Wild", "Domestic"), , drop = FALSE]
eu_aiv_df$observation.date <- lubridate::ymd_hms(eu_aiv_df$observation.date)
eu_aiv_df <- eu_aiv_df[!is.na(eu_aiv_df$observation.date), ]
eu_aiv_df$observation_year <- lubridate::year(eu_aiv_df$observation.date)
eu_aiv_df$observation_month <- lubridate::month(eu_aiv_df$observation.date)

id_outbreak_freq <- build_outbreak_frequency(eu_map_df, eu_aiv_df, neighbors_df, opts$outbreak_type)
sum_stage1_freq <- sum(id_outbreak_freq$stage1_freq)
sum_stage2_freq <- sum(id_outbreak_freq$stage2_freq)
sum_stage3_freq <- sum(id_outbreak_freq$stage3_freq)
stage_dates <- get_month_weight_tables(eu_aiv_df, sum_stage1_freq, sum_stage2_freq, sum_stage3_freq)

id_df <- data.frame(Id = eu_map_df$Id)
chicken_livestock_df <- data.frame(chicken_livestock_dat)[, c("Id", "sum_VALUE")]
duck_livestock_df <- data.frame(duck_livestock_dat)[, c("Id", "sum_VALUE")]
names(chicken_livestock_df)[2] <- "chicken_sum_VALUE"
names(duck_livestock_df)[2] <- "duck_sum_VALUE"

all_df <- merge(id_df, chicken_livestock_df, by = "Id", all = TRUE)
all_df <- merge(all_df, duck_livestock_df, by = "Id", all = TRUE)
all_df <- merge(all_df, id_outbreak_freq, by = "Id", all = TRUE)
all_df$stage1_outbreak_binary <- ifelse(all_df$stage1_freq > 0, 1, 0)
all_df$stage2_outbreak_binary <- ifelse(all_df$stage2_freq > 0, 1, 0)
all_df$stage3_outbreak_binary <- ifelse(all_df$stage3_freq > 0, 1, 0)

save_overview_maps(
  eu_map,
  chicken_livestock_dat,
  duck_livestock_dat,
  eu_aiv_2021to2022_sf,
  overview_maps_output_folder,
  opts$outbreak_type
)
save_outbreak_monthly_distribution_plot(
  all_aiv_df = all_aiv_df,
  outbreak_type_df = eu_aiv_df,
  outbreak_type = opts$outbreak_type,
  overview_dir = overview_maps_output_folder
)

all_birds_all_stages_stage_table <- data.frame()
all_birds_all_stages_density_table <- data.frame()
all_birds_all_stages_stageabundance_table <- data.frame()
all_birds_single_stage_abundance_table <- data.frame()
all_birds_single_stage_density_table <- data.frame()

for (birdname in selected_species) {
  cat("Processing species:", birdname, "\n")
  abundance_path <- file.path(opts$bird_abundance_folder, paste0(birdname, ".csv"))
  abundance <- read.csv(abundance_path, stringsAsFactors = FALSE)
  weighted_abundance_result <- get_weighted_abundance(
    abundance,
    stage_dates = stage_dates,
    stage1_aiv_month_weighted = stage_dates$stage1_aiv_month_weighted,
    stage2_aiv_month_weighted = stage_dates$stage2_aiv_month_weighted,
    stage3_aiv_month_weighted = stage_dates$stage3_aiv_month_weighted
  )

  stage_merge_df <- merge_all_stages_data(all_df, weighted_abundance_result, eu_center_coordinate_df)
  stage_rbind_df <- rbind_all_stages_data(all_df, weighted_abundance_result, eu_center_coordinate_df)
  num_grid <- nrow(stage_merge_df)
  stage_rbind_df$stage <- stats::relevel(stage_rbind_df$stage, ref = "stage2")

  gam_stage4 <- get_all_stage_model(birdname, stage_rbind_df, kval_init = 50, max_k = 150, step_k = 20)
  single_stage_gam_model <- get_single_stage_gam_model(birdname, stage_merge_df)

  all_birds_single_stage_abundance_table <- rbind(
    all_birds_single_stage_abundance_table,
    get_single_stage_abundance_table(birdname, single_stage_gam_model, num_grid)
  )
  all_birds_single_stage_density_table <- rbind(
    all_birds_single_stage_density_table,
    get_single_stage_density_table(birdname, single_stage_gam_model, num_grid)
  )
  all_birds_all_stages_stage_table <- rbind(
    all_birds_all_stages_stage_table,
    get_all_stages_stage_table(birdname, gam_stage4, num_grid)
  )
  all_birds_all_stages_density_table <- rbind(
    all_birds_all_stages_density_table,
    get_all_stages_density_table(birdname, gam_stage4, num_grid)
  )
  all_birds_all_stages_stageabundance_table <- rbind(
    all_birds_all_stages_stageabundance_table,
    get_all_stages_stageabundance_table(birdname, gam_stage4, num_grid)
  )

  save_outbreak_and_abundance_image(eu_map, stage_rbind_df, "stage1", birdname, weighted_abundance_output_folder)
  save_outbreak_and_abundance_image(eu_map, stage_rbind_df, "stage2", birdname, weighted_abundance_output_folder)
  save_outbreak_and_abundance_image(eu_map, stage_rbind_df, "stage3", birdname, weighted_abundance_output_folder)
}

rownames(all_birds_all_stages_stage_table) <- NULL
rownames(all_birds_all_stages_density_table) <- NULL
rownames(all_birds_all_stages_stageabundance_table) <- NULL
rownames(all_birds_single_stage_abundance_table) <- NULL
rownames(all_birds_single_stage_density_table) <- NULL

write.csv(
  all_birds_all_stages_stage_table,
  file.path(aiv_analysis_output_folder, paste0(opts$outbreak_type, "_all_birds_all_stages_stage_", opts$write_date, ".csv")),
  row.names = FALSE
)
write.csv(
  all_birds_all_stages_density_table,
  file.path(aiv_analysis_output_folder, paste0(opts$outbreak_type, "_all_birds_all_stages_density_", opts$write_date, ".csv")),
  row.names = FALSE
)
write.csv(
  all_birds_all_stages_stageabundance_table,
  file.path(aiv_analysis_output_folder, paste0(opts$outbreak_type, "_all_birds_all_stages_stageabundance_", opts$write_date, ".csv")),
  row.names = FALSE
)
write.csv(
  all_birds_single_stage_abundance_table,
  file.path(aiv_analysis_output_folder, paste0(opts$outbreak_type, "_all_birds_single_stage_abundance_", opts$write_date, ".csv")),
  row.names = FALSE
)
write.csv(
  all_birds_single_stage_density_table,
  file.path(aiv_analysis_output_folder, paste0(opts$outbreak_type, "_all_birds_single_stage_density_", opts$write_date, ".csv")),
  row.names = FALSE
)

cat("Done.\n")
cat("AIV analysis folder:", aiv_analysis_output_folder, "\n")
cat("Weighted abundance plot folder:", weighted_abundance_output_folder, "\n")
cat("Overview maps folder:", overview_maps_output_folder, "\n")
