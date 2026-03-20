# 1. Load Package ----------
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


# 2-1. Input path ----------
EU_shp_path <- "D:/aiv_project_2025/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp"
env_folder <- "D:/aiv_project_2025/gee data/EU_combine_avg_result"
lc_path <- "D:/aiv_project_2025/gee data/2016_to_2022_land_cover_interporlation_by_xgboost.csv"
checklist_folder <- "D:/aiv_project_2025/ebird filtered checklist2/"


# 2-2. Output path ----------
CCI_output_folder <- "D:/aiv_project_table/glmm_performance"
validtion_data_folder <- "D:/aiv_project_table/validation_data/"
abundance_model_output_folder <- "D:/aiv_project_table/abundance_spatiotemporal_sampling_method/"


# 3-1. Define fuction part1 (data reading and pre-merging processing) ----------
# read EU shp file and calculate grid neighbors
get_EU_map_and_neighbors <- function(shp_path) {
  EU.map <- st_read(shp_path, quiet = TRUE)
  EU_map <- as.data.frame(EU.map)
  neighbors <- st_touches(EU.map)
  
  # create adjacency matrix
  matrix_data <- matrix(nrow = length(EU.map$Id), ncol = 9)
  colnames(matrix_data) <- c('id', paste0('neighbor', 1:8))
  matrix_data[, 1] <- EU.map$Id
  
  for (i in 1:length(EU.map$Id)) {
    neighbors_indices <- neighbors[[i]]
    if (length(neighbors_indices) > 0) {
      neighbors_id <- EU.map[neighbors_indices, ]
      for (j in 1:min(length(neighbors_id$Id), 8)) {
        matrix_data[i, j+1] <- neighbors_id$Id[j]
      }
    }
  }
  neighbors_df <- as.data.frame(matrix_data)
  return(list(EU.map = EU.map, EU_map = EU_map, neighbors_df = neighbors_df))
}

# read environmental data
load_environmental_data <- function(start_year, end_year, base_path) {
  en_df <- NULL
  
  for (year_num in start_year:end_year) {
    path <- paste0(base_path, "/", year_num, "_median_combined_result.csv")
    part_df <- read.csv(path)
    
    if (is.null(en_df)) {
      en_df <- part_df
    } else {
      en_df <- rbind(en_df, part_df)
    }
  }
  
  # separate month and year
  en_df <- en_df %>%
    separate(Month, into = c("year_number", "month_number"), sep = "-", convert = TRUE)
  
  en_df$month_number <- as.factor(en_df$month_number)
  en_df$year_number <- as.factor(en_df$year_number)
  
  return(en_df)
}

# read land cover data
load_land_cover_data <- function(file_path, year_range) {
  
  land_cover_df <- read.csv(file_path)
  
  land_cover_df <- land_cover_df %>%
    arrange(Id, year_number, month_number)
  
  land_cover_df$month_number <- as.factor(land_cover_df$month_number)
  land_cover_df$year_number <- as.factor(land_cover_df$year_number)
  
  land_cover_df <- land_cover_df %>%
    filter(year_number %in% year_range)
  
  return(land_cover_df)
}


# 3-2. Load shp file, environmental data, and land cover data ----------
EU_map_info_list <- get_EU_map_and_neighbors(EU_shp_path)
EU.map <- EU_map_info_list$EU.map
EU_map <- EU_map_info_list$EU_map
neighbors_df <- EU_map_info_list$neighbors_df

en_df <- load_environmental_data(
  start_year = 2021,
  end_year = 2022,
  base_path = env_folder
)

lc_df <- load_land_cover_data(
  file_path =lc_path,
  year_range = c('2021', '2022')
)


# 4-1. Define fuction part2 (data filtering and abundance model process settings) ----------
#1.
process_observation_data <- function(data_Id, quantile_threshold = 0.99, year_range, protocol_range) {
  
  checklist_feature_col <- c("observation_count", "locality_type", "observation_date", "observer_id",
                             "protocol_type", "duration_mins", "effort_km", "Id", "month_number")
  #data_Id <- oringinal_data_Id
  data_Id <- data_Id[, checklist_feature_col]
  
  data_Id$effort_km[is.na(data_Id$effort_km)] <- 0
  data_Id$duration_mins[is.na(data_Id$duration_mins)] <- 0
  
  #year_range <- c(2021, 2022)
  data_Id <- data_Id %>%
    mutate(year_number = year(as.Date(observation_date)))
  data_Id <- data_Id %>%
    filter(year_number %in% year_range)
  
  X_count <- nrow(data_Id[which(data_Id$observation_count == "X"), ])
  X_proportion <- X_count / nrow(data_Id)
  
  data_Id$observation_count <- ifelse(data_Id$observation_count == "X", "0", data_Id$observation_count)
  data_Id$observation_count <- as.numeric(data_Id$observation_count)
  
  # Keep only specific observation modes
  #protocol_range <- c("Traveling" , "Stationary")
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
  
  # remove outliers with too high number of observations
  ob_count_threshold <- quantile(data_Id$observation_count, probs = quantile_threshold, na.rm = TRUE)
  print(paste0(quantile_threshold*100, '% quantile = ', ob_count_threshold))
  data_Id <- data_Id %>%
    filter(observation_count <= ob_count_threshold)
  
  return(list(
    X_proportion = X_proportion,
    observation_summary = obs_summary,
    ob_count_threshold = ob_count_threshold,
    num_row_data = num_row_data,
    filtered_data = data_Id
  ))
}

#2.
filter_observers_by_quantile <- function(data_Id, X = 0.7) {
  
  # count the number of observations for each observer
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
  
  # statistical chart parameter settings
  retained_plot <- ggplot(retained_stats, aes(x = Percentile, y = Retained_Ratio)) +
    geom_line(group = 1, color = "steelblue", size = 1.2) +
    geom_point(size = 3, color = "darkred") +
    geom_text(aes(label = paste0(round(Retained_Ratio * 100, 1), "%")),
              vjust = -0.5, size = 4.5) +
    labs(
      title = "sample retention ratio vs observer observation number quantile threshold",
      x = "quantile threshold of observer observation times",
      y = "sample retention ratio"
    ) +
    ylim(0, 1) +
    theme_minimal()
  
  # get the filter threshold for X quantile
  threshold_X <- quantile(observer_freq$observation_times, probs = X)
  
  top_observers_X <- observer_freq %>%
    filter(observation_times >= threshold_X) %>%
    pull(observer_id)
  
  data_top_X <- data_Id %>%
    filter(observer_id %in% top_observers_X)
  
  retained_ratio <- nrow(data_top_X) / nrow(data_Id)
  
  return(list(
    plot = retained_plot,
    threshold = threshold_X,
    retained_ratio = retained_ratio,
    data_top_X = data_top_X
  ))
}

#3.
filter_grids_with_plot <- function(data_Id, quantile_threshold = 0.5) {
  
  # count the number of observations for each geographic grid
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
  
  # statistical chart parameter settings
  retained_plot <- ggplot(retained_stats, aes(x = Percentile, y = Retained_Ratio)) +
    geom_line(group = 1, color = "steelblue", size = 1.2) +
    geom_point(size = 3, color = "darkred") +
    geom_text(aes(label = paste0(round(Retained_Ratio * 100, 1), "%")),
              vjust = -0.5, size = 4.5) +
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
    rename(id_observation_times = observation_times)  # new variable: total number of observations for each grid(Id)
  
  # calculate statistics for filter results
  num_grids_above <- length(retained_ids)
  total_grids <- nrow(id_freq)
  proportion_above <- num_grids_above / total_grids
  
  cat("Quantile threshold =", quantile_threshold, "\n")
  cat("Observation count threshold =", threshold, "\n")
  cat("The number of grids that meet the threshold =", num_grids_above, "/", total_grids, "\n")
  cat("The ratio =", round(proportion_above, 4), "(", round(proportion_above * 100, 2), "%¡^\n")
  
  return(list(
    threshold = threshold,
    valid_grid_count = num_grids_above,
    total_grids = total_grids,
    proportion = proportion_above,
    filtered_data = filtered_data,
    plot = retained_plot
  ))
}

#4.
calculate_observer_cci <- function(data_Id, birdname, cci_output_folder) {
  #data_Id$duration_mins_scaled <- scale(data_Id$duration_mins)
  #data_Id$effort_km_scaled <- scale(data_Id$effort_km)
  
  cci_randomeffect_output_folder <- paste0(cci_output_folder, '/cci_random_effect')
  glmm_performance_output_folder <- paste0(cci_output_folder, '/glmm_performance')
  
  if (!dir.exists(cci_randomeffect_output_folder)) {
    dir.create(cci_randomeffect_output_folder, recursive = TRUE)
  }
  if (!dir.exists(glmm_performance_output_folder)) {
    dir.create(glmm_performance_output_folder, recursive = TRUE)
  }
  birdname=bird_name
  
  # log-transform the effort variable
  data_Id <- data_Id %>%
    mutate(
      duration_mins_log = log(duration_mins + 1),
      effort_km_log = log(effort_km + 1),
      observer_id = as.factor(observer_id)
    )
  colnames(data_Id)
  
  # perform different glmm modeling based on protocol_type ("Traveling", "Stationary")
  if (all(c(1,0) %in% unique(data_Id$protocol_type))) {
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
  
  # compare AIC/BIC
  cat("AIC:\n")
  aic_values <- AIC(model_poisson, model_nbinom)
  cat("BIC:\n")
  bic_values <- BIC(model_poisson, model_nbinom)
  
  cat("=== AIC performance ===\n")
  print(aic_values)
  cat("=== BIC performance===\n")
  print(bic_values)
  
  glmm_performance_df <- data.frame(birdname = birdname,
                                    AIC_poisson = aic_values$AIC[1],
                                    AIC_nbinom = aic_values$AIC[2],
                                    BIC_poisson = bic_values$BIC[1],
                                    BIC_nbinom =bic_values$BIC[2])
  
  write.csv(glmm_performance_df,
            file.path(glmm_performance_output_folder, paste0(birdname, '.csv')),
            row.names = FALSE)
  
  # select the random effect of Negative Binomial as CCI
  random_effects <- ranef(model_nbinom)
  observer_effects <- random_effects$cond$observer_id
  
  observer_effects_df <- data.frame(
    observer_id = rownames(observer_effects),
    cci = scale(observer_effects[, "(Intercept)"]),
    row.names = NULL
  )
  
  write.csv(observer_effects_df,
            file = file.path(cci_randomeffect_output_folder, paste0(birdname, '.csv')),
            row.names = FALSE)
  
  cat("=== Observer CCI summary ===\n")
  cat("Median:", median(observer_effects_df$cci), "\n")
  cat("Mean:", mean(observer_effects_df$cci), "\n")
  cat("SD:", sd(observer_effects_df$cci), "\n")
  
  # merge back to original data
  data_Id$observer_id <- as.character(data_Id$observer_id)
  data_Id <- merge(data_Id, observer_effects_df, by = "observer_id", all.x = TRUE)
  
  return(data_Id)
}

#6.
# compute Poisson deviance explained (pseudo R^2) to evaluate count model performance
calculate_poisson_deviance_explained <- function(y_true, y_pred) {
  
  # avoid numerical issues (e.g., log(0)) by enforcing lower bound
  y_pred_safe <- pmax(y_pred, 1e-15)
  y_mean <- mean(y_true)
  y_mean_safe <- pmax(y_mean, 1e-15)
  
  
  # calculate model deviance
  # special handling when y_true == 0 (log term becomes 0)
  D_model <- 2 * sum(ifelse(y_true == 0,
                            y_pred_safe,
                            y_true * log(y_true / y_pred_safe) - (y_true - y_pred_safe)),
                     na.rm = TRUE)
  
  # calculate null deviance (baseline model using mean)
  D_null <- 2 * sum(ifelse(y_true == 0,
                           y_mean_safe,
                           y_true * log(y_true / y_mean_safe) - (y_true - y_mean_safe)),
                    na.rm = TRUE)
  
  return(1 - (D_model / D_null))
}


#7-2.
spatiotemporal_sampling <- function(data_Id, validation_data, unique_data_Id, en_df, lc_df,
                                    n_iter, seed = 123, bird_name, output_folder) {
  #seed<-123
  set.seed(seed)
  
  #n_iter<-3
  train_rmse_list <- numeric(n_iter)
  oob_rmse_list <- numeric(n_iter)
  train_size_list <- numeric(n_iter)
  oob_size_list <- numeric(n_iter)
  train_spearman_list <- numeric(n_iter)
  oob_spearman_list <- numeric(n_iter)
  train_pde_list <- numeric(n_iter)
  oob_pde_list <- numeric(n_iter)
  
  # create output folder
  abundance_output_folder <- paste0(output_folder, 'abundance_prediction/', bird_name)
  performance_output_folder <- paste0(output_folder, 'model_performance')
  validation_output_folder <- paste0(output_folder, 'validation_prediction/', bird_name)
  importance_output_folder <- paste0(output_folder, 'feature_importance/', bird_name)
  
  if (!dir.exists(abundance_output_folder)) {
    dir.create(abundance_output_folder, recursive = TRUE)
  }
  if (!dir.exists(performance_output_folder)) {
    dir.create(performance_output_folder, recursive = TRUE)
  }
  if (!dir.exists(validation_output_folder)) {
    dir.create(validation_output_folder, recursive = TRUE)
  }
  if (!dir.exists(importance_output_folder)) {
    dir.create(importance_output_folder, recursive = TRUE)
  }
  
  for (j in 1:12) validation_data[[paste0("month_number", j)]] <- ifelse(validation_data$month_number == j, 1, 0)
  for (j in 2021:2022) validation_data[[paste0("year_number", j)]] <- ifelse(validation_data$year_number == j, 1, 0)
  
  for (i in 1:n_iter) {
    #i<-1
    cat("\nIteration", i, "\n")
    
    # perform stratified sampling while preserving the original zero proportion
    # step 1: calculate zero proportion in the full data set
    zero_proportion <- sum(data_Id$observation_count == 0) / nrow(data_Id)
    
    # step 2: define target sample size (balanced across years and Ids)
    sample_size <- ceiling(mean(table(data_Id$year_number)))
    sample_Id_size <- ceiling(sample_size/length(table(data_Id$Id)))
    
    # step 3: Stratified sampling by Id and year (with replacement if needed)
    Id_list <- table(data_Id$Id)
    length(Id_list)
    sampled_data <- data_Id %>%
      filter(Id %in% names(Id_list)) %>%  
      group_by(Id, year_number) %>%       
      group_modify(~ {
        if (nrow(.x) == 0) return(tibble())                                     # skip empty groups
        slice_sample(.x, n = sample_Id_size, replace = TRUE)                    # sample within each (Id, year) group
      }) %>%
      ungroup() 
    
    
    # step 4: rebalance sampled data to maintain original zero proportion
    sample_zero_count <- ceiling(length(sampled_data$observation_count)*zero_proportion)
    sample_nonzero_count <- length(sampled_data$observation_count) - sample_zero_count
    
    # sample zero observations
    zero_sample <- sampled_data %>%
      filter(observation_count == 0) %>%
      slice_sample(n = sample_zero_count, replace = TRUE)
    
    # sample non-zero observations
    nonzero_sample <- sampled_data %>%
      filter(observation_count > 0) %>%
      slice_sample(n = sample_nonzero_count, replace = TRUE)
    
    # combine final dataset with preserved zero/non-zero ratio
    sampled_data <- rbind(zero_sample, nonzero_sample)
    
    
    for (j in 1:12) sampled_data[[paste0("month_number", j)]] <- ifelse(sampled_data$month_number == j, 1, 0)
    for (j in 2021:2022) sampled_data[[paste0("year_number", j)]] <- ifelse(sampled_data$year_number == j, 1, 0)
    
    #, "protocol_type"
    if (all(c(0,1) %in% unique(data_Id$protocol_type))){
      exclude_cols <- c("observation_date", "observation_count", "Id", "observer_id", "binary_count",
                        "year_number", "month_number", 
                        "duration_mins_log", "effort_km_log", "sample_key")
    }else{
      exclude_cols <- c("observation_date", "observation_count", "Id", "observer_id", "binary_count",
                        "year_number", "month_number", 
                        "duration_mins_log", "effort_km_log", "sample_key", "protocol_type")
    }
    
    feature_cols <- setdiff(colnames(data_Id), exclude_cols)
    filtered_data <- as.data.frame(sampled_data)[, feature_cols]
    X_train <- model.matrix(~ . -1, data = filtered_data)
    y_train <- sampled_data$observation_count
    dtrain <- xgb.DMatrix(data = X_train, label = y_train)
    
    model <- xgboost(data = dtrain, objective = "count:poisson", eta = 0.3, nrounds = 500,
                     eval_metric = "rmse", colsample_bytree = 0.8, nthread = 12, verbose = 0)
    importance_matrix <- xgb.importance(model = model)
    write.csv(importance_matrix,
              file.path(importance_output_folder, paste0('iteration', i, '.csv')),
              row.names = FALSE)
    
    full_data <- data_Id
    for (j in 1:12) full_data[[paste0("month_number", j)]] <- ifelse(full_data$month_number == j, 1, 0)
    for (j in 2021:2022) full_data[[paste0("year_number", j)]] <- ifelse(full_data$year_number == j, 1, 0)
    
    full_X <- model.matrix(~ . -1, data = as.data.frame(full_data)[, feature_cols])
    pred <- predict(model, full_X)
    
    dvalid <- model.matrix(~ . -1, data = as.data.frame(validation_data)[, feature_cols])
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
    write.csv(valid_pred_df,
              file.path(validation_output_folder, sprintf("abundance_iteration%d.csv", i)),
              row.names = FALSE)
    
    
    residual <- pred - full_data$observation_count
    full_data$set <- ifelse(full_data$sample_key %in% sampled_data$sample_key, "train", "oob")
    
    train_rmse <- sqrt(mean(residual[full_data$set == "train"]^2, na.rm = TRUE))
    oob_rmse   <- sqrt(mean(residual[full_data$set == "oob"]^2, na.rm = TRUE))
    valid_rmse   <- sqrt(mean(valid_residual^2, na.rm = TRUE))
    
    train_spearman <- cor(pred[full_data$set == "train"],
                          full_data$observation_count[full_data$set == "train"],
                          method = "spearman")
    oob_spearman   <- cor(pred[full_data$set == "oob"],
                          full_data$observation_count[full_data$set == "oob"],
                          method = "spearman")
    valid_spearman   <- cor(valid_pred,
                            validation_data$observation_count,
                            method = "spearman")
    
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
    
    
    cat("Train RMSE:", round(train_rmse, 4), " | OOB RMSE:", round(oob_rmse, 4),
        " | validation RMSE:", round(valid_rmse, 4), "\n")
    cat("Train size:", sum(full_data$set == "train"), " | OOB size:", sum(full_data$set == "oob"),
        " | validation size:", nrow(validation_data), "\n")
    cat("Train Spearman:", round(train_spearman, 4), " | OOB Spearman:", round(oob_spearman, 4),
        " | validation Spearman:", round(valid_spearman, 4), "\n")
    cat("Train P-DE:", round(train_pde_list[i], 4), " | OOB P-DE:", round(oob_pde_list[i], 4),
        " | validation P-DE:", round(valid_pde, 4), "\n")
    
    train_rmse_list[i] <- train_rmse
    oob_rmse_list[i] <- oob_rmse
    train_size_list[i] <- sum(full_data$set == "train")
    oob_size_list[i] <- sum(full_data$set == "oob")
    train_spearman_list[i] <- train_spearman
    oob_spearman_list[i] <- oob_spearman
    
    # predict 12 months
    Id_list1 <- unique(unique_data_Id$Id)
    abundance_sublist <- vector("list", length(Id_list1))
    
    for (k in seq_along(Id_list1)) {
      gid <- Id_list1[k]
      #print(gid)
      gid1 <- unique_data_Id$id_observation_times[unique_data_Id$Id == gid][1]
      #print(gid1)
      
      predict_data <- expand.grid(Id = gid, year_number = 2021:2022, month_number = 1:12,
                                  id_observation_times = gid1)
      #print(head(predict_data))
      predict_data$locality_type <- 1
      predict_data$protocol_type <- 1
      predict_data$duration_mins <- 60
      predict_data$effort_km <- 2
      predict_data$cci <- 0.5
      
      for (j in 1:12) predict_data[[paste0("month_number", j)]] <- ifelse(predict_data$month_number == j, 1, 0)
      for (j in 2021:2022) predict_data[[paste0("year_number", j)]] <- ifelse(predict_data$year_number == j, 1, 0)
      
      predict_data <- merge(predict_data, en_df, by = c("Id", "year_number", "month_number"), all = FALSE)
      predict_data <- merge(predict_data, lc_df, by = c("Id", "year_number", "month_number"), all = FALSE)
      
      predict_data_X <- predict_data[, intersect(feature_cols, colnames(predict_data))]
      pred_features <- model.matrix(~ . -1, data = predict_data_X)
      dpred <- xgb.DMatrix(data = pred_features)
      predictions_poisson <- predict(model, dpred)
      
      
      result_df <- data.frame(Id = predict_data$Id,
                              year_number = predict_data$year_number,
                              month_number = predict_data$month_number,
                              abundance = predictions_poisson)
      abundance_sublist[[k]] <- result_df
    }
    
    write.csv(bind_rows(abundance_sublist),
              file.path(abundance_output_folder, sprintf("abundance_iteration%d.csv", i)),
              row.names = FALSE)
  }
  
  result_df <- data.frame(
    iteration = 1:length(train_rmse_list),
    train_rmse = train_rmse_list,
    oob_rmse = oob_rmse_list,
    train_size = train_size_list,
    oob_size = oob_size_list,
    train_spearman_list = train_spearman_list,
    oob_spearman_list = oob_spearman_list,
    train_pde_list = train_pde_list,
    oob_pde_list = oob_pde_list
  )
  
  write.csv(result_df, 
            file.path(performance_output_folder, paste0(bird_name, '.csv')),
            row.names = FALSE)
  return(result_df)
}


# 4-2. Start species abundance model training and prediction (each bird) ----------

csv_paths <- list.files(
  path = checklist_folder,
  pattern = "\\.csv$",
  full.names = TRUE
)

for (path_csv in csv_paths){
  
  #path_csv <- csv_paths[1]
  print('Load species checklist...')
  oringinal_data_Id <- read.csv(path_csv)
  
  bird_name <- oringinal_data_Id$scientific_name[1]
  print(paste0('species name : ', bird_name))
  
  env_vars <- colnames(en_df)[grepl("_MEDIAN$", colnames(en_df))]
  lc_vars <- colnames(lc_df)[grepl("_MEAN$", colnames(lc_df))]
  en_df[env_vars] <- scale(en_df[env_vars])
  lc_df[lc_vars] <- scale(lc_df[lc_vars])
  
  print('Data filtering step1 ...')
  filtered_result <- process_observation_data(oringinal_data_Id, quantile_threshold = 0.99,     #Use 99% threshold
                                              year_range = c(2021, 2022),
                                              protocol_range = c("Traveling"))                  #"Stationary"
  #filtered_result$X_proportion
  #filtered_result$observation_summary
  #filtered_result$ob_count_threshold
  num_row_data <- filtered_result$num_row_data
  data_Id <- filtered_result$filtered_data
  print(table(data_Id$protocol_type))
  
  print('Data filtering step2 ...')
  filtered_result <- filter_observers_by_quantile(data_Id, X = 0.7)             #Use 70% threshold
  #filtered_result$threshold
  #filtered_result$retained_ratio
  data_Id <- filtered_result$data_top_X
  print(filtered_result$plot)
  
  print('Data filtering step3 ...')
  filtered_result <- filter_grids_with_plot(data_Id, quantile_threshold = 0.5)  #Use 50% threshold
  #filtered_result$threshold
  #filtered_result$proportion
  data_Id <- filtered_result$filtered_data
  print(filtered_result$plot)
  
  
  print(paste0('The ratio of the number of filtered samples to the number of original samples = ',
               nrow(data_Id)/num_row_data))
  print('Data filtering completed')
  
  
  # convert variable encoding
  data_Id$locality_type <- ifelse(data_Id$locality_type == "H", 1, 0)
  data_Id$binary_count <- ifelse(data_Id$observation_count > 0, 1, 0)
  data_Id$protocol_type <- ifelse(data_Id$protocol_type == "Traveling", 1, 0)
  cat("protocol_type:\n")
  print(table(data_Id$protocol_type))
  cat("locality_type:\n")
  print(table(data_Id$locality_type))
  cat("observation count binary format:\n")
  print(table(data_Id$binary_count))
  
  # Calculate the proportion of sample observations with zero number
  zero_proportion <- sum(data_Id$observation_count==0)/length(data_Id$observation_count)
  print(paste0("sample zero proportion = ", zero_proportion))
  
  #sample_size <- ceiling(mean(table(data_Id$year_number)))
  #sample_size
  #sample_Id_size <- ceiling(sample_size/length(table(data_Id$Id)))
  #sample_Id_size
  
  print('Calculate the observation bias correction index...')
  #CCI_output_folder <- "D:/aiv_project_table/glmm_performance"
  data_Id <- calculate_observer_cci(data_Id, birdname=bird_name, cci_output_folder=CCI_output_folder)
  
  colnames(data_Id)
  data_Id_only_effort <- data_Id[,c("observer_id", "Id", "year_number", "month_number", "observation_count", "locality_type", 
                                    "protocol_type", "duration_mins", "effort_km", "id_observation_times", "cci")]
  head(data_Id_only_effort)
  max(data_Id_only_effort$cci)
  
  # find the observer with the largest and smallest CCI
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
  
  
  print('Split training set and testing set...')
  #data_Id$sample_key <- paste0(data_Id$Id, "_", data_Id$observation_date)
  data_Id$sample_key <- 1:nrow(data_Id)
  unique_data_Id <- data_Id[,c('Id', 'year_number', 'month_number', "id_observation_times")]
  unique_data_Id <- unique_data_Id %>%
    distinct()
  print(paste0('number of pair (grid, time) : ', nrow(unique_data_Id)))
  
  set.seed(123)  
  n <- nrow(data_Id)
  val_indices <- sample(1:n, size = floor(0.1 * n))
  validation_data <- data_Id[val_indices, ]
  data_Id <- data_Id[-val_indices, ]
  validation_observation_count_data <- validation_data[,c('sample_key', 'observation_count')]
  
  write.csv(validation_observation_count_data,
            paste0(validtion_data_folder, bird_name, '.csv'),
            row.names = FALSE)
  
  
  print('start training the model and predicting results...')  
  result <- spatiotemporal_sampling(
    data_Id = data_Id,
    validation_data = validation_data,
    unique_data_Id = unique_data_Id,
    en_df = en_df,
    lc_df = lc_df,
    n_iter = 3,
    seed = 123,
    bird_name = bird_name,
    output_folder = abundance_model_output_folder
  )
  
}
