# 1. Load packages ----------
library(mgcv)
library(parallel)
library(pdp)
library(viridis)
library(fields)
library(rnaturalearth)
library(sf)
library(terra)
library(ebirdst)
library(ggplot2)
library(raster)
library(dplyr)
library(gbm)
library(ROSE)
library(xgboost)
library(caret)
library(gbm)
library(lubridate)
library(tidyr)
library(zoo)
library(forecast)
library(lme4)
library(MASS)
library(brms)
library(glmmTMB)
library(car)
library(DHARMa)
library(data.table)
library(pbapply)

# 2-1. Input path ----------
EU_shp_path <- "D:/aiv_project_2025/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp"
origin_lc_path <-"D:/aiv_project_2025/gee data/EU_2016_to_2022_land_cover_and_climate_data_containing_missing_values.csv"

# 2-2. Output path ----------
sampling_output_folder <- "D:/aiv_project_2025/google_land_cover_sampling_data/"
ML_imputation_output_folder <- "D:/aiv_project_2025/google_land_cover_sampling_data/ML_method_prediction_output/"
two_method_performance_output_folder <- "D:/aiv_project_2025/google_land_cover_sampling_data/two_method_performance/"
two_method_test_imputation_output_folder <- "D:/aiv_project_2025/google_land_cover_sampling_data/two_method_test_output/"


# 3. Read European shp files & land cover data ----------
EU.map <- st_read(EU_shp_path, quiet = TRUE)
EU_centroids <- st_centroid(EU.map)
EU_center_coordinate_df <- data.frame(
  Id = EU_centroids$Id,                                
  lat = st_coordinates(EU_centroids)[, "Y"],     # extract latitude
  lon = st_coordinates(EU_centroids)[, "X"]      # extract longitude
)


land_cover_data <- read.csv(origin_lc_path)
land_cover_data <- merge(land_cover_data, EU_center_coordinate_df, by="Id")
land_cover_data <- land_cover_data[order(land_cover_data$Id, land_cover_data$year_number, land_cover_data$month_number), ] 
grid_Id_table <- table(land_cover_data$Id)


# 4. Define sampling function -----------
run_sampling_interpolation <- function(seed_num, land_cover_data, grid_Id_table, sampling_out_folder) {
  land_cover_interpolation <- data.frame()
  
  #seed_num <- 123
  for (j in seq_along(grid_Id_table)) {
    # j <- 6
    set.seed(seed_num)
    Id_num <- names(grid_Id_table)[j]
    single_grid_data <- land_cover_data[land_cover_data$Id == Id_num, ]
    grass_MEAN <- single_grid_data$grass_MEAN
    non_na <- !is.na(grass_MEAN)                                                
        
    
    # apply run-length encoding to the logical vector 'non_na'
    # rle_result contains:
    # $values  : the value (TRUE/FALSE) of each consecutive segment
    # $lengths : the length of each segment (how many times it repeats consecutively)
    rle_result <- rle(non_na)                                                   
    segments <- which(rle_result$values & rle_result$lengths >= 3)              
    
    
    # For each selected segment index, compute its start and end positions
    continuous_segments <- lapply(segments, function(idx) {
      
      # Calculate the start position of the segment
      # If it is the first segment, it starts at position 1
      # Otherwise, start = sum of lengths of all previous segments + 1
      
      start <- if (idx == 1) 1 else sum(rle_result$lengths[1:(idx - 1)]) + 1
      end <- start + rle_result$lengths[idx] - 1
      list(start = start, end = end)
    })
    
    
    # For each continuous segment, sample positions based on segment length
    sample_positions <- lapply(continuous_segments, function(segment) {
      
      start <- segment$start
      end <- segment$end
      length <- end - start + 1
      
      # Sampling strategy based on segment length
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
        
        #enforce adjacency (make them consecutive if not already)
        if (sampled[2] - sampled[1] != 1) sampled[2] <- sampled[1] + 1
        other_candidates <- setdiff(candidates, sampled)
        sampled <- c(sampled, sample(other_candidates, 1))
        
      } else {
        sampled <- NA
      }
      sampled
    })
    
    
    # Combine sampled positions from all segments into a single data.frame
    # For each segment i:
    # segment_id: index of the segment
    # position: sampled positions within that segment
    sample_positions_df <- do.call(rbind, lapply(seq_along(sample_positions), function(i) {
      data.frame(segment_id = i, position = sample_positions[[i]])
    }))
    
    
    # If no sampled positions exist, skip this iteration
    if (is.null(sample_positions_df) || nrow(sample_positions_df) == 0) next
    
    
    #colnames(single_grid_data)[4:12]
    land_cover_types <- c("bare", "built", "crops", "flooded_vegetation", "grass",
                          "shrub_and_scrub", "snow_and_ice", "trees", "water")
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
    
    
    # Perform interpolation validation for each segment
    final_interpolated_results <- do.call(rbind, lapply(unique(sample_positions_df$segment_id), function(seg_id) {
      
      # extract sampled positions for this segment
      segment_positions <- sample_positions_df[sample_positions_df$segment_id == seg_id, "position"]
      start_idx <- continuous_segments[[seg_id]]$start
      end_idx <- continuous_segments[[seg_id]]$end
      
      result_rows <- data.frame(segment_id = seg_id, position = segment_positions)
      for (lc_type in land_cover_types) {
        segment_values <- land_cover_values[[lc_type]]
        current_segment_values <- segment_values[start_idx:end_idx]
        
        # Step 1: mask sampled positions (set them to NA)
        # Convert global index → local segment index
        current_segment_values[segment_positions - start_idx + 1] <- NA
        
        # Step 2: perform interpolation (linear interpolation)
        interpolated_segment <- na.approx(current_segment_values, na.rm = FALSE)
        
        # Step 3: extract interpolated values at sampled positions
        interpolated_values <- interpolated_segment[segment_positions - start_idx + 1]
        
        # Step 4: Get original true values at sampled positions
        original_values <- segment_values[segment_positions]
        
        # Step 5: Compute residuals (absolute error)
        residuals <- abs(interpolated_values - original_values)
        
        result_rows[[paste0(lc_type, "_original")]] <- original_values
        result_rows[[paste0(lc_type, "_interpolated")]] <- interpolated_values
        result_rows[[paste0(lc_type, "_residual")]] <- residuals
      }
      result_rows
    }))
    
    single_grid_result <- single_grid_data[, 1:3]
    index <- 3
    for (lc_type in land_cover_types) {
      index=index+1
      single_grid_result[[paste0(lc_type, "_original")]] <- single_grid_data[,index]
      single_grid_result[[paste0(lc_type, "_interpolated")]] <- NA
      single_grid_result[[paste0(lc_type, "_residual")]] <- NA
    }
    
   
    
    for (i in 1:nrow(final_interpolated_results)) {
      index <- final_interpolated_results$position[i]
      for (lc_type in land_cover_types) {
        single_grid_result[[paste0(lc_type, "_interpolated")]][index] <- final_interpolated_results[[paste0(lc_type, "_interpolated")]][i]
        single_grid_result[[paste0(lc_type, "_residual")]][index] <- final_interpolated_results[[paste0(lc_type, "_residual")]][i]
      }
    }
    
    land_cover_interpolation <- rbind(land_cover_interpolation, single_grid_result)
  }
  
  # output file
  out_path <- paste0(sampling_output_folder, "land_cover_seed", seed_num, ".csv")
  write.csv(land_cover_interpolation, out_path, row.names = FALSE)
  return(paste("Seed", seed_num, " finished"))
}



run_imputation <- function(seed_num, sampling_out_folder, land_cover_data) {
  #seed_num <- 123
  lc_path <- paste0(sampling_output_folder, "land_cover_seed", seed_num, ".csv")
  lc_imputation <- read.csv(lc_path)
  
  #ncol(land_cover_interpolation)
  imputation_cols_to_check <- 4:30
  lc_imputation_train_test_data <- lc_imputation[!apply(lc_imputation[, imputation_cols_to_check],
                                                        1, function(row) all(is.na(row))), ]
  
  lc_imputation_predict_data <- lc_imputation[apply(lc_imputation[, imputation_cols_to_check],
                                                    1, function(row) all(is.na(row))), ]
  
  lc_imputation_test_data <- lc_imputation_train_test_data[!is.na(lc_imputation_train_test_data$grass_interpolated), ]
  lc_imputation_train_data <- lc_imputation_train_test_data[is.na(lc_imputation_train_test_data$grass_interpolated), ]
  
  # take out "Id", "year_number" and "month_number" respectively to create train_data, test_data, predict_data
  lc_test_data <- lc_imputation_test_data[,c(1:3)]
  lc_train_data <- lc_imputation_train_data[,c(1:3)]
  lc_predict_data <- lc_imputation_predict_data[,c(1:3)]
  
  # merge the divided train_data, test_data, predict_data according to "Id", "year_number", "month_number" and merge the land_cover_data 
  lc_train_data <- merge(lc_train_data, land_cover_data, 
                                 by=c('Id','year_number','month_number'), all.x = FALSE )
  lc_test_data <- merge(lc_test_data, land_cover_data, 
                                by=c('Id','year_number','month_number'), all.x = FALSE )
  lc_predict_data <- merge(lc_predict_data, land_cover_data, 
                                   by=c('Id','year_number','month_number'), all.x = FALSE )
  
  lc_train_data$Id <- as.factor(lc_train_data$Id)
  lc_train_data$year_number <- as.factor(lc_train_data$year_number)
  lc_train_data$month_number <- as.factor(lc_train_data$month_number)
  
  lc_test_data$Id <- as.factor(lc_test_data$Id)
  lc_test_data$year_number <- as.factor(lc_test_data$year_number)
  lc_test_data$month_number <- as.factor(lc_test_data$month_number)
  
  lc_predict_data$Id <- as.factor(lc_predict_data$Id)
  lc_predict_data$year_number <- as.factor(lc_predict_data$year_number)
  lc_predict_data$month_number <- as.factor(lc_predict_data$month_number)
  
  exclude_columns <- c("Id", "year_number", "month_number", "grass_MEAN", 
                       "shrub_and_scrub_MEAN", "snow_and_ice_MEAN", "trees_MEAN", 
                       "water_MEAN", "bare_MEAN", "built_MEAN", "crops_MEAN", 
                       "flooded_vegetation_MEAN") 
  remaining_columns <- setdiff(names(lc_train_data), exclude_columns)
  train_features <- lc_train_data[, remaining_columns] 
  train_features <- as.matrix(train_features)
  
  test_features <- lc_test_data[, remaining_columns]
  test_features <- as.matrix(test_features)
  
  predict_features <- lc_predict_data[, remaining_columns]
  predict_features <- as.matrix(predict_features)
  
  land_cover_types <- c("bare", "built", "crops", "flooded_vegetation", "grass", 
                        "shrub_and_scrub", "snow_and_ice", "trees", "water")
  
  output_two_method_df <- lc_imputation_test_data[,1:3]
  output_ML_predict_df <- lc_predict_data[,1:3]
  output_validation_df <- {}
  
  for (i in 1:9){
    lc_type <- land_cover_types[i]
    set.seed(123)
    train_target <- log(lc_train_data[[paste0(lc_type, "_MEAN")]]/
                          (1-lc_train_data[[paste0(lc_type, "_MEAN")]]))
    train_target <- as.numeric(train_target)
    
    dtrain <- xgb.DMatrix(data = train_features, label = train_target)
    dtest <- xgb.DMatrix(data = test_features)
    dpred <- xgb.DMatrix(data = predict_features)
    
    xgb_model <- xgboost(data = dtrain,
                         objective = "reg:squarederror", 
                         max_depth = 15,
                         eta = 0.1,
                         nrounds = 300,
                         subsample = 1,           
                         colsample_bytree = 1,
                         nthread = 5)
    
    predictions_test <- 1 / (1 + exp(-predict(xgb_model, dtest)))
    na_predictions <- 1 / (1 + exp(-predict(xgb_model, dpred)))
    
    ML_rmse_test <- sqrt(mean((predictions_test - lc_test_data[[paste0(lc_type, "_MEAN")]])^2))
    ML_mae_test <- mean(abs(predictions_test - lc_test_data[[paste0(lc_type, "_MEAN")]]))
    linear_rmse_test <- sqrt(mean((lc_imputation_test_data[[paste0(lc_type, "_residual")]])^2))
    linear_mae_test <- mean(lc_imputation_test_data[[paste0(lc_type, "_residual")]]) 
    
    print(paste("Seed", seed_num, "-", lc_type, "ML_RMSE:", ML_rmse_test))
    print(paste("Seed", seed_num, "-", lc_type, "ML_mae:", ML_mae_test))
    
    ML_column_name <- paste0(lc_type, "_ML_interpolation")
    linear_column_name <- paste0(lc_type, "_linear_interpolation")
    original_column_name <- paste0(lc_type, "_original")
    predict_ML_column_name <- paste0(lc_type, "_MEAN")
    
    ML_comparison_df <- cbind(lc_test_data[,1:3],
                              setNames(data.frame(predictions_test), ML_column_name))
    linear_comparison_df <- cbind(lc_imputation_test_data[,1:3],
                                  setNames(data.frame(lc_imputation_test_data[[paste0(lc_type, "_original")]]), original_column_name),
                                  setNames(data.frame(lc_imputation_test_data[[paste0(lc_type, "_interpolated")]]), linear_column_name))
    two_method_df <- merge(linear_comparison_df, ML_comparison_df,
                           by=c('Id', 'year_number', 'month_number'), all.x=FALSE)                              
    two_method_df[[paste0(lc_type, "_linear_residual")]] <- abs(two_method_df[[paste0(lc_type, "_original")]] -
                                                                  two_method_df[[paste0(lc_type, "_linear_interpolation")]]) 
    two_method_df[[paste0(lc_type, "_ML_residual")]] <- abs(two_method_df[[paste0(lc_type, "_original")]] -
                                                              two_method_df[[paste0(lc_type, "_ML_interpolation")]])
    
    ML_predict_df <- cbind(lc_predict_data[,1:3],
                           setNames(data.frame(na_predictions), predict_ML_column_name))
    validation_df <- cbind(linear_rmse=linear_rmse_test, ML_rmse=ML_rmse_test,
                           linear_mae=linear_mae_test, ML_mae=ML_mae_test)
    
    output_two_method_df <- merge(output_two_method_df, two_method_df,
                                  by=c('Id', 'year_number', 'month_number'), all.x=FALSE)
    output_validation_df <- rbind(output_validation_df, validation_df)
    rownames(output_validation_df)[nrow(output_validation_df)] <- lc_type
    output_ML_predict_df <- merge(output_ML_predict_df, ML_predict_df,
                                  by=c('Id', 'year_number', 'month_number'), all.x=FALSE)
    
    #lc_type <- 'bare'
    write.csv(output_ML_predict_df,
              paste0(ML_imputation_output_folder, lc_type, "/seed", seed_num, ".csv"),
              row.names = F)
  }
  write.csv(output_validation_df,
            paste0(two_method_performance_output_folder, "/seed", seed_num, ".csv"))
  write.csv(output_two_method_df,
            paste0(two_method_test_imputation_output_folder, "/seed", seed_num, ".csv"),
            row.names = F)
  write.csv(output_ML_predict_df,
            paste0(ML_imputation_output_folder, "/seed", seed_num, ".csv"),
            row.names = F)
  
}


# 5-1. Run parallel land cover sampling ----------
#seed_list <- 123:(123 + 99)  
#n_cores <- detectCores() - 6
n_cores <- 2

cl <- makeCluster(n_cores)
clusterExport(cl, c("land_cover_data", "grid_Id_table", "run_sampling_interpolation", "sampling_out_folder"))
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  library(zoo)
})


results <- pbapply::pblapply(seed_list, function(s) {
  run_sampling_interpolation(s, land_cover_data, grid_Id_table, sampling_out_folder)
}, cl = cl)

stopCluster(cl)
cat(unlist(results), sep = "\n")


# 5-2. Run parallel land cover imputation ----------
n_cores <- 2
cl <- makeCluster(n_cores)
clusterExport(cl, c("run_imputation", "land_cover_data", 'sampling_out_folder'))

# 每個核心載入需要的套件
clusterEvalQ(cl, {
  library(xgboost)
  library(dplyr)
  library(data.table)
  library(zoo)
})

# 執行平行化＋進度條
results <- pbapply::pblapply(127:128, function(seed_num) {
  run_seed_process(seed_num)
}, cl = cl)

stopCluster(cl)





# 6. Aggregate multi-seed ML imputation results & compute mean value ----------
all_lc_df <- land_cover_data[,c("Id", "year_number", "month_number")]
for (lc_type in land_cover_types){
  
  single_lc_df <- land_cover_data[,c("Id", "year_number", "month_number", paste0(lc_type, '_MEAN'))]
  temp_lc_df <- single_lc_df[,c("Id", "year_number", "month_number")]
  for (seed_num in 123:222){
    #seed_num <- 123
    lc_ip_df <- read.csv(paste0(ML_imputation_output_folder, "/seed", seed_num, ".csv"))
    lc_ip_df <- lc_ip_df[,c("Id", "year_number", "month_number", paste0(lc_type, '_MEAN'))]
    colnames(lc_ip_df)[4] <- paste0('ML_', lc_type, '_MEAN_', seed_num)
    temp_lc_df <- merge(temp_lc_df, lc_ip_df, by=c("Id", "year_number", "month_number"), all=T)
    
  }
  
  result_df <- temp_lc_df %>%
    group_by(Id, year_number, month_number) %>%
    summarise(
      mean_value = mean(c_across(starts_with(paste0('ML_', lc_type, '_MEAN_'))), na.rm = TRUE),
      sd_value   = sd(c_across(starts_with(paste0('ML_', lc_type, '_MEAN_'))), na.rm = TRUE),
      .groups = "drop"
    )
  
  single_lc_df <- merge(single_lc_df, result_df[,1:4], by=c("Id", "year_number", "month_number"), all=T)
  single_lc_df[, 4][is.na(single_lc_df[, 4]) & !is.na(single_lc_df[, 5])] <- 
    single_lc_df[, 5][is.na(single_lc_df[, 4]) & !is.na(single_lc_df[, 5])]
  single_lc_df <- single_lc_df[,1:4]
  single_lc_df <- single_lc_df[!is.na(single_lc_df[, 4]), ]
  single_lc_df <- single_lc_df[order(single_lc_df$Id,
                                     single_lc_df$year_number,
                                     single_lc_df$month_number), ] 
  
  all_lc_df <- merge(all_lc_df, single_lc_df, by=c("Id", "year_number", "month_number"), all=T)
}

# create a clean copy and retain only rows with complete land cover information
all_lc_df_copy <- all_lc_df
all_lc_df_copy <- all_lc_df_copy[complete.cases(all_lc_df_copy[, 4:12]), ]
write.csv(all_lc_df_copy,
          "D:/aiv_project_2025/gee data/EU_2016_to_2022_land_cover_imputation_by_xgboost.csv",
          row.names = F)




