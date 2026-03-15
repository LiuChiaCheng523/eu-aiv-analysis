# 1. Load Package ----------
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
library(gridExtra)  
library(rnaturalearthdata)
library(ggrepel)
library(ggspatial)
library(RColorBrewer)
library(scales)
library(geosphere)
library(cowplot)

# 2-1. Input path ----------
EU_shp_path <- "D:/aiv_project_2025/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp"
chicken_density_path <- "D:/aiv_project_2025/livestock density 10km/chicken livestock density 10km.csv"
duck_density_path <- "D:/aiv_project_2025/livestock density 10km/duck livestock density_2015_10km.csv"
EU_aiv_2022_path <- "D:/aiv_project_2025/aiv fixed data/EU aiv fixed data 2022.csv"
EU_aiv_2021_path <- "D:/aiv_project_2025/aiv fixed data/EU aiv fixed data 2021.csv"
birdname_folder_path <- "D:/aiv_project_2025/abundance_estimation_results/abundance spatial method/MAD filter abundance"


# 2-2. Output path ----------
chicken_density_output_path <- "D:/aiv_project_plot/chicken_density.png"
duck_density_output_path <- "D:/aiv_project_plot/duck_density.png"
domestic_outbreak_output_path <- "D:/aiv_project_plot/Domestic_outbreak_map.png"
wild_outbreak_output_path <- "D:/aiv_project_plot/Wild outbreak map.png"
write_csv_date <- '20260315'
aiv_analysis_output_folder <- 'D:/aiv_project_table/aiv_analysis/'
weighted_abundance_output_folder <- 'D:/aiv_project_plot/weighted_abundance/'


# 3-1. Read European shp files ----------
EU.map = st_read(EU_shp_path)
EU_map <- as.data.frame(EU.map)
head(EU.map)
neighbors <- st_touches(EU.map)                                                 # Compute geographic grid neighbors


#calculate center point
EU_centroids <- st_centroid(EU.map)
EU_center_coordinate_df <- data.frame(
  Id = EU_centroids$Id,                                                         # extract geographic grid Ids
  latitude = st_coordinates(EU_centroids)[, "Y"],                               # extract latitude
  longitude = st_coordinates(EU_centroids)[, "X"]                               # extract longitude
)


# 3-2. Build a grid adjacency matrix ----------
matrix_data <- matrix(nrow = length(EU.map$Id), ncol = 9)
colnames(matrix_data) <- c('Id', paste0('neighbor', 1:8))
matrix_data[, 1] <- EU.map$Id


debug_neighbor_count <- 0
for (i in 1:length(EU.map$Id)) {
  neighbors_indices <- neighbors[[i]]                                           # get the index of adjacent grid
  
  if (length(neighbors_indices) > 0) {                                          # check if there are adjacent grids
    neighbors_id <- EU.map[neighbors_indices, ]                                 # extract neighboring grids based on index
    
    for (j in 1:length(neighbors_id$Id)) {
      matrix_data[i, j+1] <- neighbors_id$Id[j]                                 # fill the matrix with the IDs of adjacent grids
    }
  }
  
  debug_neighbor_count <- debug_neighbor_count + 1
}
neighbors_df<-as.data.frame(matrix_data)


# 3-3. Read FAO simulated livestock density ----------
chicken_livestock_density_df <- read.csv(chicken_density_path)
duck_livestock_density_df <- read.csv(duck_density_path)


# remove rows with missing Id and calculate total livestock density (VALUE) per grid
chicken_livestock_sum_density_df <- chicken_livestock_density_df %>%
  filter(!is.na(Id)) %>%        
  group_by(Id) %>%             
  summarise(sum_VALUE = sum(VALUE, na.rm = TRUE)) 
duck_livestock_sum_density_df <- duck_livestock_density_df %>%
  filter(!is.na(Id)) %>%
  group_by(Id) %>%
  summarise(sum_VALUE = sum(VALUE, na.rm = TRUE))
print(paste0('chicken livestock sum density maximun : ', max(chicken_livestock_sum_density_df$sum_VALUE)))
print(paste0('chicken livestock sum density minimun : ', min(chicken_livestock_sum_density_df$sum_VALUE)))
print(paste0('duck livestock sum density maximun : ', max(duck_livestock_sum_density_df$sum_VALUE)))
print(paste0('duck livestock sum density minimun : ', min(duck_livestock_sum_density_df$sum_VALUE)))


# merge livestock density with grid map and replace missing values with 0
chicken_livestock_dat <- merge(EU.map, chicken_livestock_sum_density_df, by='Id', all=TRUE)
chicken_livestock_dat$sum_VALUE[is.na(chicken_livestock_dat$sum_VALUE)] <- 0
duck_livestock_dat <- merge(EU.map, duck_livestock_sum_density_df, by='Id', all=TRUE)
duck_livestock_dat$sum_VALUE[is.na(duck_livestock_dat$sum_VALUE)] <- 0


# 3-4. Plot spatial distribution of chicken density by quartile ----------
# a. chicken livestock
chicken_livestock_dat <- chicken_livestock_dat %>%
  mutate(sum_VALUE_quartile = ntile(sum_VALUE, 4))

quartile_breaks <- quantile(chicken_livestock_dat$sum_VALUE, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
quartile_labels <- c(
  paste0(formatC(quartile_breaks[1], format = "d", big.mark = ","), " - ", formatC(quartile_breaks[2], format = "d", big.mark = ",")),
  paste0(formatC(quartile_breaks[2]+1, format = "d", big.mark = ","), " - ", formatC(quartile_breaks[3], format = "d", big.mark = ",")),
  paste0(formatC(quartile_breaks[3]+1, format = "d", big.mark = ","), " - ", formatC(quartile_breaks[4], format = "d", big.mark = ",")),
  paste0(formatC(quartile_breaks[4]+1, format = "d", big.mark = ","), " - ", formatC(quartile_breaks[5], format = "d", big.mark = ","))
)

chicken_livestock_map <- ggplot() +
  geom_sf(data = chicken_livestock_dat, 
          aes(fill = as.factor(sum_VALUE_quartile)), 
          color = "gray30", size = 0.1) +  # 加上網格邊界
  scale_fill_brewer(palette = "Reds", name = "range", labels = quartile_labels) +
  coord_sf(xlim = c(-25, 70), ylim = c(33, 79), expand = FALSE) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  theme_dark(base_size = 12) +
  ggtitle("chicken livestock density")
chicken_livestock_map


# b. duck livestock
duck_livestock_dat <- duck_livestock_dat %>%
  mutate(sum_VALUE_quartile = ntile(sum_VALUE, 4))

quartile_breaks <- quantile(duck_livestock_dat$sum_VALUE, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
quartile_labels <- c(
  paste0(formatC(quartile_breaks[1], format = "d", big.mark = ","), " - ", formatC(quartile_breaks[2], format = "d", big.mark = ",")),
  paste0(formatC(quartile_breaks[2]+1, format = "d", big.mark = ","), " - ", formatC(quartile_breaks[3], format = "d", big.mark = ",")),
  paste0(formatC(quartile_breaks[3]+1, format = "d", big.mark = ","), " - ", formatC(quartile_breaks[4], format = "d", big.mark = ",")),
  paste0(formatC(quartile_breaks[4]+1, format = "d", big.mark = ","), " - ", formatC(quartile_breaks[5], format = "d", big.mark = ","))
)

duck_livestock_map <- ggplot() +
  geom_sf(data = duck_livestock_dat, 
          aes(fill = as.factor(sum_VALUE_quartile)), 
          color = "gray30", size = 0.1) +  # 加上網格邊界
  scale_fill_brewer(palette = "Reds", name = "range", labels = quartile_labels) +
  coord_sf(xlim = c(-25, 70), ylim = c(33, 79), expand = FALSE) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  theme_dark(base_size = 12) +
  ggtitle("duck livestock density")
duck_livestock_map

ggsave(chicken_density_output_path,
       plot = chicken_livestock_map, width = 12, height = 6, dpi = 600, device = "png")
ggsave(duck_density_output_path,
       plot = duck_livestock_map, width = 12, height = 6, dpi = 600, device = "png")


# chicken and duck density distribution displayed side by side
plot_grid(chicken_livestock_map, duck_livestock_map,
          labels = c("A", "B"), 
          ncol = 2)  


# 3-5. Read FAO avian influenza outbreak case data ----------
aiv_2022_df <- read.csv(EU_aiv_2022_path)
aiv_2021_df <- read.csv(EU_aiv_2021_path)
EU_aiv_2022_df <- aiv_2022_df[!is.na(aiv_2022_df$Id),]
EU_aiv_2021_df <- aiv_2021_df[!is.na(aiv_2021_df$Id),]
EU_aiv_2022_df$Species <- sub(",.*", "", EU_aiv_2022_df$Species)
EU_aiv_2021_df$Species <- sub(",.*", "", EU_aiv_2021_df$Species)


# select outbreak case type (Domestic, Wild)
outbreak_type <- 'Wild'     
EU_aiv_2022_df <- EU_aiv_2022_df[which(EU_aiv_2022_df$Species==outbreak_type),]     
EU_aiv_2021_df <- EU_aiv_2021_df[which(EU_aiv_2021_df$Species==outbreak_type),]
EU_aiv_2022_sf <- st_as_sf(EU_aiv_2022_df, coords = c("longitude", "latitude"), crs = st_crs(EU.map))
EU_aiv_2021_sf <- st_as_sf(EU_aiv_2021_df, coords = c("longitude", "latitude"), crs = st_crs(EU.map))
EU_aiv_2021to2022_sf <- rbind(EU_aiv_2021_sf, EU_aiv_2022_sf)


# 3-6. Plot the spread of the 2021-2022 avian influenza outbreak ----------
outbreak_map <- ggplot() +
  geom_sf(data = EU.map, fill = "white", color = "black", size = 0.1) +
  geom_sf(data = EU_aiv_2021to2022_sf, 
          color = "red", size = 0.7, alpha = 0.9) +
  coord_sf(xlim = c(-25, 70), ylim = c(33, 80), expand = FALSE) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  labs(
    title = paste0("AIV Outbreak (", outbreak_type, " case) on European grid map (2021–2022)"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_dark(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  annotate("text", x = 68, y = 34, 
           label = "Source: FAO Empress-i ", size = 3.5, hjust = 1, color = "black")
outbreak_map

ggsave(domestic_outbreak_output_path,
       plot = outbreak_map, width = 8, height = 6, dpi = 300)
ggsave(wild_outbreak_output_path,
       plot = outbreak_map, width = 8, height = 6, dpi = 600)


# 3-7. Combine AIV outbreak records (2021–2022), extract year-month information ----------
EU_aiv_df <- rbind(EU_aiv_2022_df, EU_aiv_2021_df)
head(EU_aiv_df$observation.date)

EU_aiv_df$observation.date <- ymd_hms(EU_aiv_df$observation.date)
EU_aiv_df <-EU_aiv_df[!is.na(EU_aiv_df$observation.date),]
EU_aiv_df$observation_year <- year(EU_aiv_df$observation.date)
EU_aiv_df$observation_month <- month(EU_aiv_df$observation.date)
head(EU_aiv_df[, c("observation.date", "observation_year", "observation_month")])
table(EU_aiv_df$observation_year)
table(EU_aiv_df$observation_month)

EU_aiv_df$observation_year_month <- paste(EU_aiv_df$observation_year, 
                                          sprintf("%02d", EU_aiv_df$observation_month), sep = "-")


# 3-8. plot monthly outbreak counts ----------
outbreak_time_series_data <- table(EU_aiv_df$observation_year_month)
outbreak_time_series_df <- as.data.frame(outbreak_time_series_data)
colnames(outbreak_time_series_df) <- c("observation_year_month", "outbreak_count")

ggplot(outbreak_time_series_df, aes(x = observation_year_month, y = outbreak_count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = paste0('Monthly Distribution of Avian Influenza Events(', outbreak_type, ') (2021-2022)'),
       x = "Year-Month",
       y = "Number of Events") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# 3-9. Statistics of the total number of avian influenza outbreaks in each geographical grid from 2021 to 2022 ----------
Id_outbreak_freq <- as.data.frame(table(EU_aiv_df$Id))
colnames(Id_outbreak_freq) <- c("Id", "freq21to22")
nrow(Id_outbreak_freq)     

# filter avian influenza outbreak records in different periods 
# stage1 (2021/05-2021/10), stage2 (2021/11-2022/04) and stage3 (2022/05-2022/10)
Id_outbreak_freq$stage1_freq <- 0
Id_outbreak_freq$stage2_freq <- 0
Id_outbreak_freq$stage3_freq <- 0

# pre-establish the mapping from Id to row_index to improve efficiency
# Find the row index in Id_outbreak_freq that corresponds to each Id in EU_aiv_df
id_to_row <- match(EU_aiv_df$Id, Id_outbreak_freq$Id) 

# classify AIV outbreaks into stage1, stage2, and stage3 according to observation year and month,
# and accumulate outbreak frequencies for each grid
if(outbreak_type=='Domestic'){
  for (i in 1:nrow(EU_aiv_df)){
    row_index <- id_to_row[i]
    if (EU_aiv_df$observation_year[i] == 2021){
      if (EU_aiv_df$observation_month[i] >= 5 & EU_aiv_df$observation_month[i] <= 10) {
        Id_outbreak_freq$stage1_freq[row_index] <- Id_outbreak_freq$stage1_freq[row_index] + 1
      }else if(EU_aiv_df$observation_month[i] >= 11){
        Id_outbreak_freq$stage2_freq[row_index] <- Id_outbreak_freq$stage2_freq[row_index] + 1
      }
    }else if (EU_aiv_df$observation_year[i] == 2022){
      if (EU_aiv_df$observation_month[i] <= 4) {
        Id_outbreak_freq$stage2_freq[row_index] <- Id_outbreak_freq$stage2_freq[row_index] + 1
      }else if(5 <= EU_aiv_df$observation_month[i] & EU_aiv_df$observation_month[i] <= 10){
        Id_outbreak_freq$stage3_freq[row_index] <- Id_outbreak_freq$stage3_freq[row_index] + 1
      }
    }
  }
  print('stages are seperated by outbreak of Domestic')
}else{
  for (i in 1:nrow(EU_aiv_df)){
    row_index <- id_to_row[i]
    if (EU_aiv_df$observation_year[i] == 2021){
      if (EU_aiv_df$observation_month[i] >= 5 & EU_aiv_df$observation_month[i] <= 10) {
        Id_outbreak_freq$stage1_freq[row_index] <- Id_outbreak_freq$stage1_freq[row_index] + 1
      }else if(EU_aiv_df$observation_month[i] >= 11){
        Id_outbreak_freq$stage2_freq[row_index] <- Id_outbreak_freq$stage2_freq[row_index] + 1
      }
    }else if (EU_aiv_df$observation_year[i] == 2022){
      if (EU_aiv_df$observation_month[i] <= 4) {
        Id_outbreak_freq$stage2_freq[row_index] <- Id_outbreak_freq$stage2_freq[row_index] + 1
      }else if(5 <= EU_aiv_df$observation_month[i] & EU_aiv_df$observation_month[i] <= 10){
        Id_outbreak_freq$stage3_freq[row_index] <- Id_outbreak_freq$stage3_freq[row_index] + 1
      }
    }
  }
  print('stages are seperated by outbreak of Wild')
}

nrow(Id_outbreak_freq)
nrow(Id_outbreak_freq[which(Id_outbreak_freq$stage1_freq > 0),])
nrow(Id_outbreak_freq[which(Id_outbreak_freq$stage2_freq > 0),])
nrow(Id_outbreak_freq[which(Id_outbreak_freq$stage3_freq > 0),])
nrow(Id_outbreak_freq[which(Id_outbreak_freq$stage1_freq > 0 | Id_outbreak_freq$stage2_freq > 0),])
nrow(Id_outbreak_freq[which(Id_outbreak_freq$stage1_freq > 0 | Id_outbreak_freq$stage3_freq > 0),])
nrow(Id_outbreak_freq[which(Id_outbreak_freq$stage2_freq > 0 | Id_outbreak_freq$stage3_freq > 0),])

# 3-10. Calculate the monthly weight of stage1 ~ stage3 ----------
sum_stage1_freq <- sum(Id_outbreak_freq$stage1_freq)
sum_stage2_freq <- sum(Id_outbreak_freq$stage2_freq)
sum_stage3_freq <- sum(Id_outbreak_freq$stage3_freq)

aiv_year_month_df <- EU_aiv_df[,c('Event.ID', 'observation.date', 'observation_year', 'observation_month')]

stage1_start <- as.Date("2021-05-01")
stage1_end <- as.Date("2021-10-31")
stage2_start <- as.Date("2021-11-01")
stage2_end <- as.Date("2022-04-30")
stage3_start <- as.Date("2022-05-01")
stage3_end <- as.Date("2022-10-31")

aiv_year_month_df$observation.date <- as.Date(aiv_year_month_df$observation.date)

stage1_aiv_month <- aiv_year_month_df %>%
  filter(observation.date >= stage1_start & observation.date <= stage1_end)
stage2_aiv_month <- aiv_year_month_df %>%
  filter(observation.date >= stage2_start & observation.date <= stage2_end)
stage3_aiv_month <- aiv_year_month_df %>%
  filter(observation.date >= stage3_start & observation.date <= stage3_end)

stage1_aiv_month_weighted <- stage1_aiv_month %>%
  group_by(observation_year, observation_month) %>%
  summarise(sample_count = n(), .groups = "drop")
stage2_aiv_month_weighted <- stage2_aiv_month %>%
  group_by(observation_year, observation_month) %>%
  summarise(sample_count = n(), .groups = "drop")
stage3_aiv_month_weighted <- stage3_aiv_month %>%
  group_by(observation_year, observation_month) %>%
  summarise(sample_count = n(), .groups = "drop")

stage1_aiv_month_weighted$month_weighted <- stage1_aiv_month_weighted$sample_count/sum_stage1_freq
stage2_aiv_month_weighted$month_weighted <- stage2_aiv_month_weighted$sample_count/sum_stage2_freq
stage3_aiv_month_weighted$month_weighted <- stage3_aiv_month_weighted$sample_count/sum_stage3_freq
print(stage1_aiv_month_weighted)
print(stage2_aiv_month_weighted)
print(stage3_aiv_month_weighted)

# equal month weight
#stage1_aiv_month_weighted$month_weighted <- 1/6
#stage2_aiv_month_weighted$month_weighted <- 1/6
#stage3_aiv_month_weighted$month_weighted <- 1/6

# 3-11. Calculate the number of neighboring grids with AIV outbreaks for each stage ----------
Id_df <- data.frame(Id=EU_map$Id)
Id_outbreak_freq <- merge(Id_df , Id_outbreak_freq, by='Id', all=TRUE)
Id_outbreak_freq$Id <- as.numeric(Id_outbreak_freq$Id)

Id_outbreak_freq$stage1_outbreak_nei <- 0
Id_outbreak_freq$stage2_outbreak_nei <- 0
Id_outbreak_freq$stage3_outbreak_nei <- 0
Id_outbreak_freq$stage1_freq[is.na(Id_outbreak_freq$stage1_freq)] <-0
Id_outbreak_freq$stage2_freq[is.na(Id_outbreak_freq$stage2_freq)] <-0
Id_outbreak_freq$stage3_freq[is.na(Id_outbreak_freq$stage3_freq)] <-0

for (i in 1:nrow(Id_outbreak_freq)){
  id_number <- as.numeric(Id_outbreak_freq$Id[i])
  print(id_number)
  id_neighbors <- neighbors_df[neighbors_df$Id == id_number, 2:9]
  id_neighbors
  id_neighbors <- id_neighbors[!is.na(id_neighbors)]
  id_neighbors
  if(length(id_neighbors) > 0){
    for (j in 1:length(id_neighbors)){
      nei_id_num <- id_neighbors[j]
      nei_id_num
      if(nei_id_num %in% Id_outbreak_freq$Id){
        if(Id_outbreak_freq$stage1_freq[which(Id_outbreak_freq$Id==nei_id_num)] > 0){
          Id_outbreak_freq$stage1_outbreak_nei[i] <- Id_outbreak_freq$stage1_outbreak_nei[i] + 1
        }
        if(Id_outbreak_freq$stage2_freq[which(Id_outbreak_freq$Id==nei_id_num)] > 0){
          Id_outbreak_freq$stage2_outbreak_nei[i] <- Id_outbreak_freq$stage2_outbreak_nei[i] + 1
        }
        if(Id_outbreak_freq$stage3_freq[which(Id_outbreak_freq$Id==nei_id_num)] > 0){
          Id_outbreak_freq$stage3_outbreak_nei[i] <- Id_outbreak_freq$stage3_outbreak_nei[i] + 1
        }
      }
    }
  }
}

# 3-12. Merge outbreak and livestock density data ----------
Id_df <- data.frame(Id=EU_map$Id)
chicken_livestock_df <- data.frame(chicken_livestock_dat)
duck_livestock_df <- data.frame(duck_livestock_dat)
chicken_livestock_df <- chicken_livestock_df[,c('Id', 'sum_VALUE')]   
duck_livestock_df <- duck_livestock_df[,c('Id', 'sum_VALUE')]   
names(chicken_livestock_df)[2] <- 'chicken_sum_VALUE'
names(duck_livestock_df)[2] <- 'duck_sum_VALUE'

all_df <- merge(Id_df, chicken_livestock_df, by='Id', all=TRUE)
all_df <- merge(all_df, duck_livestock_df, by='Id', all=TRUE)
all_df <- merge(all_df, Id_outbreak_freq, by='Id', all=TRUE)
all_df$stage1_outbreak_binary <- ifelse(all_df$stage1_freq >0, 1, 0)
all_df$stage2_outbreak_binary <- ifelse(all_df$stage2_freq >0, 1, 0)
all_df$stage3_outbreak_binary <- ifelse(all_df$stage3_freq >0, 1, 0)
head(all_df)


# 4. Define function ----------

get_weighted_abundance <- function(abundance){
  weighted_abundance <- abundance
  weighted_abundance$unweighted_abundance <- weighted_abundance$abundance_filtered_mean
  
  # create date field
  weighted_abundance <- weighted_abundance %>%
    mutate(date = as.Date(paste(year_number, month_number, "01", sep = "-")))
  
  # weighted aggregated data
  weighted_abundance_result <- weighted_abundance %>%
    mutate(stage = case_when(
      date >= stage1_start & date <= stage1_end ~ "stage1",
      date >= stage2_start & date <= stage2_end ~ "stage2",
      date >= stage3_start & date <= stage3_end ~ "stage3",
      TRUE ~ NA_character_  
    )) %>%
    filter(!is.na(stage)) %>%  # remove data that does not belong to stage1 ~ stage3
    mutate(
      observation_year = as.numeric(format(date, "%Y")),
      observation_month = as.numeric(format(date, "%m"))
    ) %>%
    left_join(
      bind_rows(
        stage1_aiv_month_weighted %>% mutate(stage = "stage1"),
        stage2_aiv_month_weighted %>% mutate(stage = "stage2"),
        stage3_aiv_month_weighted %>% mutate(stage = "stage3")
      ),
      by = c("observation_year", "observation_month", "stage")
    ) %>%
    mutate(weighted_abundance = unweighted_abundance * month_weighted) %>%
    group_by(Id, stage) %>%
    summarise(
      sum_weighted_abundance = sum(weighted_abundance, na.rm = TRUE),  # 按加權值聚合
      .groups = "drop"
    ) %>%
    complete(Id, stage = c("stage1", "stage2", "stage3"), fill = list(sum_weighted_abundance = 0))
  
  return(weighted_abundance_result)
}


merge_all_stages_data <- function(all_df, weighted_abundance_result, EU_center_coordinate_df){
  stage1_result <- weighted_abundance_result[,c(1,3)][which(weighted_abundance_result$stage=='stage1'),]
  stage2_result <- weighted_abundance_result[,c(1,3)][which(weighted_abundance_result$stage=='stage2'),]
  stage3_result <- weighted_abundance_result[,c(1,3)][which(weighted_abundance_result$stage=='stage3'),]
  
  colnames(stage1_result)[2] <- 'weighted_abundance_stage1'   
  colnames(stage2_result)[2] <- 'weighted_abundance_stage2'   
  colnames(stage3_result)[2] <- 'weighted_abundance_stage3'   
  
  # merge data of all stages
  all_info_df <- merge(all_df, stage1_result, by='Id', all=TRUE)
  all_info_df <- merge(all_info_df, stage2_result, by='Id', all=TRUE)
  all_info_df <- merge(all_info_df, stage3_result, by='Id', all=TRUE)
  
  all_info_df <- merge(all_info_df, EU_center_coordinate_df, by='Id', all=TRUE)
  all_info_df <- all_info_df[!is.na(all_info_df$weighted_abundance_stage1),]   #weighted_expected_abundance_stage1
  
  all_info_df$stage1_freq[is.na(all_info_df$stage1_freq)] <- 0
  all_info_df$stage2_freq[is.na(all_info_df$stage2_freq)] <- 0
  all_info_df$stage3_freq[is.na(all_info_df$stage3_freq)] <- 0
  all_info_df$standard_chicken_density <- scale(all_info_df$chicken_sum_VALUE)
  all_info_df$standard_duck_density <- scale(all_info_df$duck_sum_VALUE)
  all_info_df$stage1_outbreak_nei[is.na(all_info_df$stage1_outbreak_nei)] <- 0
  all_info_df$stage2_outbreak_nei[is.na(all_info_df$stage2_outbreak_nei)] <- 0
  all_info_df$stage3_outbreak_nei[is.na(all_info_df$stage3_outbreak_nei)] <- 0
  
  return(all_info_df)
}


rbind_all_stages_data <- function(all_df, weighted_abundance_result, EU_center_coordinate_df){
  all_stage1_df <- all_df[,c('Id', 'chicken_sum_VALUE', 'duck_sum_VALUE', 'stage1_freq', 'stage1_outbreak_binary', 'stage1_outbreak_nei')]
  all_stage2_df <- all_df[,c('Id', 'chicken_sum_VALUE', 'duck_sum_VALUE', 'stage2_freq', 'stage2_outbreak_binary', 'stage2_outbreak_nei')]
  all_stage3_df <- all_df[,c('Id', 'chicken_sum_VALUE', 'duck_sum_VALUE', 'stage3_freq', 'stage3_outbreak_binary', 'stage3_outbreak_nei')]
  colnames(all_stage1_df)[c(4,5,6)] <- c('outbreak_freq', 'outbreak_binary', 'outbreak_nei')
  colnames(all_stage2_df)[c(4,5,6)] <- c('outbreak_freq', 'outbreak_binary', 'outbreak_nei')
  colnames(all_stage3_df)[c(4,5,6)] <- c('outbreak_freq', 'outbreak_binary', 'outbreak_nei')
  
  # add stage category variables
  all_stage1_df$stage <- 'stage1'
  all_stage2_df$stage <- 'stage2'
  all_stage3_df$stage <- 'stage3'
  
  # each stage combines outbreak (response variable), weighted abundance (explanatory variable))
  all_stage1_df <- merge(all_stage1_df, weighted_abundance_result, by=c('Id', 'stage'), all=FALSE)
  all_stage2_df <- merge(all_stage2_df, weighted_abundance_result, by=c('Id', 'stage'), all=FALSE)
  all_stage3_df <- merge(all_stage3_df, weighted_abundance_result, by=c('Id', 'stage'), all=FALSE)
  
  all_stage1_df <- merge(all_stage1_df, EU_center_coordinate_df, by='Id', all=FALSE)
  all_stage2_df <- merge(all_stage2_df, EU_center_coordinate_df, by='Id', all=FALSE)
  all_stage3_df <- merge(all_stage3_df, EU_center_coordinate_df, by='Id', all=FALSE)
  
  # merge data from stage1–stage3
  # replace missing outbreak values, define stage as an ordered factor
  # standardize chicken and duck livestock densities
  all_stage_df <- rbind(all_stage1_df, all_stage2_df, all_stage3_df)
  all_stage_df$outbreak_freq[is.na(all_stage_df$outbreak_freq)] <- 0
  all_stage_df$outbreak_binary[is.na(all_stage_df$outbreak_binary)] <- 0
  all_stage_df$stage <- factor(all_stage_df$stage, levels = c("stage1", "stage2", "stage3"))
  all_stage_df$standard_chicken_density <- scale(all_stage_df$chicken_sum_VALUE)
  all_stage_df$standard_duck_density <- scale(all_stage_df$duck_sum_VALUE)
  
  return(all_stage_df)
}


get_single_stage_gam_model <- function(birdname, stage_merge_df) {
  
  # Automatically fit a GAM model and iteratively increase the spatial smooth basis (k)
  # until mgcv k.check indicates an adequate basis dimension
  auto_gam_fit <- function(stage_var, data, kval_init = 30, max_k = 150, step_k = 20) {
    kval <- kval_init
    acceptable <- FALSE
    gam_fit <- NULL
    
    while (!acceptable && kval <= max_k) {
      
      formula_expr <- as.formula(
        paste0(stage_var, "_freq ~ weighted_abundance_", stage_var,
               " + standard_chicken_density + standard_duck_density + ",
               "s(latitude, longitude, k = ", kval, ", bs = 'tp')")
      )
      
      gam_fit <- tryCatch(
        {
          gam(
            formula = formula_expr,
            family = poisson(link = "log"),
            method = "REML",
            data = data
          )
        },
        warning = function(w) {
          message(paste0(stage_var, " model warning for bird ", birdname, " at k = ", kval, ":\n", conditionMessage(w)))
          return(NULL)
        },
        error = function(e) {
          message(paste0(stage_var, " model failed for bird ", birdname, " at k = ", kval, ":\n", conditionMessage(e)))
          return(NULL)
        }
      )
      
      if (is.null(gam_fit)) {
        return(NULL)
      }
      
      # Use mgcv:::k.check() to determine whether the basis dimension is sufficient
      kcheck <- tryCatch(
        {
          mgcv:::k.check(gam_fit)
        },
        error = function(e) {
          message(paste0("k.check failed for bird ", birdname, " (", stage_var, ") at k = ", kval, ":\n", conditionMessage(e)))
          return(NULL)
        }
      )
      
      if (is.null(kcheck)) {
        return(NULL)
      }
      
      pvals <- kcheck[, "p-value"]
      
      if (all(pvals >= 0.05)) {
        acceptable <- TRUE
      } else {
        kval <- kval + step_k
      }
    }
    
    return(gam_fit)
  }
  
  # perform model fitting for each stage
  gam_stage1 <- auto_gam_fit("stage1", stage_merge_df)
  gam_stage2 <- auto_gam_fit("stage2", stage_merge_df)
  gam_stage3 <- auto_gam_fit("stage3", stage_merge_df)
  
  return(list(
    gam_stage1 = gam_stage1,
    gam_stage2 = gam_stage2,
    gam_stage3 = gam_stage3
  ))
}


get_all_stage_model <- function(birdname, stage_rbind_df, kval_init = 50, max_k = 150, step_k = 20){
  
  kval <- kval_init
  acceptable <- FALSE
  gam_stage4 <- NULL
  
  while (!acceptable && kval <= max_k) {
    
    gam_stage4 <- tryCatch(
      {
        gam(
          outbreak_freq ~ stage*sum_weighted_abundance + standard_chicken_density + standard_duck_density +
            s(latitude, longitude, k = kval, bs = "tp") , 
          family = poisson(link = "log"),  
          method = "REML",
          data = stage_rbind_df
        )
      },
      warning = function(w) {
        message(paste("stage all satge model issued a warning for bird:", birdname, "\n", conditionMessage(w)))
        return(NULL)
      },
      error = function(e) {
        message(paste("stage all satge model failed for bird:", birdname, "\n", conditionMessage(e)))
        return(NULL)
      }
    )
    
    if (is.null(gam_stage4)) {
      return(NULL)
    }
    
    kcheck <- tryCatch(
      {
        mgcv:::k.check(gam_stage4)
      },
      error = function(e) {
        message(paste0("k.check failed for bird ", birdname, " (all stages at k = ", kval, ":\n", conditionMessage(e)))
        return(NULL)
      }
    )
    
    if (is.null(kcheck)) {
      return(NULL)
    }
    
    pvals <- kcheck[, "p-value"]
    
    if (all(pvals >= 0.05)) {
      acceptable <- TRUE
    } else {
      kval <- kval + step_k
    }
  }
  
  return(gam_stage4)
}


get_all_stages_stage_table <- function(birdname, gam_stage4, num_grid){
  p_coeff <- p_pv <- rep(NA, 8)
  
  if (!is.null(gam_stage4)) {
    p_coeff <- summary(gam_stage4)$p.coeff
    p_pv <- summary(gam_stage4)$p.pv
  }
  
  all_stages_stage_table <- data.frame(
    bird_name = birdname,
    coef_stagestage1 = p_coeff[2],
    p_stagestage1 = p_pv[2],
    coef_stagestage3 = p_coeff[3],
    p_stagestage3 = p_pv[3],
    num_grid=num_grid
  )
  return(all_stages_stage_table)
}


get_all_stages_density_table <- function(birdname, gam_stage4, num_grid){
  p_coeff <- p_pv <- rep(NA, 8)
  
  if (!is.null(gam_stage4)) {
    p_coeff <- summary(gam_stage4)$p.coeff
    p_pv <- summary(gam_stage4)$p.pv
  }
  
  all_stages_stage_table <- data.frame(
    bird_name = birdname,
    coef_chicken_density = p_coeff[5],
    p_chicken_density = p_pv[5],
    coef_duck_density = p_coeff[6],
    p_duck_density = p_pv[6],
    num_grid=num_grid
  )
  return(all_stages_stage_table)
}


get_all_stages_stageabundance_table <- function(birdname, gam_stage4, num_grid){
  p_coeff <- p_pv <- rep(NA, 8)
  
  if (!is.null(gam_stage4)) {
    p_coeff <- summary(gam_stage4)$p.coeff
    p_pv <- summary(gam_stage4)$p.pv
  }
  
  all_stages_stageabundance_table <- data.frame(
    bird_name = birdname,
    coef_abundance =  p_coeff[4],
    p_abundance = p_pv[4],
    coef_stage1abundance =  p_coeff[7],
    p_stage1abundance = p_pv[7],
    coef_stage3abundance =  p_coeff[8],
    p_stage3abundance = p_pv[8],
    num_grid=num_grid
  )
  return(all_stages_stageabundance_table)
}


get_single_stage_abundance_table <- function(birdname, single_stage_gam_model, num_grid){
  gam_stage1 <- single_stage_gam_model$gam_stage1
  gam_stage2 <- single_stage_gam_model$gam_stage2
  gam_stage3 <- single_stage_gam_model$gam_stage3
  
  # initialize variables
  p_coeff_stage1 <- p_pv_stage1 <- rep(NA, 4)
  p_coeff_stage2 <- p_pv_stage2 <- rep(NA, 4)
  p_coeff_stage3 <- p_pv_stage3 <- rep(NA, 4)
  
  # if gam_stage is neither NULL nor NA, execute summary(gam_stage)
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
  
  
  single_stage_abundance_table <- data.frame(
    bird_name = birdname,
    coef_abundance_stage1 = ifelse(length(p_coeff_stage1) >= 2, p_coeff_stage1[2], NA),
    p_abundance_stage1 = ifelse(length(p_pv_stage1) >= 2, p_pv_stage1[2], NA),
    coef_abundance_stage2 = ifelse(length(p_coeff_stage2) >= 2, p_coeff_stage2[2], NA),
    p_abundance_stage2 = ifelse(length(p_pv_stage2) >= 2, p_pv_stage2[2], NA),
    coef_abundance_stage3 = ifelse(length(p_coeff_stage3) >= 2, p_coeff_stage3[2], NA),
    p_abundance_stage3 = ifelse(length(p_pv_stage3) >= 2, p_pv_stage3[2], NA),
    num_grid = num_grid
  )
  return(single_stage_abundance_table)
}


get_single_stage_density_table <- function(birdname, single_stage_gam_model, num_grid){
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
  
  
  single_stage_density_table <- data.frame(
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
  return(single_stage_density_table)
}


save_outbreak_and_abundance_image <- function(EU.map, stage_rbind_df, stage, birdname, output_folder){
  
  # make sure the main folder exists
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # create subfolder (named with birdname)
  bird_folder <- file.path(output_folder, birdname)
  if (!dir.exists(bird_folder)) {
    dir.create(bird_folder)
  }
  
  # filter stage_rbind_df and merge with EU.map
  stage_df <- stage_rbind_df[which(stage_rbind_df$stage==stage),]
  stage_df <- stage_df[,c('Id','outbreak_freq','outbreak_binary','sum_weighted_abundance','stage')]
  outbreak_abundance_dat <- merge(EU.map, stage_df, by="Id", all=TRUE)
  
  # calculate quartiles of outcome and weighted abundance
  outbreak_quartiles <- quantile(outbreak_abundance_dat$outbreak_freq, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  abundance_quartiles <- quantile(outbreak_abundance_dat$sum_weighted_abundance, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  outbreak_Q1 <- outbreak_quartiles[1]  
  outbreak_Q2 <- outbreak_quartiles[2]  
  outbreak_Q3 <- outbreak_quartiles[3]  
  abundance_Q1 <- abundance_quartiles[1]  
  abundance_Q2 <- abundance_quartiles[2]  
  abundance_Q3 <- abundance_quartiles[3]  
  
  # if the outbreak is NA, encode as -1
  # if the weighted abundance is NA, encode as -1
  outbreak_abundance_dat$outbreak_freq[is.na(outbreak_abundance_dat$outbreak_freq)] <- -1
  outbreak_abundance_dat$sum_weighted_abundance[is.na(outbreak_abundance_dat$sum_weighted_abundance)] <- -1
  
  outbreak_abundance_dat$outbreak_freq_range <- 0
  outbreak_abundance_dat$abundance_range <- 0
  
  # recode the outbreak and weighted abundance (according to the quartile range)
  for (i in 1:nrow(outbreak_abundance_dat)){
    if(outbreak_abundance_dat$outbreak_freq[i]==-1){
      outbreak_abundance_dat$outbreak_freq_range[i] <- 0
    }else if(0 <= outbreak_abundance_dat$outbreak_freq[i] & outbreak_abundance_dat$outbreak_freq[i] <= outbreak_Q1){
      outbreak_abundance_dat$outbreak_freq_range[i] <- 1
    }else if(outbreak_Q1 < outbreak_abundance_dat$outbreak_freq[i] & outbreak_abundance_dat$outbreak_freq[i] <= outbreak_Q2){
      outbreak_abundance_dat$outbreak_freq_range[i] <- 2
    }else if(outbreak_Q2 < outbreak_abundance_dat$outbreak_freq[i] & outbreak_abundance_dat$outbreak_freq[i] <= outbreak_Q3){
      outbreak_abundance_dat$outbreak_freq_range[i] <- 3
    }else{
      outbreak_abundance_dat$outbreak_freq_range[i] <- 4
    }
  }
  
  for (i in 1:nrow(outbreak_abundance_dat)){
    if(outbreak_abundance_dat$sum_weighted_abundance[i]==-1){
      outbreak_abundance_dat$abundance_range[i] <- 0
    }else if(0 <= outbreak_abundance_dat$sum_weighted_abundance[i] & outbreak_abundance_dat$sum_weighted_abundance[i] <= abundance_Q1){
      outbreak_abundance_dat$abundance_range[i] <- 1
    }else if(abundance_Q1 < outbreak_abundance_dat$sum_weighted_abundance[i] & outbreak_abundance_dat$sum_weighted_abundance[i] <= abundance_Q2){
      outbreak_abundance_dat$abundance_range[i] <- 2
    }else if(abundance_Q2 < outbreak_abundance_dat$sum_weighted_abundance[i] & outbreak_abundance_dat$sum_weighted_abundance[i] <= abundance_Q3){
      outbreak_abundance_dat$abundance_range[i] <- 3
    }else{
      outbreak_abundance_dat$abundance_range[i] <- 4
    }
  }
  
  # round to the fourth decimal place
  outbreak_label1 <- paste0('minimum-', round(outbreak_Q1, 4))
  outbreak_label2 <- paste0(round(outbreak_Q1, 4), '-', round(outbreak_Q2, 4))
  outbreak_label3 <- paste0(round(outbreak_Q2, 4), '-', round(outbreak_Q3, 4))
  outbreak_label4 <- paste0(round(outbreak_Q3, 4), '-maximum')
  
  abundance_label1 <- paste0('minimum-', round(abundance_Q1, 4))
  abundance_label2 <- paste0(round(abundance_Q1, 4), '-', round(abundance_Q2, 4))
  abundance_label3 <- paste0(round(abundance_Q2, 4), '-', round(abundance_Q3, 4))
  abundance_label4 <- paste0(round(abundance_Q3, 4), '-maximum')
  
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
  
  # draw the distribution diagram of outbreak and weighted abundance
  outbreak_plot <-ggplot(data = outbreak_abundance_dat) +
    geom_sf(aes(fill = outbreak_freq_range), color = "black", size = 0.1) +  
    scale_fill_manual(
      values = c("#FFEEEE", "#FFCCCC", "#FFAAAA", "#FF6666", "#FF0000"),  
      name = "outbreak Range"
    ) +  
    labs(
      title = paste0('aiv outbreak in ', stage),
      subtitle = "Colored with Red-White Gradient",
      caption = "Data source: FAO epress-i +"
    ) +
    theme_minimal() +
    theme(legend.position = "right", legend.title = element_text(size = 10))
  
  abundance_plot <-ggplot(data = outbreak_abundance_dat) +
    geom_sf(aes(fill = abundance_range), color = "black", size = 0.1) +  
    scale_fill_manual(
      values = c("#FFEEEE", "#FFCCCC", "#FFAAAA", "#FF6666", "#FF0000"),  
      name = "abundance range"
    ) +  
    labs(
      title = paste0(birdname, ' weighted abundance in ', stage),
      subtitle = "Colored with Red-White Gradient",
      caption = "Data source: ebird checklist data"
    ) +
    theme_minimal() +
    theme(legend.position = "right", legend.title = element_text(size = 10))
  
  #combined_plot <- grid.arrange(
  #  outbreak_plot, abundance_plot, 
  #  ncol = 2, 
  #  widths = c(0.9, 1.1)  
  #)
  #dev.off()
  
  output_path <- file.path(bird_folder, paste0(stage, "_plot.png"))
  ggsave(output_path, plot = abundance_plot, width = 12, height = 6, dpi = 600, device = "png")
  
  if (file.exists(output_path)) {
    message("Graphic saved successfully！")
  } else {
    message("The graphic was not saved successfully, please check for errors!")
  }
  
  #return(list(outbreak_plot = outbreak_plot, abundance_plot = abundance_plot))
  return(abundance_plot)
}



# 5. Run the aiv analysis & save result -----
csv_files <- list.files(path = birdname_folder_path, pattern = "\\.csv$", full.names = FALSE)
all_bird_name_vector <- sub("\\.csv$", "", csv_files)                           # delete the .csv file extension, leaving only the bird name
print(all_bird_name_vector)

all_birds_all_stages_stage_table <- data.frame()
all_birds_all_stages_density_table <- data.frame()
all_birds_all_stages_stageabundance_table <- data.frame()
all_birds_single_stage_abundance_table <- data.frame()
all_birds_single_stage_density_table <- data.frame()

sub_bird_name_vector <- all_bird_name_vector[1:5]
for (birdname in sub_bird_name_vector){
  
  abundance <- read.csv(file.path(birdname_folder_path, paste0(birdname, '.csv')))
  weighted_abundance_result <- get_weighted_abundance(abundance)
  
  stage_merge_df <- merge_all_stages_data(all_df, weighted_abundance_result, EU_center_coordinate_df)
  stage_rbind_df <- rbind_all_stages_data(all_df, weighted_abundance_result, EU_center_coordinate_df)
  num_grid <- nrow(stage_merge_df)
  
  stage_rbind_df$stage <- relevel(stage_rbind_df$stage, ref = "stage2")
  gam_stage4 <- get_all_stage_model(birdname, stage_rbind_df, kval_init = 50, max_k = 150, step_k = 20)
  
  if (!is.null(gam_stage4)) {
    print(summary(gam_stage4))
    print(gam_stage4$method)
    print(gam.check(gam_stage4))
  }
  
  single_stage_gam_model <- get_single_stage_gam_model(birdname, stage_merge_df)
  #summary(single_stage_gam_model$gam_stage1)
  #summary(single_stage_gam_model$gam_stage2)
  #summary(single_stage_gam_model$gam_stage3)
  
  #single_stage_gam_model$gam_stage1$method
  #gam.check(single_stage_gam_model$gam_stage1)
  #k_check_result <- mgcv:::k.check(single_stage_gam_model$gam_stage1)
  #pval <- k_check_result[4]
  
  single_stage_abundance_table <- get_single_stage_abundance_table(birdname, single_stage_gam_model, num_grid)
  single_stage_density_table <- get_single_stage_density_table(birdname, single_stage_gam_model, num_grid) 
  all_birds_single_stage_abundance_table <- rbind(all_birds_single_stage_abundance_table,
                                                  single_stage_abundance_table)
  all_birds_single_stage_density_table <- rbind(all_birds_single_stage_density_table,
                                                single_stage_density_table)
  
  single_bird_all_stages_stage_table <- get_all_stages_stage_table(birdname, gam_stage4, num_grid)
  single_bird_all_stages_density_table <- get_all_stages_density_table(birdname, gam_stage4, num_grid)
  single_bird_all_stages_stageabundance_table <- get_all_stages_stageabundance_table(birdname, gam_stage4, num_grid)
  all_birds_all_stages_stage_table <- rbind(all_birds_all_stages_stage_table,
                                            single_bird_all_stages_stage_table)
  all_birds_all_stages_density_table <- rbind(all_birds_all_stages_density_table,
                                              single_bird_all_stages_density_table)
  all_birds_all_stages_stageabundance_table <- rbind(all_birds_all_stages_stageabundance_table,
                                                     single_bird_all_stages_stageabundance_table)
  
  print(paste0(birdname,' is done'))
}

rownames(all_birds_all_stages_stage_table) <- NULL
rownames(all_birds_all_stages_density_table) <- NULL
rownames(all_birds_all_stages_stageabundance_table) <- NULL
rownames(all_birds_single_stage_abundance_table) <- NULL
rownames(all_birds_single_stage_density_table) <- NULL

write.csv(all_birds_all_stages_stage_table,
          file = paste0(aiv_analysis_output_folder, outbreak_type, '_all_birds_all_stages_stage_', write_csv_date, '.csv'),
          row.names = FALSE)
write.csv(all_birds_all_stages_density_table,
          file = paste0(aiv_analysis_output_folder, outbreak_type, '_all_birds_all_stages_density_', write_csv_date, '.csv'),
          row.names = FALSE)
write.csv(all_birds_all_stages_stageabundance_table,
          file = paste0(aiv_analysis_output_folder, outbreak_type, '_all_birds_all_stages_stageabundance_', write_csv_date, '.csv'),
          row.names = FALSE)
write.csv(all_birds_single_stage_abundance_table,
          file = paste0(aiv_analysis_output_folder, outbreak_type, '_all_birds_single_stage_abundance_', write_csv_date, '.csv'),
          row.names = FALSE)
write.csv(all_birds_single_stage_density_table,
          file = paste0(aiv_analysis_output_folder, outbreak_type, '_all_birds_single_stage_density_', write_csv_date, '.csv'),
          row.names = FALSE)


# 6. Plot weighted abundance distribution -----
for (birdname in sub_bird_name_vector){
  abundance <- read.csv(file.path(birdname_folder_path, paste0(birdname, '.csv')))
  weighted_abundance_result <- get_weighted_abundance(abundance)
  
  stage_merge_df <- merge_all_stages_data(all_df, weighted_abundance_result, EU_center_coordinate_df)
  stage_rbind_df <- rbind_all_stages_data(all_df, weighted_abundance_result, EU_center_coordinate_df)
  stage_rbind_df$stage <- relevel(stage_rbind_df$stage, ref = "stage2")
  
  output_folder <- paste0(weighted_abundance_output_folder, outbreak_type)
  stage1_plots <- save_outbreak_and_abundance_image(EU.map, stage_rbind_df, stage='stage1', birdname, output_folder)
  stage2_plots <- save_outbreak_and_abundance_image(EU.map, stage_rbind_df, stage='stage2', birdname, output_folder)
  stage3_plots <- save_outbreak_and_abundance_image(EU.map, stage_rbind_df, stage='stage3', birdname, output_folder)
  
  print(paste0(birdname,' is done'))
}
