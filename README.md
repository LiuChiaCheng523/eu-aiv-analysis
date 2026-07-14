## 專案摘要

本研究專案透過機器學習預測 2021 年至 2022 年部分歐洲禽流感高風險帶原鳥種的時空相對豐度，共計 66 種。最終預測結果使用 Python Plotly 進行互動式可視化。使用者可點選不同鳥類名稱，比較不同空間單位的平均相對豐度，以及個別空間單位的月趨勢曲線變化。

整體研究結果重現了候鳥遷徙模式：冬季往南遷移過冬，夏季往北遷移避暑。模型同時校正觀察者偏誤，並以環境因子驅動預測結果。後續分析結合已知禽流感爆發案例與廣義加性模型，驗證多種鳥類，尤其是水鳥，在疫情爆發高峰期的冬季與春季具有顯著較高的風險。

## Project Summary

This research project uses machine learning to predict the spatiotemporal relative abundance of 66 high-risk avian influenza carrier bird species across parts of Europe from 2021 to 2022. The final prediction outputs are visualized through an interactive Python Plotly dashboard. Users can select different bird species to compare average relative abundance across spatial units and examine monthly trend curves within individual spatial units.

The overall results reproduce broad migratory bird movement patterns, with birds moving southward in winter for overwintering and northward in summer. The modeling workflow also corrects for observer bias and uses environmental variables to drive abundance predictions. Subsequent analyses combine known avian influenza outbreak records with generalized additive models to evaluate seasonal risk, showing that many bird species, especially waterbirds, have significantly higher risk during the winter and spring outbreak peaks.

![Plotly Dashboard example](docs/Plotly_Dashboard_example.png)

## 資料來源

整合多個公開資料來源，建立歐洲禽流感高風險鳥種之時空豐度與疫情風險分析資料庫。主要資料來源包含：

1. **eBird checklist**
   - Website: https://science.ebird.org/en/use-ebird-data/download-ebird-data-products
   - 作為鳥類觀測紀錄與相對豐度建模的主要資料來源。
   - 資料包含鳥種觀察紀錄、觀測時間、地理位置、觀察者資訊與觀測努力量等欄位。
   - 本研究使用 2021 至 2022 年期間之 eBird 鳥類觀測資料進行模型訓練與預測。

2. **eBird Status and Trends**
   - Website: https://science.ebird.org/status-and-trends
   - 用於驗證本研究所估計之鳥類相對豐度時空分布趨勢。
   - 該資料提供鳥種在不同時間與空間下的相對豐度、出現機率與分布範圍等資訊。

3. **FAO EMPRES-i**
   - Website: https://empres-i.apps.fao.org/
   - 作為禽流感疫情案例資料來源。
   - 使用 2021 至 2022 年歐洲地區家禽與野鳥之禽流感通報紀錄，分析鳥類豐度與疫情爆發風險之關聯。

4. **FAO Livestock Systems**
   - Website: https://data.apps.fao.org/
   - 作為家禽養殖密度資料來源。
   - 納入雞與鴨之單位面積養殖密度，作為禽流感風險模型中的控制變數與潛在風險因子。

5. **Google Earth Engine**
   - Website: https://earthengine.google.com/
   - API access: https://developers.google.com/earth-engine/guides/access
   - ERA5-Land: https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR
   - Dynamic World: https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_DYNAMICWORLD_V1
   - 本研究透過 Google Earth Engine 取得衛星遙測與氣候環境資料。
   - 使用 ERA5-Land Daily Aggregated 氣候資料，例如地表溫度、降雨量、風速與濕度等環境變數。
   - 使用 Dynamic World 土地覆蓋資料，取得水域、農地、森林、人造建地等土地覆蓋比例。
   - 本研究使用 **Google Earth Engine Python API** 下載與處理遙測資料；使用前需先申請並啟用 GEE API。
   - 下載資料時需準備對應研究範圍之 `.shp` 空間邊界檔，用於計算各網格內的多種遙測與氣候統計量。

6. **European 100 km × 100 km grid**
   - QGIS: https://www.qgis.org/
   - 本研究以歐洲地區 100 公里 × 100 公里網格作為主要空間分析單位。
   - 網格資料是利用 QGIS 對公開歐洲地圖 `.shp` 檔進行加工處理後建立，包括切割固定解析度網格、處理邊界區域、平滑化邊緣，以及降低 shapefile 檔案大小。
   - 所有鳥類觀測、環境因子、土地覆蓋、家禽密度與禽流感案例資料，皆依月份與網格進行時空配對。

## 研究方法架構

本研究方法架構可分為五個主要階段：資料整合與前處理、土地覆蓋缺失值補值、鳥類相對豐度建模、禽流感風險分析，以及互動式視覺化。

```text
1. Data Integration and Preprocessing
   資料整合與前處理
   │
   ├── eBird checklist bird observation records
   │   eBird 鳥類觀測紀錄
   │
   ├── ERA5-Land climate variables
   │   ERA5-Land 氣候環境因子
   │
   ├── Dynamic World land-cover variables
   │   Dynamic World 土地覆蓋因子
   │
   ├── FAO poultry density data
   │   FAO 家禽養殖密度資料
   │
   ├── FAO EMPRES-i avian influenza records
   │   FAO EMPRES-i 禽流感通報紀錄
   │
   └── European 100 km × 100 km grid
       歐洲 100 公里 × 100 公里網格
      
    - 本研究將不同來源、不同時間尺度與空間尺度的資料，統一整理成「月份 × 100 公里 × 100 公里網格」的分析資料表。
    - eBird 鳥類觀測紀錄依照觀測日期與座標，對應到特定月份與網格。
    - ERA5-Land 氣候資料與 Dynamic World 土地覆蓋資料，透過 Google Earth Engine Python API 計算各月份、各網格內的環境統計量。
    - FAO 家禽養殖密度資料與 FAO EMPRES-i 禽流感通報紀錄，依空間位置配對至相同網格，並依月份彙整。
    - 整合後的資料表以每一列代表一個「月份 × 網格」單位，作為後續土地覆蓋補值、鳥類相對豐度建模與禽流感風險分析的基礎。
        ↓

2. Dynamic World Land-cover Imputation
   Dynamic World 土地覆蓋缺失值補值
   │
   ├── Logit transformation for land-cover proportions
   │   土地覆蓋比例 logit 轉換
   │
   ├── XGBoost regression using climate variables and coordinates
   │   使用氣候因子與經緯度座標進行 XGBoost 迴歸補值
   │
   └── Validation using RMSE and MAE
       使用 RMSE 與 MAE 進行補值結果驗證
        ↓

3. Bird Relative Abundance Modeling
   鳥類相對豐度建模
   │
   ├── eBird checklist filtering
   │   eBird 觀測資料篩選
   │
   ├── Observer bias correction using GLMM
   │   使用 GLMM 校正觀察者偏誤
   │
   ├── Spatiotemporal sampling strategy
   │   時空抽樣策略
   │
   ├── XGBoost Poisson model for each risk bird species
   │   針對各風險鳥種建立 XGBoost Poisson 模型
   │   - 模型特徵包含 59 種 ERA5-Land 氣候因子、9 種 Dynamic World 土地覆蓋因子、eBird 觀測努力量，以及觀察者偏誤校正指標等。
   │
   ├── Model evaluation using SRC
   │   使用 SRC（Spearman Rank Correlation）評估模型預測表現
   │
   └── Monthly grid-level relative abundance prediction
       產生月份 × 網格層級之相對豐度預測
        ↓

4. Avian Influenza Risk Analysis
   禽流感風險分析
   │
   ├── Weighted relative abundance by outbreak stage
   │   依疫情階段計算加權相對豐度
   │
   ├── Integration with poultry and wild bird outbreak records
   │   結合家禽與野鳥禽流感通報紀錄
   │
   ├── Adjustment for chicken and duck density
   │   控制雞與鴨養殖密度
   │
   └── Generalized Additive Model analysis
       使用廣義加性模型分析疫情風險
        ↓

5. Interactive Visualization
   互動式視覺化
   │
   ├── Plotly Dash dashboard
   │   Plotly Dash 儀表板
   │
   ├── Species-level abundance map
   │   鳥種層級相對豐度地圖
   │
   └── Monthly trend curve by spatial unit
       個別空間單位月趨勢曲線
```

# CLI-based workflow

This repository provides five command-line entry points. Click each item to jump to the detailed CLI section below.

1. [Land-cover imputation](#1-land-cover-imputation-cli): optional workflow for regenerating the processed Dynamic World land-cover table.
2. [Bird abundance model training and prediction](#2-abundance-prediction-cli): trains species-level abundance models and exports monthly grid predictions.
3. [Validation / ensemble post-processing](#3-validation--prediction-summary-cli): combines iterations, computes validation metrics, applies MAD filtering, and creates summary outputs.
4. [AIV analysis](#4-aiv-analysis-cli): links bird abundance outputs with Wild or Domestic avian influenza outbreak records.
5. [Interactive Plotly dashboard](#5-plotly-dashboard): launches a local browser dashboard for exploring abundance maps and monthly trends.

## Scripts Download

To download the analysis scripts, clone this GitHub repository using Git:

```bash
git clone https://github.com/LiuChiaCheng523/eu-aiv-analysis.git
```

## Data Download

Download the required data from Google Drive:

[Google Drive folder](https://drive.google.com/drive/folders/1YX5561Yos3P4PsK1OfChdaQ7mAzx2Ihm)

Including:

1. `aiv_fixed_data`
2. `ebird_filtered_checklist`
3. `EU_100km_fishnet_simple_by_distance`
4. `gee_data`
5. `livestock_density_10km`

After downloading:

1. Extract the zip files.
2. Place the extracted folders under the project root.
3. Keep the internal relative folder structure unchanged.

## Project Structure

You may rename the repository folder after cloning or downloading it.

What must stay unchanged is the **internal relative folder structure**. The CLI scripts do not require the root folder to be named `abundance_r_test`, but they do expect sibling folders such as `gee_data/`, `EU_100km_fishnet_simple_by_distance/`, and `ebird_filtered_checklist/` to remain in the expected relative locations.

```text
eu-aiv-analysis/
├─ run_abundance_prediction.R
├─ run_validation_and_prediction_summary.R
├─ run_aiv_analysis.R
├─ run_land_cover_imputation.R
├─ Plotly_Dashboard.py
├─ bird_species_list.csv
├─ bird_type_lookup.csv
├─ README_plotly_dashboard.md
├─ gee_data/
│  ├─ era5_2016_2022/
│  ├─ land_cover_2016_2022.csv
│  └─ EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv
├─ EU_100km_fishnet_simple_by_distance/
├─ ebird_filtered_checklist/
├─ abundance_spatiotemporal_sampling_method/
├─ validation_data/
├─ validation_prediction_summary/
├─ aiv_fixed_data/
├─ livestock_density_10km/
├─ aiv_analysis/
└─ land_cover_imputation/
```

Folder notes:

- `gee_data/`: GEE-derived climate and land-cover data. `land_cover_2016_2022.csv` is the processed land-cover input used by the abundance model.
- `EU_100km_fishnet_simple_by_distance/`: European 100 km by 100 km grid shapefile.
- `ebird_filtered_checklist/`: filtered species-level eBird checklist CSV files.
- `abundance_spatiotemporal_sampling_method/`, `validation_data/`, `validation_prediction_summary/`, `aiv_analysis/`, and `land_cover_imputation/`: output folders generated by the CLI workflows.
- `.venv_plotly_dashboard/`: local Python virtual environment. It should usually not be uploaded to GitHub.

## Recommended Environments

### R

- Recommended version: `R 4.3.3`
- Tested version: `R 4.3.3`

All command examples in this README use `Rscript`.

If `Rscript` is not available in your terminal `PATH`, use the full executable path instead:

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" ".\run_land_cover_imputation.R" --help
```

Required R packages:

```r
install.packages(c("sf", "dplyr", "zoo", "xgboost", "pbapply", "data.table"))
```

Additional packages are required by the abundance / AIV scripts depending on the workflow. If a script reports a missing package, install it in the same R environment.

### Python

- Recommended version for dashboard: `Python 3.12`

Suggested virtual environment setup:

```powershell
py -3.12 -m venv ".venv_plotly_dashboard"
& ".\.venv_plotly_dashboard\Scripts\Activate.ps1"
python -m pip install --upgrade pip setuptools wheel
python -m pip install numpy pandas plotly dash geopandas shapely pyogrio fiona
```

## Workflow Order

The usual workflow is:

1. Optional: run [`run_land_cover_imputation.R`](#1-land-cover-imputation-cli) only if you need to regenerate `gee_data/land_cover_2016_2022.csv`.
2. Run [`run_abundance_prediction.R`](#2-abundance-prediction-cli) to train abundance models and export prediction iterations.
3. Run [`run_validation_and_prediction_summary.R`](#3-validation--prediction-summary-cli) to combine model iterations and export MAD-filtered abundance data.
4. Run [`run_aiv_analysis.R`](#4-aiv-analysis-cli) to analyze Wild or Domestic AIV associations.
5. Run [`Plotly_Dashboard.py`](#5-plotly-dashboard) to explore the processed abundance outputs interactively.

## 1. Land-cover Imputation CLI

Script: `run_land_cover_imputation.R`

This CLI is an optional preprocessing workflow. It takes a GEE-exported land-cover and climate CSV with missing Dynamic World land-cover values, creates validation samples, runs XGBoost imputation for selected land-cover classes, compares imputation performance, and aggregates multi-seed predictions into a final land-cover table.

You do **not** need to run this step if `gee_data/land_cover_2016_2022.csv` is already available, because that file is the processed land-cover table used by `run_abundance_prediction.R`.

Can it run without previous steps?

- Yes. This workflow is independent of the abundance, validation, and AIV CLIs.
- It still requires the GEE input CSV and the EU grid shapefile.

Input locations:

- `gee_data/EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv`: GEE-exported land-cover and climate table containing missing land-cover values.
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: EU 100 km grid shapefile.

Output locations and file patterns:

- `land_cover_imputation/sampling_data/`: validation samples, pattern `land_cover_seed*.csv`.
- `land_cover_imputation/ml_prediction_output/`: XGBoost prediction outputs by seed and land-cover type, pattern `seed*.csv`.
- `land_cover_imputation/two_method_performance/`: method comparison metrics, pattern `seed*.csv`.
- `land_cover_imputation/two_method_test_output/`: validation-set prediction tables, pattern `seed*.csv`.
- `land_cover_imputation/final_output/`: final imputed land-cover CSV.
- `land_cover_imputation/final_output/{land_cover_type}/`: final CSV for a single selected land-cover type, pattern `EU_2016_2022_land_cover_imputation_by_xgboost_{type}.csv`.

CLI parameters:

### Operation parameters

Operation parameters are the options users most commonly adjust for controlling the workflow.

- `--mode`: Select which part of the workflow to run.
  - Format / default: `sampling`, `imputation`, `aggregate`, or `all`. Default: `all`.

- `--seed-start`: First seed number.
  - Format / default: Integer. Default: `123`.

- `--seed-end`: Last seed number.
  - Format / default: Integer. Default: `222`.

- `--n-cores`: Number of CPU cores for parallel processing.
  - Format / default: Integer. Default: `2`.

- `--land-cover-types`: Select land-cover types to impute.
  - Format / default: Comma-separated values. Default: all 9 types: `bare,built,crops,flooded_vegetation,grass,shrub_and_scrub,snow_and_ice,trees,water`.

- `--help`: Show help message.
  - Format / default: No value.


### Path parameters

Path parameters are usually left unchanged. The CLI is designed to work with the default relative project structure, so changing these paths is only recommended when your data are stored outside the project folder.

- `--eu-shp-path`: Set the EU grid shapefile path.
  - Format / default: File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`.

- `--input-csv`: Set the input CSV with missing land-cover values.
  - Format / default: File path. Default: `gee_data/EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv`.

- `--output-folder`: Set the project root used for output folders.
  - Format / default: Folder path. Default: script folder / project root.

Simplest usage:

```powershell
Rscript ".\run_land_cover_imputation.R" --mode all
```

Example usage:

```powershell
Rscript ".\run_land_cover_imputation.R" `
  --mode all `
  --seed-start 123 `
  --seed-end 124 `
  --n-cores 2 `
  --land-cover-types bare
```

## 2. Abundance Prediction CLI

Script: `run_abundance_prediction.R`

This CLI trains bird relative abundance models for one, multiple, or all species. For each species, it reads the filtered eBird checklist, joins ERA5-Land climate variables, Dynamic World land-cover variables, the EU grid, temporal variables, geographic coordinates, eBird effort variables, and a GLMM-based observer correction index. It then trains XGBoost Poisson models over multiple iterations and exports abundance predictions, validation predictions, model performance, and feature importance outputs.

Can it run without previous steps?

- Yes, if the downloaded data already include `gee_data/land_cover_2016_2022.csv`.
- No, if the processed land-cover table is missing. In that case, run the optional land-cover imputation workflow first or provide a custom `--lc-path`.

Input locations:

- `ebird_filtered_checklist/`: species-level filtered eBird checklist CSVs. File pattern: `^{Scientific name}_filtered_2019to2022.csv` or `^{Scientific name}_filtered2019to2022.csv`.
- `gee_data/era5_2016_2022/`: yearly ERA5-Land environmental CSVs. File pattern: `{year}_median_combined_result.csv`.
- `gee_data/land_cover_2016_2022.csv`: processed Dynamic World land-cover table.
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: EU 100 km grid shapefile.

Output locations and file patterns:

- `glmm_performance/cci_random_effect/`: observer correction random effects, pattern `{species}.csv`.
- `glmm_performance/glmm_performance/`: GLMM performance summaries, pattern `{species}.csv`.
- `validation_data/`: held-out validation records, pattern `{species}.csv`.
- `abundance_spatiotemporal_sampling_method/abundance_prediction/{species}/`: monthly grid-level prediction iterations, pattern `abundance_iteration*.csv`.
- `abundance_spatiotemporal_sampling_method/validation_prediction/{species}/`: validation prediction iterations, pattern `abundance_iteration*.csv`.
- `abundance_spatiotemporal_sampling_method/feature_importance/{species}/`: feature importance by iteration, pattern `iteration*.csv`.
- `abundance_spatiotemporal_sampling_method/model_performance/`: species-level model performance, pattern `{species}.csv`.

CLI parameters:

### Operation parameters

Operation parameters are the options users most commonly adjust for choosing species, model settings, filtering rules, and reproducibility.

- `--species`: Run selected species.
  - Format / default: Scientific name, comma-separated names, or repeated values. Example: `"Anas crecca"` or `"Anas crecca,Ardea alba"`.

- `--species-file`: Run species listed in a TXT or CSV file.
  - Format / default: File path. CSV can contain `scientific_name`, `birdname`, or use the first column.

- `--all-species`: Run all species found in `ebird_filtered_checklist/`.
  - Format / default: No value.

- `--list-species`: List available species and exit.
  - Format / default: No value.

- `--start-year`: Start year for filtering.
  - Format / default: Integer. Default: `2021`.

- `--end-year`: End year for filtering.
  - Format / default: Integer. Default: `2022`.

- `--protocol`: Select eBird protocol(s).
  - Format / default: Comma-separated values. Default: `Traveling`.

- `--obs-quantile-cutoff`: Observation count outlier cutoff.
  - Format / default: Numeric quantile. Default: `0.99`.

- `--observer-quantile-cutoff`: Observer-level cutoff.
  - Format / default: Numeric quantile. Default: `0.7`.

- `--grid-quantile-cutoff`: Grid-level cutoff.
  - Format / default: Numeric quantile. Default: `0.5`.

- `--validation-ratio`: Validation split ratio.
  - Format / default: Numeric value between 0 and 1. Default: `0.1`.

- `--n-iter`: Number of model iterations.
  - Format / default: Integer. Default: `3`.

- `--seed`: Random seed.
  - Format / default: Integer. Default: `123`.

- `--help`: Show help message.
  - Format / default: No value.


### Path parameters

Path parameters are usually left unchanged. The CLI is designed to work with the default relative project structure, so changing these paths is only recommended when your data are stored outside the project folder.

- `--eu-shp-path`: Set EU grid shapefile path.
  - Format / default: File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`.

- `--env-folder`: Set ERA5-Land environmental CSV folder.
  - Format / default: Folder path. Default: `gee_data/era5_2016_2022/`.

- `--lc-path`: Set processed land-cover CSV path.
  - Format / default: File path. Default: `gee_data/land_cover_2016_2022.csv`.

- `--checklist-folder`: Set filtered eBird checklist folder.
  - Format / default: Folder path. Default: `ebird_filtered_checklist/`.

- `--output-folder`: Set output root folder.
  - Format / default: Folder path. Default: script folder / project root.

Simplest usage:

```powershell
Rscript ".\run_abundance_prediction.R" --all-species
```

Example usage:

```powershell
Rscript ".\run_abundance_prediction.R" `
  --species "Anas crecca" `
  --n-iter 5 `
  --obs-quantile-cutoff 0.98 `
  --observer-quantile-cutoff 0.7 `
  --grid-quantile-cutoff 0.5
```

## 3. Validation / Prediction Summary CLI

Script: `run_validation_and_prediction_summary.R`

This CLI post-processes abundance model outputs. It combines multiple prediction iterations, computes validation metrics such as SRC (Spearman Rank Correlation), applies MAD filtering to remove unstable prediction values, exports a cleaned abundance table for each species, summarizes land-cover feature importance, and generates validation and map figures.

Can it run without previous steps?

- No, in the usual workflow. It requires `run_abundance_prediction.R` outputs under `abundance_spatiotemporal_sampling_method/` and `validation_data/`.
- It can run only if equivalent abundance prediction folders and validation files already exist or are provided through custom paths.

Input locations:

- `abundance_spatiotemporal_sampling_method/abundance_prediction/{species}/`: model prediction iterations from the abundance CLI.
- `abundance_spatiotemporal_sampling_method/validation_prediction/{species}/`: validation prediction iterations.
- `abundance_spatiotemporal_sampling_method/feature_importance/{species}/`: feature importance files.
- `validation_data/`: held-out validation data.
- `ebird_filtered_checklist/`: checklist files used for sample summaries and species matching.
- `bird_type_lookup.csv`: bird group lookup table.
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: EU 100 km grid shapefile.

Output locations and file patterns:

- `validation_prediction_summary/validation_summary.csv`: species-level validation summary.
- `validation_prediction_summary/mad_filter_abundance/`: MAD-filtered abundance outputs, pattern `{species}.csv`.
- `validation_prediction_summary/bird_type_lookup_intersected.csv`: bird type lookup after intersecting with processed species.
- `validation_prediction_summary/feature_importance_land_cover_summary.csv`: land-cover feature-importance summary.
- `validation_prediction_summary/plots/validation_summary_scatter.png`: validation performance plot.
- `validation_prediction_summary/plots/feature_importance_land_cover_radar.png`: land-cover feature-importance radar chart.
- `validation_prediction_summary/plots/species_maps/{species}/`: species abundance maps.

CLI parameters:

### Operation parameters

Operation parameters are the options users most commonly adjust for choosing species, ensemble settings, MAD filtering, and map outputs.

- `--species`: Process selected species.
  - Format / default: Scientific name, comma-separated names, or repeated values.

- `--species-file`: Process species listed in a TXT or CSV file.
  - Format / default: File path. CSV can contain `scientific_name`, `birdname`, or use the first column.

- `--all-species`: Process all species available under `abundance_prediction/`.
  - Format / default: No value.

- `--list-species`: List species currently available for post-processing and exit.
  - Format / default: No value.

- `--n-iter`: Number of iterations to combine.
  - Format / default: Integer. Default: auto-detect.

- `--mad-k`: MAD filter multiplier.
  - Format / default: Numeric value. Default: `1.5`.

- `--map-years`: Years used for species map output.
  - Format / default: Comma-separated years. Default: auto-detect all years.

- `--map-months`: Months used for species map output.
  - Format / default: Comma-separated month numbers. Default: `1,4,7,10`.

- `--help`: Show help message.
  - Format / default: No value.


### Path parameters

Path parameters are usually left unchanged. The CLI is designed to work with the default relative project structure, so changing these paths is only recommended when your data are stored outside the project folder.

- `--abundance-root`: Set root folder created by `run_abundance_prediction.R`.
  - Format / default: Folder path. Default: `abundance_spatiotemporal_sampling_method/`.

- `--validation-folder`: Set validation data folder.
  - Format / default: Folder path. Default: `validation_data/`.

- `--checklist-folder`: Set checklist folder.
  - Format / default: Folder path. Default: `ebird_filtered_checklist/`.

- `--bird-type-csv`: Set bird type lookup CSV.
  - Format / default: File path. Default: `bird_type_lookup.csv`.

- `--eu-shp-path`: Set EU grid shapefile path.
  - Format / default: File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`.

- `--output-folder`: Set output root folder.
  - Format / default: Folder path. Default: script folder / project root.

Simplest usage:

Before running all species, you can first list species that already have outputs from `run_abundance_prediction.R`.

List available species:

```powershell
Rscript ".\run_validation_and_prediction_summary.R" --list-species
```

Process all available species:

```powershell
Rscript ".\run_validation_and_prediction_summary.R" --all-species
```

Example usage:

```powershell
Rscript ".\run_validation_and_prediction_summary.R" `
  --species "Vanellus vanellus,Branta canadensis" `
  --n-iter 50 `
  --mad-k 1.5
```

### Example Validation SRC Performance

![Validation SRC Performance example](docs/Validation_SRC_Performance_example.png)

### Example Land Cover Feature Importance Radar Chart

![Land Cover Feature Importance Radar Chart example](docs/Land_Cover_FI_Radar_Chart_example.png)

## 4. AIV Analysis CLI

Script: `run_aiv_analysis.R`

This CLI analyzes associations between predicted bird abundance and avian influenza outbreak records. It reads MAD-filtered abundance outputs, computes stage-level weighted abundance, joins Wild or Domestic outbreak counts, adjusts for chicken and duck density, runs GAM-based association analyses, and exports result tables plus overview and stage-level abundance maps.

Can it run without previous steps?

- No, in the usual workflow. It requires `validation_prediction_summary/mad_filter_abundance/`, which is generated by the validation / prediction summary CLI.
- It can run only if equivalent MAD-filtered abundance CSVs already exist or are provided through `--bird-abundance-folder`.

Input locations:

- `validation_prediction_summary/mad_filter_abundance/`: MAD-filtered abundance CSVs. File pattern: `{species}.csv`.
- `aiv_fixed_data/EU aiv fixed data 2021.csv`: 2021 AIV outbreak records.
- `aiv_fixed_data/EU aiv fixed data 2022.csv`: 2022 AIV outbreak records.
- `livestock_density_10km/chicken livestock density 10km.csv`: chicken density covariate.
- `livestock_density_10km/duck livestock density_2015_10km.csv`: duck density covariate.
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: EU 100 km grid shapefile.

Output locations and file patterns:

- `aiv_analysis/{Wild_or_Domestic}/{YYYYMMDD}/`: date-stamped AIV analysis output folder.
- `aiv_analysis/{Wild_or_Domestic}/{YYYYMMDD}/*.csv`: GAM result tables.
- `aiv_analysis/{Wild_or_Domestic}/{YYYYMMDD}/overview_maps/`: chicken density, duck density, outbreak map, and monthly outbreak distribution plots.
- `aiv_analysis/{Wild_or_Domestic}/{YYYYMMDD}/weighted_abundance/{species}/`: stage-level weighted abundance maps, pattern `stage*_plot.png`.

CLI parameters:

### Operation parameters

Operation parameters are the options users most commonly adjust for choosing outbreak type, species, and output date tags.

- `--outbreak-type`: Select outbreak type.
  - Format / default: Required. One of `Wild` or `Domestic`.

- `--species`: Analyze selected species.
  - Format / default: Scientific name, comma-separated names, or repeated values.

- `--species-file`: Analyze species listed in a TXT or CSV file.
  - Format / default: File path. CSV can contain `scientific_name`, `birdname`, or use the first column.

- `--all-species`: Analyze all species available in `mad_filter_abundance/`.
  - Format / default: No value.

- `--list-species`: List species currently available for AIV analysis and exit.
  - Format / default: No value.

- `--write-date`: Set output date tag.
  - Format / default: Format `YYYYMMDD`. Default: today's date.

- `--help`: Show help message.
  - Format / default: No value.


### Path parameters

Path parameters are usually left unchanged. The CLI is designed to work with the default relative project structure, so changing these paths is only recommended when your data are stored outside the project folder.

- `--eu-shp-path`: Set EU grid shapefile path.
  - Format / default: File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`.

- `--chicken-density-path`: Set chicken density CSV path.
  - Format / default: File path. Default: `livestock_density_10km/chicken livestock density 10km.csv`.

- `--duck-density-path`: Set duck density CSV path.
  - Format / default: File path. Default: `livestock_density_10km/duck livestock density_2015_10km.csv`.

- `--aiv-2021-path`: Set 2021 AIV fixed data path.
  - Format / default: File path. Default: `aiv_fixed_data/EU aiv fixed data 2021.csv`.

- `--aiv-2022-path`: Set 2022 AIV fixed data path.
  - Format / default: File path. Default: `aiv_fixed_data/EU aiv fixed data 2022.csv`.

- `--bird-abundance-folder`: Set MAD-filtered abundance folder.
  - Format / default: Folder path. Default: `validation_prediction_summary/mad_filter_abundance/`.

- `--output-folder`: Set output root folder.
  - Format / default: Folder path. Default: script folder / project root.

Simplest usage:

Before running all species, you can first list species that already have outputs from `run_abundance_prediction.R`.

List available species:

```powershell
Rscript ".\run_aiv_analysis.R" --list-species
```

Process all available species:

```powershell
Rscript ".\run_aiv_analysis.R" --outbreak-type Wild --all-species
```

Example usage:

```powershell
Rscript ".\run_aiv_analysis.R" `
  --outbreak-type Domestic `
  --species "Vanellus vanellus,Branta canadensis"
```

### Example GAM Model Abundance Effect And P-value Bar Chart

![GAM model abundance effect and p-value bar chart example](docs/abundance_effect_pvalue_plot.png)

### Example Vanellus vanellus Stage 2 Weighted Abundance Distribution

![Vanellus vanellus stage 2 weighted abundance distribution example](docs/stage2_plot.png)

## 5. Plotly Dashboard

Script: `Plotly_Dashboard.py`

This Python script starts a local Plotly Dash app for exploring the MAD-filtered abundance outputs. The dashboard builds a species dropdown from the intersection of available checklist files and abundance files, displays mean abundance across the EU 100 km grid, and allows users to click individual grid cells to inspect monthly abundance trends and summary metrics.

Can it run without previous steps?

- No, in the usual workflow. It requires `validation_prediction_summary/mad_filter_abundance/`, which is generated by the validation / prediction summary CLI.
- It can run only if equivalent MAD-filtered abundance CSVs already exist or are provided through `--abundance-folder`.

Input locations:

- `validation_prediction_summary/mad_filter_abundance/`: MAD-filtered abundance CSVs. File pattern: `{species}.csv`.
- `ebird_filtered_checklist/`: checklist CSVs used for species matching and sample counts.
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: EU 100 km grid shapefile.

Output / runtime behavior:

- Local dashboard URL: `http://127.0.0.1:8050`.
- No permanent output files are written by default.
- The dashboard runs only while the Python process is active.

CLI parameters:

### Operation parameters

Operation parameters are the options users most commonly adjust for choosing the initial species and local dashboard server settings.

- `--species`: Set initial species shown in the dashboard.
  - Format / default: Scientific name. Default: first available species.

- `--host`: Set Dash host.
  - Format / default: Host string. Default: `127.0.0.1`.

- `--port`: Set Dash port.
  - Format / default: Integer. Default: `8050`.

- `--debug`: Enable Dash debug mode.
  - Format / default: No value. Default: off.


### Path parameters

Path parameters are usually left unchanged. The dashboard is designed to work with the default relative project structure, so changing these paths is only recommended when your data are stored outside the project folder.

- `--eu-shp-path`: Set EU grid shapefile path.
  - Format / default: File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`.

- `--checklist-folder`: Set checklist folder.
  - Format / default: Folder path. Default: `ebird_filtered_checklist/`.

- `--abundance-folder`: Set MAD-filtered abundance folder.
  - Format / default: Folder path. Default: `validation_prediction_summary/mad_filter_abundance/`.

Simplest usage:

```powershell
python ".\Plotly_Dashboard.py"
```

Example usage:

```powershell
python ".\Plotly_Dashboard.py" --species "Vanellus vanellus" --host 127.0.0.1 --port 8051
```
