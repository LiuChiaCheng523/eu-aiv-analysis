## 專案摘要

本研究專案透過機器學習預測 2021 年至 2022 年部分歐洲禽流感高風險帶原鳥種的時空相對豐度，共計 67 種。最終預測結果使用 Python Plotly 進行互動式可視化。使用者可點選不同鳥類名稱，比較不同空間單位的平均相對豐度，以及個別空間單位的月趨勢曲線變化。

整體研究結果重現了候鳥遷徙模式：冬季往南遷移過冬，夏季往北遷移避暑。模型同時校正觀察者偏誤，並以環境因子驅動預測結果。後續分析結合已知禽流感爆發案例與廣義加性模型，驗證多種鳥類，尤其是水鳥，在疫情爆發高峰期的冬季與春季具有顯著較高的風險。

## Project Summary

This research project uses machine learning to predict the spatiotemporal relative abundance of 67 high-risk avian influenza carrier bird species across parts of Europe from 2021 to 2022. The final prediction outputs are visualized through an interactive Python Plotly dashboard. Users can select different bird species to compare average relative abundance across spatial units and examine monthly trend curves within individual spatial units.

The overall results reproduce broad migratory bird movement patterns, with birds moving southward in winter for overwintering and northward in summer. The modeling workflow also corrects for observer bias and uses environmental variables to drive abundance predictions. Subsequent analyses combine known avian influenza outbreak records with generalized additive models to evaluate seasonal risk, showing that many bird species, especially waterbirds, have significantly higher risk during the winter and spring outbreak peaks.

![Plotly Dashboard example](docs/Plotly_Dashboard_example.png)

## 資料來源

整合多個公開資料來源，建立歐洲禽流感高風險鳥種之時空豐度與疫情風險分析資料庫。主要資料來源包含：

1. **eBird checklist**
   - Website: https://science.ebird.org/en/use-ebird-data/download-ebird-data-products
   - 作為鳥類觀測紀錄與相對豐度建模的主要資料來源。
   - 資料包含鳥種觀察紀錄、觀測時間、地理位置、觀察者資訊與觀測努力量等欄位。
   - 本研究使用 2021 至 2022 年期間之 eBird 鳥類觀測資料進行模型擬合。

2. **eBird Status and Trends**
   - Website: https://science.ebird.org/status-and-trends
   - 作為外部資料，可用於驗證本研究所估計之鳥類相對豐度時空分布趨勢。
   - 該資料提供鳥種在不同時間與空間下的相對豐度、出現機率與分布範圍等資訊。

3. **FAO EMPRES-i**
   - Website: https://empres-i.apps.fao.org/
   - 作為禽流感疫情案例資料來源。
   - 使用 2021 至 2022 年歐洲地區家禽與野鳥之禽流感通報紀錄，分析估計的鳥類豐度與疫情爆發風險之關聯。

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
   - 本研究使用 **Python** 下載與處理遙測資料；使用前需先使用 Google 帳號申請 GEE 專案並啟用。
   - 下載資料時需準備對應研究範圍之 `.shp` 空間邊界檔，用於計算各網格內的多種遙測與氣候統計量。

6. **European 100 km × 100 km grid**
   - QGIS: https://www.qgis.org/
   - 本研究以歐洲地區 100 公里 × 100 公里網格作為主要空間分析單位。
   - 網格資料是利用 QGIS 對歐洲地圖 `.shp` 檔進行加工處理後建立，包括切割固定解析度網格、處理邊界區域、平滑化邊緣，以及降低 shapefile 檔案大小。
   - 所有鳥類觀測、環境因子、土地覆蓋、家禽密度與禽流感案例資料，皆依月份與網格進行時空配對。

## 研究方法架構

本研究方法架構可分為五個主要階段：資料整合與前處理、土地覆蓋缺失值補值、鳥類相對豐度建模、禽流感風險分析，以及模型結果視覺化。

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
    - ERA5-Land 氣候資料與 Dynamic World 土地覆蓋資料，分別透過 Python 腳本[ERA5-Land Climate Data Download](#1-era5-land-climate-data-download) 與 [Dynamic World Land-cover Data Download](#2-dynamic-world-land-cover-data-download)下載，並計算 Google Earth Engine 所提供之各月份、各網格內的環境統計量。
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
   ├── GPBoost model for each risk bird species
   │   針對各風險鳥種建立 GPBoost Negative Binomial 模型
   │   - 模型特徵包含 59 種 ERA5-Land 氣候因子、9 種 Dynamic World 土地覆蓋因子、eBird 觀察者與觀測努力量等。
   │
   ├── Model evaluation using SRC
   │   使用 SRC（Spearman Rank Correlation）評估模型表現
   │
   └── Monthly grid-level relative abundance prediction
       重建月份 × 網格層級之相對豐度
       - 將 observer random effect 與努力量標準化為相同參考水準，使不同月份與網格的相對豐度數值可相互比較。
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
   ├── Adjustment for chicken and duck livestock density and geographic coordinates as covariates
   │   將雞與鴨養殖密度及經緯度座標作為共變數進行校正
   │
   └── Generalized Additive Model analysis
       使用廣義加性模型分析疫情風險
        ↓

5. Model Results Visualization
   模型結果視覺化
   │
   ├── Relative abundance maps and model performance plots
   │   相對豐度地圖與模型表現圖
   │
   ├── Land-cover feature importance visualization
   │   土地覆蓋變數重要性視覺化
   │
   ├── AIV outbreak and weighted abundance maps
   │   禽流感通報與加權相對豐度分布圖
   │
   └── GAM effect and poultry-density result plots
       GAM 效應與家禽密度結果圖
```

# Analysis workflow and scripts

This repository contains notebooks, command-line scripts, and fixed-path analysis scripts. The workflow is not entirely CLI-based, so this section describes each script by its role, adjustable parameters, required inputs, outputs, and model usage.

The usual order is:

1. Download climate data with `01-ERA5-download-colab.ipynb`.
2. Download land-cover data with `02-DynamicWorld-download-colab.ipynb`.
3. If needed, impute missing land-cover values with `03-run_land_cover_imputation.R`.
4. Fitting of bird relative abundance models with `05-gpboost-cli.py`.
5. Summarize and visualize model outputs with `06-model-visualization-cli.py`.
6. Run AIV outbreak association analysis with `07-aiv-outbreak-analysis.R`.
7. Summarize AIV model results with `08-analysis-visualization-cli.py`.

`04-gpboost-val.py` is a Colab-derived GPBoost validation / batch experiment script.

## Scripts Download

To download the analysis scripts, clone this GitHub repository using Git:

```bash
git clone https://github.com/LiuChiaCheng523/eu-aiv-analysis.git
```

## Data Download

Download the required data from Google Drive:

[Google Drive folder](https://drive.google.com/drive/folders/1YX5561Yos3P4PsK1OfChdaQ7mAzx2Ihm)

Download and extract these zip files into the same project folder as the scripts:

1. `gee_data.zip`
   - Contains Google Earth Engine outputs used by the later modeling steps.
   - Includes ERA5-Land climate tables from `01-ERA5-download-colab.ipynb`.
   - Includes Dynamic World land-cover tables from `02-DynamicWorld-download-colab.ipynb`.
   - Includes land-cover imputation outputs from `03-run_land_cover_imputation.R`, which fills missing Dynamic World land-cover values.
   - Includes the processed / imputed `land_cover_2016_2022.csv` used by `05-gpboost-cli.py`, so users can start from step 05 without rerunning steps 01-03.
2. `ebird_filtered_checklist.zip`
   - Contains species-level filtered eBird checklist CSV files used for GPBoost abundance modeling.
3. `EU_100km_fishnet_simple_by_distance.zip`
   - Contains the European 100 km by 100 km grid shapefile and its sidecar files.
4. `aiv_fixed_data.zip`
   - Contains cleaned 2021 and 2022 avian influenza outbreak records.
5. `livestock_density_10km.zip`
   - Contains chicken and duck livestock density covariates.

After downloading:

1. Extract the zip files.
2. Place the extracted folders under the project root.
3. Keep the internal relative folder structure unchanged, as shown in [Project Structure](#project-structure).

In the default input examples below, `your path/eu-aiv-analysis/` means the full location of this project on your computer, for example `C:/Users/yourname/eu-aiv-analysis/`.

## Project Structure

The repository folder can be renamed, but the internal data-folder names should stay the same because several scripts read files by relative folder names.

```text
eu-aiv-analysis/
├─ 01-ERA5-download-colab.ipynb
├─ 02-DynamicWorld-download-colab.ipynb
├─ 03-run_land_cover_imputation.R
├─ 04-gpboost-val.py
├─ 05-gpboost-cli.py
├─ 06-model-visualization-cli.py
├─ 07-aiv-outbreak-analysis.R
├─ 08-analysis-visualization-cli.py
├─ install-r-packages.R
├─ r-packages.txt
├─ bird_type_lookup.csv
├─ gee_data/
│  ├─ era5_2016_2022/
│  ├─ land_cover_2016_2022.csv
│  └─ EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv
├─ EU_100km_fishnet_simple_by_distance/
├─ ebird_filtered_checklist/
├─ livestock_density_10km/
├─ aiv_fixed_data/
├─ land_cover_imputation/          # generated by 03-run_land_cover_imputation.R
├─ gpboost_abundance_outputs/      # generated by 05-gpboost-cli.py
├─ all-birds-abd/                  # generated by 06-model-visualization-cli.py
├─ weighted_abundance/             # generated by 07-aiv-outbreak-analysis.R
├─ aiv_analysis/                   # generated by 07-aiv-outbreak-analysis.R
├─ plot/                           # generated by 06 or 08 visualization scripts
└─ table/                          # generated by 06 or 08 visualization scripts
```

Folder notes:

- `install-r-packages.R`: R package installation script. It reads package names from `r-packages.txt`.
- `bird_type_lookup.csv`: bird species lookup table used to label waterbird and non-waterbird groups for model-result summaries and visualizations.
- `gee_data/`: Google Earth Engine outputs, including ERA5-Land climate tables and Dynamic World land-cover tables.
- `gee_data/EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv`: GEE table with missing land-cover values, used only if rerunning imputation.
- `EU_100km_fishnet_simple_by_distance/`: European 100 km by 100 km grid shapefile.
- `ebird_filtered_checklist/`: species-level filtered eBird checklist CSV files.
- `aiv_fixed_data/`: cleaned FAO EMPRES-i outbreak records for 2021 and 2022.
- `livestock_density_10km/`: chicken and duck density tables used in the AIV association analysis.


## Recommended Environments

### Python

Recommended version: `Python 3.10+`

Main Python packages:

```powershell
python -m pip install pandas numpy matplotlib seaborn geopandas shapely scikit-learn scipy xgboost gpboost
```

The Google Earth Engine notebooks also require notebook / Colab access and Earth Engine authentication.

### R

Recommended version: `R 4.3.3`

Install the required R packages with `install-r-packages.R`. The script reads package names from `r-packages.txt` and installs only missing packages.

```powershell
cd "your path/eu-aiv-analysis"
Rscript ".\install-r-packages.R"
```

## Quick Start

If the required Google Drive data folders are already placed under the project root, the local workflow can start from `05-gpboost-cli.py`. This is because `01-ERA5-download-colab.ipynb`, `02-DynamicWorld-download-colab.ipynb`, and `03-run_land_cover_imputation.R` have already been run, and their outputs are provided in `gee_data.zip`.

A small local test can start from one species:

```powershell
cd "your path/eu-aiv-analysis"
python ".\05-gpboost-cli.py" --birds "Anas crecca"
```

To run all checklist species:

```powershell
cd "your path/eu-aiv-analysis"
python ".\05-gpboost-cli.py" --birds all
python ".\06-model-visualization-cli.py"
```

Run `07-aiv-outbreak-analysis.R` after checking and editing the hard-coded input and output paths inside the script. Then use the same `write_csv_date` value in `08-analysis-visualization-cli.py`:

```powershell
cd "your path/eu-aiv-analysis"
Rscript ".\07-aiv-outbreak-analysis.R"
python ".\08-analysis-visualization-cli.py" --csv-date 20260811
```

## 1. ERA5-Land Climate Data Download

Script: `01-ERA5-download-colab.ipynb`

This notebook downloads ERA5-Land climate variables from Google Earth Engine and summarizes them by month and European 100 km grid cell. The output is used as climate input for bird abundance modeling.

Inputs:

- Google Earth Engine account and authentication.
- European 100 km grid shapefile.
- Study years and month range configured inside the notebook.

Default input path if placed with the notebook:

- `your path/eu-aiv-analysis/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: 100 km grid used to summarize ERA5-Land climate values by grid cell.

Adjustable settings:

- Study years and months.
- ERA5-Land variables to export.
- Google Drive / local output location.
- Grid or region shapefile path.

Outputs:

- ERA5-Land CSV files under `your path/eu-aiv-analysis/gee_data/era5_2016_2022/`.
- Expected file pattern: `your path/eu-aiv-analysis/gee_data/era5_2016_2022/{year}_median_combined_result.csv`.
- These outputs are already included in `gee_data.zip`.

Model used:

- No statistical model is trained here. This step extracts and summarizes climate covariates from Google Earth Engine.

## 2. Dynamic World Land-cover Data Download

Script: `02-DynamicWorld-download-colab.ipynb`

This notebook downloads Dynamic World land-cover probabilities from Google Earth Engine and summarizes land-cover proportions by month and European 100 km grid cell.

Inputs:

- Google Earth Engine account and authentication.
- European 100 km grid shapefile.
- Study years and months configured inside the notebook.

Default input path if placed with the notebook:

- `your path/eu-aiv-analysis/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: 100 km grid used to summarize Dynamic World land-cover probabilities by grid cell.

Adjustable settings:

- Study years and months.
- Dynamic World land-cover classes.
- Google Drive / local output location.
- Grid or region shapefile path.

Outputs:

- Raw or processed land-cover CSV files under `your path/eu-aiv-analysis/gee_data/`.
- If missing values remain, use `03-run_land_cover_imputation.R` to create the final `your path/eu-aiv-analysis/gee_data/land_cover_2016_2022.csv`.
- These outputs are already included in `gee_data.zip`.

Model used:

- No statistical model is trained in the download notebook. It prepares land-cover covariates for later modeling.

## 3. Land-cover Imputation

Script: `03-run_land_cover_imputation.R`

This optional CLI script fills missing Dynamic World land-cover values. It uses climate variables and coordinates to predict missing land-cover proportions, compares imputation performance, and writes final imputed land-cover tables.

Can it run independently?

- Yes. It does not require the abundance model outputs.
- It does require the GEE land-cover / climate CSV and the EU grid shapefile.

Inputs:

- `gee_data/EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv`
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`

Default input paths if the data folders are placed with the script:

- `your path/eu-aiv-analysis/gee_data/EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv`: GEE-exported land-cover and climate table with missing land-cover values.
- `your path/eu-aiv-analysis/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: European 100 km grid shapefile.

Adjustable parameters:

- `--mode`: choose `sampling`, `imputation`, `aggregate`, or `all`.
- `--seed-start` and `--seed-end`: control repeated sampling / imputation seeds.
- `--n-cores`: number of CPU cores.
- `--land-cover-types`: land-cover classes to impute, for example `water,crops,trees`.
- `--eu-shp-path`, `--input-csv`, and `--output-folder`: custom paths.

Outputs:

- `your path/eu-aiv-analysis/land_cover_imputation/sampling_data/`: validation samples.
- `your path/eu-aiv-analysis/land_cover_imputation/ml_prediction_output/`: model prediction outputs.
- `your path/eu-aiv-analysis/land_cover_imputation/two_method_performance/`: RMSE / MAE comparison tables.
- `your path/eu-aiv-analysis/land_cover_imputation/final_output/`: final imputed land-cover CSVs.
- The final processed land-cover table used later is provided in `gee_data.zip` as `your path/eu-aiv-analysis/gee_data/land_cover_2016_2022.csv`.

Model used:

- XGBoost regression is used for land-cover imputation.
- RMSE and MAE are used to compare imputation performance.

Example:

```powershell
cd "your path/eu-aiv-analysis"
Rscript ".\03-run_land_cover_imputation.R" --mode all --seed-start 123 --seed-end 124 --n-cores 2 --land-cover-types water
```

## 4. GPBoost Abundance Modeling

Script: `05-gpboost-cli.py`

This CLI script trains bird relative-abundance models. It combines eBird checklist records with ERA5-Land climate data, Dynamic World land-cover data, geographic grid information, observation effort variables, and observer IDs. For each selected species, it predicts monthly grid-level relative abundance for 2021 and 2022.

Inputs:

- `ebird_filtered_checklist/`: filtered eBird checklist CSVs.
- `gee_data/land_cover_2016_2022.csv`: processed land-cover table.
- `gee_data/era5_2016_2022/2021_median_combined_result.csv`
- `gee_data/era5_2016_2022/2022_median_combined_result.csv`

Default input paths if the data folders are placed with the script:

- `your path/eu-aiv-analysis/ebird_filtered_checklist/`: filtered eBird observation records for each species.
- `your path/eu-aiv-analysis/gee_data/land_cover_2016_2022.csv`: final Dynamic World land-cover table after preprocessing / imputation.
- `your path/eu-aiv-analysis/gee_data/era5_2016_2022/2021_median_combined_result.csv`: 2021 monthly ERA5-Land climate variables by grid.
- `your path/eu-aiv-analysis/gee_data/era5_2016_2022/2022_median_combined_result.csv`: 2022 monthly ERA5-Land climate variables by grid.

Adjustable parameters:

- `--birds`: species to run. Use scientific names such as `"Anas crecca"`, or use `all`.
- `--path-folder`: parent folder containing `gee_data/` and `ebird_filtered_checklist/`.
- `--output-folder`: parent folder for model outputs.

Important code settings that can be edited inside the script:

- eBird protocol filter currently keeps Traveling checklists.
- Observation-count outlier cutoff uses the 99th percentile.
- Observer filter uses the 70th percentile of observer checklist count.
- Grid filter uses the 50th percentile of grid checklist count.
- GPBoost settings include `learning_rate`, `num_leaves`, `min_data_in_leaf`, `feature_fraction`, `lambda_l2`, and `num_boost_round`.

Outputs:

- `your path/eu-aiv-analysis/gpboost_abundance_outputs/{species}/{species}_reconstruction_abundance.csv`
- `your path/eu-aiv-analysis/gpboost_abundance_outputs/{species}/{species}_feature_importance_gain.csv`
- `your path/eu-aiv-analysis/gpboost_abundance_outputs/{species}/{species}_model_performance.csv`
- `your path/eu-aiv-analysis/gpboost_abundance_outputs/{species}/{species}_random_effect.csv`
- `your path/eu-aiv-analysis/gpboost_abundance_outputs/{species}/{species}_random_effect_cov_pars.csv`

Model used:

- GPBoost model with a negative-binomial likelihood.
- Observer ID is modeled as a grouped random effect to reduce observer-bias effects.
- Model performance is summarized with Spearman Rank Correlation and MAE.

Example:

```powershell
cd "your path/eu-aiv-analysis"
python ".\05-gpboost-cli.py" --birds "Anas crecca" --path-folder "." --output-folder "."
```

## 5. Abundance Output Visualization and Summary

Script: `06-model-visualization-cli.py`

This CLI script reads GPBoost abundance outputs, exports simplified abundance tables, creates species maps, summarizes land-cover feature importance, and plots model performance.

Inputs:

- `gpboost_abundance_outputs/`: output folder from `05-gpboost-cli.py`.
- `gee_data/land_cover_2016_2022.csv`
- `gee_data/era5_2016_2022/`
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`
- `bird_type_lookup.csv`
- `ebird_filtered_checklist/`

Default input paths if the data folders are placed with the script:

- `your path/eu-aiv-analysis/gpboost_abundance_outputs/`: species-level abundance, feature-importance, performance, and observer-random-effect outputs from `05-gpboost-cli.py`.
- `your path/eu-aiv-analysis/gee_data/land_cover_2016_2022.csv`: land-cover covariates used for plotting and consistency checks.
- `your path/eu-aiv-analysis/gee_data/era5_2016_2022/`: ERA5-Land climate tables.
- `your path/eu-aiv-analysis/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: grid geometry used to draw maps.
- `your path/eu-aiv-analysis/bird_type_lookup.csv`: waterbird / non-waterbird lookup table.
- `your path/eu-aiv-analysis/ebird_filtered_checklist/`: checklist files used to match species names and summarize samples.

Adjustable parameters:

- `--bird-folder-path`: folder containing GPBoost outputs.
- `--main-folder`: project folder containing data folders.

Outputs:

- `your path/eu-aiv-analysis/all-birds-abd/{species}.csv`: simplified abundance table for each species.
- `your path/eu-aiv-analysis/gpboost_abundance_outputs/{species}/plot_python/*.png`: abundance maps for three seasonal periods.
- `your path/eu-aiv-analysis/plot/FI_radar.png`: land-cover feature-importance radar plot.
- `your path/eu-aiv-analysis/plot/SRC_scatter.png`: model-performance scatter plot.
- `your path/eu-aiv-analysis/table/checklist_and_model_info.csv`: checklist and model summary table.

Model used:

- No new model is trained here.
- The script summarizes outputs from the GPBoost abundance models and uses feature-importance rankings.

Example:

```powershell
cd "your path/eu-aiv-analysis"
python ".\06-model-visualization-cli.py" --bird-folder-path ".\gpboost_abundance_outputs" --main-folder "."
```

### Example Validation SRC Performance

![Validation SRC Performance example](docs/Validation_SRC_Performance_example.png)

### Example Land Cover Feature Importance Radar Chart

![Land Cover Feature Importance Radar Chart example](docs/Land_Cover_FI_Radar_Chart_example.png)

## 6. AIV Outbreak Association Analysis

Script: `07-aiv-outbreak-analysis.R`

This R script links predicted bird abundance with avian influenza outbreak records. It computes outbreak stages, neighboring-grid outbreak indicators, poultry density covariates, and weighted bird abundance by stage. It then fits association models for each bird species.

Important note:

- This script is currently not a CLI script. Input and output paths are set directly near the top of the file.
- Before running it, edit paths such as `EU_shp_path`, `chicken_density_path`, `duck_density_path`, `EU_aiv_2021_path`, `EU_aiv_2022_path`, `birdname_folder_path`, `write_csv_date`, and output folders.
- `outbreak_type` is also set inside the script as `Domestic` or `Wild`.

Inputs:

- `all-birds-abd/`: simplified abundance CSVs from `06-model-visualization-cli.py`.
- `aiv_fixed_data/EU aiv fixed data 2021.csv`
- `aiv_fixed_data/EU aiv fixed data 2022.csv`
- `livestock_density_10km/`: chicken and duck density CSVs.
- `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`

Default input paths if the data folders are placed with the script:

- `your path/eu-aiv-analysis/all-birds-abd/`: simplified monthly bird abundance tables created by `06-model-visualization-cli.py`.
- `your path/eu-aiv-analysis/aiv_fixed_data/EU aiv fixed data 2021.csv`: cleaned 2021 AIV outbreak records.
- `your path/eu-aiv-analysis/aiv_fixed_data/EU aiv fixed data 2022.csv`: cleaned 2022 AIV outbreak records.
- `your path/eu-aiv-analysis/livestock_density_10km/chicken_livestock_density_10km.csv`: chicken density by grid.
- `your path/eu-aiv-analysis/livestock_density_10km/duck_livestock_density_2015_10km.csv`: duck density by grid.
- `your path/eu-aiv-analysis/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`: grid geometry and grid IDs.

Adjustable code settings:

- `outbreak_type`: choose `Domestic` or `Wild`.
- `write_csv_date`: date tag used in output filenames.
- Input and output folder paths.
- Stage definitions currently separate 2021/05-2021/10, 2021/11-2022/04, and 2022/05-2022/10.

Outputs:

- `your path/eu-aiv-analysis/aiv_analysis/{outbreak_type}_all_birds_all_stages_stage_{write_csv_date}.csv`
- `your path/eu-aiv-analysis/aiv_analysis/{outbreak_type}_all_birds_all_stages_density_{write_csv_date}.csv`
- `your path/eu-aiv-analysis/aiv_analysis/{outbreak_type}_all_birds_all_stages_stageabundance_{write_csv_date}.csv`
- `your path/eu-aiv-analysis/aiv_analysis/{outbreak_type}_all_birds_single_stage_abundance_{write_csv_date}.csv`
- `your path/eu-aiv-analysis/aiv_analysis/{outbreak_type}_all_birds_single_stage_density_{write_csv_date}.csv`
- `your path/eu-aiv-analysis/weighted_abundance/{outbreak_type}/{species}/stage1_plot.png`
- `your path/eu-aiv-analysis/weighted_abundance/{outbreak_type}/{species}/stage2_plot.png`
- `your path/eu-aiv-analysis/weighted_abundance/{outbreak_type}/{species}/stage3_plot.png`
- Overview map paths are controlled by the variables near the top of `07-aiv-outbreak-analysis.R`, such as `chicken_density_output_path`, `duck_density_output_path`, `domestic_outbreak_output_path`, and `wild_outbreak_output_path`.

Model used:

- Generalized Additive Models from `mgcv` are used to estimate associations between outbreak outcomes, weighted bird abundance, and poultry density covariates.
- Chicken and duck densities are included as control / explanatory variables for domestic outbreak analysis.

Example:

```powershell
cd "your path/eu-aiv-analysis"
Rscript ".\07-aiv-outbreak-analysis.R"
```

### Example GAM Model Abundance Effect And P-value Bar Chart

![GAM model abundance effect and p-value bar chart example](docs/abundance_effect_pvalue_plot.png)

### Example Vanellus vanellus Stage 2 Weighted Abundance Distribution

![Vanellus vanellus stage 2 weighted abundance distribution example](docs/stage2_plot.png)

## 7. AIV Result Visualization

Script: `08-analysis-visualization-cli.py`

This CLI script summarizes AIV model output tables. It counts bird groups with positive and significant abundance effects, plots Wild and Domestic abundance-effect histograms, and visualizes chicken / duck livestock density effects.

Inputs:

- `bird_type_lookup.csv`
- `aiv_analysis/Wild_all_birds_single_stage_abundance_{csv_date}.csv`
- `aiv_analysis/Domestic_all_birds_single_stage_abundance_{csv_date}.csv`
- `aiv_analysis/Domestic_all_birds_single_stage_density_{csv_date}.csv`

Default input paths if the data folders are placed with the script:

- `your path/eu-aiv-analysis/bird_type_lookup.csv`: waterbird / non-waterbird lookup table.
- `your path/eu-aiv-analysis/aiv_analysis/Wild_all_birds_single_stage_abundance_{csv_date}.csv`: Wild outbreak abundance-effect GAM results.
- `your path/eu-aiv-analysis/aiv_analysis/Domestic_all_birds_single_stage_abundance_{csv_date}.csv`: Domestic outbreak abundance-effect GAM results.
- `your path/eu-aiv-analysis/aiv_analysis/Domestic_all_birds_single_stage_density_{csv_date}.csv`: Domestic outbreak chicken / duck density-effect GAM results.

Adjustable parameters:

- `--csv-date`: required date tag that matches the AIV result CSV filenames.
- `--main-folder`: project folder containing `aiv_analysis/` and `bird_type_lookup.csv`.

Outputs:

- `your path/eu-aiv-analysis/table/abundance_positive_signi_info.csv`
- `your path/eu-aiv-analysis/plot/Wild_abundance_effect_hist.png`
- `your path/eu-aiv-analysis/plot/Domestic_abundance_effect_hist.png`
- `your path/eu-aiv-analysis/plot/Domestic_Livestock_effect_scatter.png`

Model used:

- No new model is trained here.
- The script summarizes coefficient and p-value outputs from the AIV GAM analysis.

Example:

```powershell
cd "your path/eu-aiv-analysis"
python ".\08-analysis-visualization-cli.py" --csv-date 20260811 --main-folder "."
```

## Troubleshooting

- If a folder is not found, check that all Google Drive data folders are extracted directly under the project root.
- If `05-gpboost-cli.py` cannot find a species, confirm that the scientific name matches the checklist filename in `ebird_filtered_checklist/`.
- If `gpboost` cannot be imported, install it in the active Python environment with `python -m pip install gpboost`.
- If geospatial files fail to load, confirm that all shapefile sidecar files are present, not only the `.shp` file.
- If `07-aiv-outbreak-analysis.R` fails, first check the hard-coded paths near the top of the script.
- If `08-analysis-visualization-cli.py` cannot find CSV files, make sure `--csv-date` matches the date tag generated by `07-aiv-outbreak-analysis.R`.
