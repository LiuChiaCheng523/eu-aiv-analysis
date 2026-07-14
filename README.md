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

1. Bird abundance model training and prediction
2. Validation / ensemble post-processing
3. AIV analysis
4. Land-cover imputation
5. Interactive Plotly dashboard

## Scripts Download

To download the analysis scripts, clone this GitHub repository using Git:

```bash
git clone https://github.com/LiuChiaCheng523/eu-aiv-analysis.git
```

## Data Download

Download the required data from Google Drive:

[Google Drive folder](https://drive.google.com/drive/folders/1YX5561Yos3P4PsK1OfChdaQ7mAzx2Ihm)

Including:

1. aiv_fixed_data
2. ebird_filtered_checklist
3. EU_100km_fishnet_simple_by_distance
4. gee_data
5. livestock_density_10km

After downloading:

1. Extract the zip files.
2. Place the extracted folders under the project root.
3. Keep the internal relative folder structure unchanged.

## Project Structure

You may rename the repository folder after cloning or downloading it.

What must stay unchanged is the **internal relative folder structure**.

In other words:

- The CLI scripts do **not** require the root folder to be named `abundance_r_test`
- The CLI scripts **do** expect sibling folders such as `gee_data/`, `EU_100km_fishnet_simple_by_distance/`, `ebird_filtered_checklist/`, and output folders to remain in the expected relative locations

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

Notes:

- `abundance_spatiotemporal_sampling_method/`, `validation_data/`, `validation_prediction_summary/`, `aiv_analysis/`, and `land_cover_imputation/` are output folders.
- `.venv_plotly_dashboard/` is a local Python virtual environment and usually should not be uploaded to GitHub.

## Recommended Environments

### R

- Recommended version: `R 4.3.3`
- Tested version: `R 4.3.3`

All command examples in this README use `Rscript`.

If `Rscript` is not available in your terminal `PATH`, use the full executable path instead, for example:

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

## Workflow Overview

The main workflow is:

1. `run_abundance_prediction.R`
   Generate species-level abundance predictions.
2. `run_validation_and_prediction_summary.R`
   Combine iterations, compute validation metrics, export MAD-filter abundance, and generate maps.
3. `run_aiv_analysis.R`
   Use MAD-filter abundance outputs for Wild / Domestic AIV analysis.
4. `Plotly_Dashboard.py`
   Explore abundance results interactively.

The land-cover imputation workflow is independent:

1. `run_land_cover_imputation.R`
   Sampling -> imputation -> aggregate

## 1. Abundance Prediction CLI

Script:

- `run_abundance_prediction.R`

Purpose:

- Train abundance models for one, many, or all bird species
- Use checklist, environmental, and land-cover inputs
- Export model outputs under `abundance_spatiotemporal_sampling_method/`

### Help

```powershell
Rscript ".\run_abundance_prediction.R" --help
```

### Common examples

Single species:

```powershell
Rscript ".\run_abundance_prediction.R" --species "Anas crecca"
```

Single species with custom settings:

```powershell
Rscript ".\run_abundance_prediction.R" `
  --species "Anas crecca" `
  --n-iter 5 `
  --obs-quantile-cutoff 0.98 `
  --observer-quantile-cutoff 0.7 `
  --grid-quantile-cutoff 0.5
```

All species:

```powershell
Rscript ".\run_abundance_prediction.R" --all-species
```

List available species:

```powershell
Rscript ".\run_abundance_prediction.R" --list-species
```

### Main arguments

- `--species`
- `--species-file`
- `--all-species`
- `--list-species`
- `--eu-shp-path`
- `--env-folder`
- `--lc-path`
- `--checklist-folder`
- `--output-folder`
- `--start-year`
- `--end-year`
- `--protocol`
- `--obs-quantile-cutoff`
- `--observer-quantile-cutoff`
- `--grid-quantile-cutoff`
- `--validation-ratio`
- `--n-iter`
- `--seed`

## 2. Validation / Prediction Summary CLI

Script:

- `run_validation_and_prediction_summary.R`

Purpose:

- Combine abundance prediction iterations
- Compute validation summary metrics
- Export MAD-filter abundance CSVs
- Export summary plots and species maps
- Export feature-importance land-cover radar plot

This script expects outputs from `run_abundance_prediction.R`.

### Help

```powershell
Rscript ".\run_validation_and_prediction_summary.R" --help
```

### Common examples

Single species:

```powershell
Rscript ".\run_validation_and_prediction_summary.R" --species "Anas crecca"
```

Multiple species:

```powershell
Rscript ".\run_validation_and_prediction_summary.R" `
  --species "Vanellus vanellus,Branta canadensis" `
  --n-iter 50
```

All species:

```powershell
Rscript ".\run_validation_and_prediction_summary.R" --all-species
```

### Main arguments

- `--species`
- `--species-file`
- `--all-species`
- `--list-species`
- `--abundance-root`
- `--validation-folder`
- `--checklist-folder`
- `--bird-type-csv`
- `--eu-shp-path`
- `--output-folder`
- `--n-iter`
- `--mad-k`
- `--map-years`
- `--map-months`

### Main outputs

- `validation_prediction_summary/validation_summary.csv`
- `validation_prediction_summary/mad_filter_abundance/*.csv`
- `validation_prediction_summary/plots/validation_summary_scatter.png`
- `validation_prediction_summary/plots/feature_importance_land_cover_radar.png`
- `validation_prediction_summary/plots/species_maps/...`

## Example Validation SRC Performance

![Plotly Dashboard example](docs/Validation_SRC_Performance_example.png)

## Example Land Cover Feature Importance Radar Chart

![Plotly Dashboard example](docs/Land_Cover_FI_Radar_Chart_example.png)

## 3. AIV Analysis CLI

Script:

- `run_aiv_analysis.R`

Purpose:

- Analyze Wild / Domestic AIV outbreaks using MAD-filter abundance outputs
- Export tables, overview maps, monthly outbreak plots, and weighted-abundance stage maps

This script expects outputs from `run_validation_and_prediction_summary.R`.

### Help

```powershell
Rscript ".\run_aiv_analysis.R" --help
```

### Common examples

Domestic analysis for two species:

```powershell
Rscript ".\run_aiv_analysis.R" `
  --outbreak-type Domestic `
  --species "Vanellus vanellus,Branta canadensis"
```

Wild analysis for all available species:

```powershell
Rscript ".\run_aiv_analysis.R" `
  --outbreak-type Wild `
  --all-species
```

### Main arguments

- `--outbreak-type`
- `--species`
- `--species-file`
- `--all-species`
- `--list-species`
- `--eu-shp-path`
- `--chicken-density-path`
- `--duck-density-path`
- `--aiv-2021-path`
- `--aiv-2022-path`
- `--bird-abundance-folder`
- `--output-folder`
- `--write-date`

### Main outputs

Outputs are organized as:

```text
aiv_analysis/
└─ Domestic_or_Wild/
   └─ YYYYMMDD/
      ├─ *.csv
      ├─ overview_maps/
      └─ weighted_abundance/
```

## Example Gam Model Abundance Effect And P-value Bar Chart

![Plotly Dashboard example](docs/abundance_effect_pvalue_plot.png)

## Example Vanellus vanellus stage2 weighted abundance distribution

![Plotly Dashboard example](docs/stage2_plot.png)

## 4. Land Cover Imputation CLI

Script:

- `run_land_cover_imputation.R`

Purpose:

- Create interpolation validation samples
- Run `xgboost` imputation for land-cover types
- Aggregate multi-seed predictions

### Help

```powershell
Rscript ".\run_land_cover_imputation.R" --help
```

### Common examples

Run all land-cover types:

```powershell
Rscript ".\run_land_cover_imputation.R" `
  --mode all `
  --seed-start 123 `
  --seed-end 124 `
  --n-cores 2
```

Run one land-cover type:

```powershell
Rscript ".\run_land_cover_imputation.R" `
  --mode all `
  --seed-start 123 `
  --seed-end 124 `
  --n-cores 2 `
  --land-cover-types bare
```

Run in three steps:

```powershell
Rscript ".\run_land_cover_imputation.R" --mode sampling --seed-start 123 --seed-end 124 --n-cores 2
Rscript ".\run_land_cover_imputation.R" --mode imputation --seed-start 123 --seed-end 124 --n-cores 2
Rscript ".\run_land_cover_imputation.R" --mode aggregate --seed-start 123 --seed-end 124
```

### Main arguments

- `--mode`
- `--eu-shp-path`
- `--input-csv`
- `--output-folder`
- `--seed-start`
- `--seed-end`
- `--n-cores`
- `--land-cover-types`

### Main outputs

```text
land_cover_imputation/
├─ sampling_data/
├─ ml_prediction_output/
├─ two_method_performance/
├─ two_method_test_output/
└─ final_output/
```

If only one land-cover type is aggregated, the final output is written under a dedicated subfolder, for example:

```text
land_cover_imputation/final_output/bare/EU_2016_2022_land_cover_imputation_by_xgboost_bare.csv
```

## 5. Plotly Dashboard

Script:

- `Plotly_Dashboard.py`

Purpose:

- Explore bird abundance interactively
- Switch species from a dropdown
- Click a grid to inspect time series and summary metrics

### Start the dashboard

```powershell
python ".\Plotly_Dashboard.py"
```

Default URL:

```text
http://127.0.0.1:8050
```

### Common examples

Set initial species:

```powershell
python ".\Plotly_Dashboard.py" --species "Vanellus vanellus"
```

Set host and port:

```powershell
python ".\Plotly_Dashboard.py" --host 127.0.0.1 --port 8051
```

Debug mode:

```powershell
python ".\Plotly_Dashboard.py" --debug
```

### Main arguments

- `--eu-shp-path`
- `--checklist-folder`
- `--abundance-folder`
- `--species`
- `--host`
- `--port`
- `--debug`

### Data source for species dropdown

The dashboard uses the intersection of:

- `ebird_filtered_checklist/*.csv`
- `validation_prediction_summary/mad_filter_abundance/*.csv`
