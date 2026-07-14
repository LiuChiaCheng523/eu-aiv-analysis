## 撠???

?祉?蝛嗅?獢?璈摮貊??葫 2021 撟渲 2022 撟湧??瘣脩汗瘚?擃◢?芸葆?野蝔桃??征?詨?鞊漲嚗閮?66 蝔柴?蝯?皜祉??蝙??Python Plotly ?脰?鈭?撘閬??蝙?刻暺銝?曈仿??迂嚗?頛??征?雿?撟喳??詨?鞊漲嚗誑?蝛粹??桐???頞典?脩?霈???

?湧??弦蝯??鈭野?瑕?璅∪?嚗摮???蝘駁??穿?憭迤敺?蝘駁?芋???甇??撖?隤歹?銝虫誑?啣???撽??葫蝯???蝥????歇?亦汗瘚??獢??誨蝢拙??扳芋??撽?憭車曈仿?嚗陘?嗆瘞湧野嚗?急??擃陸???砍迤?摮??＊??擃?憸券??

## Project Summary

This research project uses machine learning to predict the spatiotemporal relative abundance of 66 high-risk avian influenza carrier bird species across parts of Europe from 2021 to 2022. The final prediction outputs are visualized through an interactive Python Plotly dashboard. Users can select different bird species to compare average relative abundance across spatial units and examine monthly trend curves within individual spatial units.

The overall results reproduce broad migratory bird movement patterns, with birds moving southward in winter for overwintering and northward in summer. The modeling workflow also corrects for observer bias and uses environmental variables to drive abundance predictions. Subsequent analyses combine known avian influenza outbreak records with generalized additive models to evaluate seasonal risk, showing that many bird species, especially waterbirds, have significantly higher risk during the winter and spring outbreak peaks.

![Plotly Dashboard example](docs/Plotly_Dashboard_example.png)

## 鞈?靘?

?游?憭????皞?撱箇?甇散蝳賣???憸券曈亦車銋?蝛箄?摨西??急?憸券??鞈?摨怒蜓閬???皞??恬?

1. **eBird checklist**
   - Website: https://science.ebird.org/en/use-ebird-data/download-ebird-data-products
   - 雿曈仿?閫皜祉????詨?鞊漲撱箸芋?蜓閬???皞?
   - 鞈??曈亦車閫撖???皜祆????蝵柴?撖?閮?閫皜砍??蝑?雿?
   - ?祉?蝛嗡蝙??2021 ??2022 撟湔??? eBird 曈仿?閫皜祈??脰?璅∪?閮毀??皜研?

2. **eBird Status and Trends**
   - Website: https://science.ebird.org/status-and-trends
   - ?冽撽??祉?蝛嗆?隡啗?銋野憿撠?摨行?蝛箏?撣隅?Ｕ?
   - 閰脰???靘野蝔桀銝????征???撠?摨艾?暹?????蝭?蝑?閮?

3. **FAO EMPRES-i**
   - Website: https://empres-i.apps.fao.org/
   - 雿蝳賣????靘???皞?
   - 雿輻 2021 ??2022 撟湔?瘣脣?摰嗥汗??曈乩?蝳賣??蝝????曈仿?鞊漲????潮◢?芯????

4. **FAO Livestock Systems**
   - Website: https://data.apps.fao.org/
   - 雿摰嗥汗擗?撖漲鞈?靘???
   - 蝝??暾其??桐??Ｙ?擗?撖漲嚗??箇汗瘚?憸券璅∪?銝剔??批霈???券◢?芸?摮?

5. **Google Earth Engine**
   - Website: https://earthengine.google.com/
   - API access: https://developers.google.com/earth-engine/guides/access
   - ERA5-Land: https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR
   - Dynamic World: https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_DYNAMICWORLD_V1
   - ?祉?蝛園? Google Earth Engine ??銵??葫?除?憓???
   - 雿輻 ERA5-Land Daily Aggregated 瘞????靘??啗”皞怠漲???券??◢??瞈漲蝑憓??詻?
   - 雿輻 Dynamic World ?閬?鞈?嚗?敺偌?噙?啜ㄝ?犖?遣?啁??閬?瘥???
   - ?祉?蝛嗡蝙??**Google Earth Engine Python API** 銝?????皜祈???雿輻???隢蒂? GEE API??
   - 銝?鞈???皞?撠??弦蝭?銋?`.shp` 蝛粹???瑼??冽閮??雯?澆??蝔桅?皜祈?瘞?絞閮???

6. **European 100 km ? 100 km grid**
   - QGIS: https://www.qgis.org/
   - ?祉?蝛嗡誑甇散?啣? 100 ?祇? ? 100 ?祇?蝬脫雿銝餉?蝛粹????桐???
   - 蝬脫鞈??臬??QGIS 撠??瘣脣??`.shp` 瑼脰??極??敺遣蝡????箏?閫??摨衣雯?潦??????像皛??楠嚗誑??雿?shapefile 瑼?憭批???
   - ??野憿?皜研憓?摮??啗??振蝳賢?摨西?蝳賣???靘??????遢?雯?潮脰??征????

## ?弦?寞??嗆?

?祉?蝛嗆瘜瑽?鈭蜓閬?畾蛛?鞈??游????????啗??撩憭勗潸??潦野憿撠?摨血遣璅～汗瘚?憸券??嚗誑????閬死??

```text
1. Data Integration and Preprocessing
   鞈??游?????
   ??
   ??? eBird checklist bird observation records
   ??  eBird 曈仿?閫皜祉???
   ??
   ??? ERA5-Land climate variables
   ??  ERA5-Land 瘞?憓?摮?
   ??
   ??? Dynamic World land-cover variables
   ??  Dynamic World ?閬???
   ??
   ??? FAO poultry density data
   ??  FAO 摰嗥汗擗?撖漲鞈?
   ??
   ??? FAO EMPRES-i avian influenza records
   ??  FAO EMPRES-i 蝳賣??蝝??
   ??
   ??? European 100 km ? 100 km grid
       甇散 100 ?祇? ? 100 ?祇?蝬脫
      
    - ?祉?蝛嗅?銝?靘??????偕摨西?蝛粹?撠箏漲????蝯曹??渡???隞?? 100 ?祇? ? 100 ?祇?蝬脫????鞈?銵具?
    - eBird 曈仿?閫皜祉????扯?皜祆??摨扳?嚗???孵??遢?雯?潦?
    - ERA5-Land 瘞???? Dynamic World ?閬?鞈?嚗? Google Earth Engine Python API 閮???隞賬?蝬脫?抒??啣?蝯梯???
    - FAO 摰嗥汗擗?撖漲鞈???FAO EMPRES-i 蝳賣??蝝??靘征??蝵桅?撠?詨?蝬脫嚗蒂靘?隞賢??氬?
    - ?游?敺?鞈?銵其誑瘥??誨銵其???隞?? 蝬脫?雿?雿敺??閬?鋆潦野憿撠?摨血遣璅∟?蝳賣??◢?芸????箇???
        ??

2. Dynamic World Land-cover Imputation
   Dynamic World ?閬?蝻箏仃?潸???
   ??
   ??? Logit transformation for land-cover proportions
   ??  ?閬?瘥? logit 頧?
   ??
   ??? XGBoost regression using climate variables and coordinates
   ??  雿輻瘞??摮?蝬楝摨血漣璅脰? XGBoost 餈湔飛鋆?
   ??
   ??? Validation using RMSE and MAE
       雿輻 RMSE ??MAE ?脰?鋆潛???霅?
        ??

3. Bird Relative Abundance Modeling
   曈仿??詨?鞊漲撱箸芋
   ??
   ??? eBird checklist filtering
   ??  eBird 閫皜祈??祟??
   ??
   ??? Observer bias correction using GLMM
   ??  雿輻 GLMM ?⊥迤閫撖?隤?
   ??
   ??? Spatiotemporal sampling strategy
   ??  ?征?賣見蝑
   ??
   ??? XGBoost Poisson model for each risk bird species
   ??  ???◢?芷野蝔桀遣蝡?XGBoost Poisson 璅∪?
   ??  - 璅∪??孵噩? 59 蝔?ERA5-Land 瘞??摮? 蝔?Dynamic World ?閬????Bird 閫皜砍??嚗誑??撖?隤斗甇??璅???
   ??
   ??? Model evaluation using SRC
   ??  雿輻 SRC嚗pearman Rank Correlation嚗?隡唳芋??皜祈”??
   ??
   ??? Monthly grid-level relative abundance prediction
       ?Ｙ??遢 ? 蝬脫撅斤?銋撠?摨阡?皜?
        ??

4. Avian Influenza Risk Analysis
   蝳賣??◢?芸???
   ??
   ??? Weighted relative abundance by outbreak stage
   ??  靘??畾菔?蝞?甈撠?摨?
   ??
   ??? Integration with poultry and wild bird outbreak records
   ??  蝯?摰嗥汗??曈亦汗瘚??蝝??
   ??
   ??? Adjustment for chicken and duck density
   ??  ?批??暾券?畾?摨?
   ??
   ??? Generalized Additive Model analysis
       雿輻撱?儔?扳芋????◢??
        ??

5. Interactive Visualization
   鈭?撘?閬箏?
   ??
   ??? Plotly Dash dashboard
   ??  Plotly Dash ?銵冽
   ??
   ??? Species-level abundance map
   ??  曈亦車撅斤??詨?鞊漲?啣?
   ??
   ??? Monthly trend curve by spatial unit
       ?蝛粹??桐??隅?Ｘ蝺?
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

| Command | Purpose | Input format / default |
|---|---|---|
| `--mode` | Select which part of the workflow to run. | One of `sampling`, `imputation`, `aggregate`, `all`. Default: `all`. |
| `--eu-shp-path` | Set the EU grid shapefile path. | File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`. |
| `--input-csv` | Set the input CSV with missing land-cover values. | File path. Default: `gee_data/EU_2016_2022_land_cover_and_climate_data_containing_missing_values.csv`. |
| `--output-folder` | Set the project root used for output folders. | Folder path. Default: script folder / project root. |
| `--seed-start` | First seed number. | Integer. Default: `123`. |
| `--seed-end` | Last seed number. | Integer. Default: `222`. |
| `--n-cores` | Number of CPU cores for parallel processing. | Integer. Default: `2`. |
| `--land-cover-types` | Select land-cover types to impute. | Comma-separated values. Default: all 9 types: `bare,built,crops,flooded_vegetation,grass,shrub_and_scrub,snow_and_ice,trees,water`. |
| `--help` | Show help message. | No value. |

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

| Command | Purpose | Input format / default |
|---|---|---|
| `--species` | Run selected species. | Scientific name, comma-separated names, or repeated values. Example: `"Anas crecca"` or `"Anas crecca,Ardea alba"`. |
| `--species-file` | Run species listed in a TXT or CSV file. | File path. CSV can contain `scientific_name`, `birdname`, or use the first column. |
| `--all-species` | Run all species found in `ebird_filtered_checklist/`. | No value. |
| `--list-species` | List available species and exit. | No value. |
| `--eu-shp-path` | Set EU grid shapefile path. | File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`. |
| `--env-folder` | Set ERA5-Land environmental CSV folder. | Folder path. Default: `gee_data/era5_2016_2022/`. |
| `--lc-path` | Set processed land-cover CSV path. | File path. Default: `gee_data/land_cover_2016_2022.csv`. |
| `--checklist-folder` | Set filtered eBird checklist folder. | Folder path. Default: `ebird_filtered_checklist/`. |
| `--output-folder` | Set output root folder. | Folder path. Default: script folder / project root. |
| `--start-year` | Start year for filtering. | Integer. Default: `2021`. |
| `--end-year` | End year for filtering. | Integer. Default: `2022`. |
| `--protocol` | Select eBird protocol(s). | Comma-separated values. Default: `Traveling`. |
| `--obs-quantile-cutoff` | Observation count outlier cutoff. | Numeric quantile. Default: `0.99`. |
| `--observer-quantile-cutoff` | Observer-level cutoff. | Numeric quantile. Default: `0.7`. |
| `--grid-quantile-cutoff` | Grid-level cutoff. | Numeric quantile. Default: `0.5`. |
| `--validation-ratio` | Validation split ratio. | Numeric value between 0 and 1. Default: `0.1`. |
| `--n-iter` | Number of model iterations. | Integer. Default: `3`. |
| `--seed` | Random seed. | Integer. Default: `123`. |
| `--help` | Show help message. | No value. |

Simplest usage:

```powershell
Rscript ".\run_abundance_prediction.R" --species "Anas crecca"
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

| Command | Purpose | Input format / default |
|---|---|---|
| `--species` | Process selected species. | Scientific name, comma-separated names, or repeated values. |
| `--species-file` | Process species listed in a TXT or CSV file. | File path. CSV can contain `scientific_name`, `birdname`, or use the first column. |
| `--all-species` | Process all species available under `abundance_prediction/`. | No value. |
| `--list-species` | List species currently available for post-processing and exit. | No value. |
| `--abundance-root` | Set root folder created by `run_abundance_prediction.R`. | Folder path. Default: `abundance_spatiotemporal_sampling_method/`. |
| `--validation-folder` | Set validation data folder. | Folder path. Default: `validation_data/`. |
| `--checklist-folder` | Set checklist folder. | Folder path. Default: `ebird_filtered_checklist/`. |
| `--bird-type-csv` | Set bird type lookup CSV. | File path. Default: `bird_type_lookup.csv`. |
| `--eu-shp-path` | Set EU grid shapefile path. | File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`. |
| `--output-folder` | Set output root folder. | Folder path. Default: script folder / project root. |
| `--n-iter` | Number of iterations to combine. | Integer. Default: auto-detect. |
| `--mad-k` | MAD filter multiplier. | Numeric value. Default: `1.5`. |
| `--map-years` | Years used for species map output. | Comma-separated years. Default: auto-detect all years. |
| `--map-months` | Months used for species map output. | Comma-separated month numbers. Default: `1,4,7,10`. |
| `--help` | Show help message. | No value. |

Simplest usage:

```powershell
Rscript ".\run_validation_and_prediction_summary.R" --species "Anas crecca"
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

| Command | Purpose | Input format / default |
|---|---|---|
| `--outbreak-type` | Select outbreak type. | Required. One of `Wild` or `Domestic`. |
| `--species` | Analyze selected species. | Scientific name, comma-separated names, or repeated values. |
| `--species-file` | Analyze species listed in a TXT or CSV file. | File path. CSV can contain `scientific_name`, `birdname`, or use the first column. |
| `--all-species` | Analyze all species available in `mad_filter_abundance/`. | No value. |
| `--list-species` | List species currently available for AIV analysis and exit. | No value. |
| `--eu-shp-path` | Set EU grid shapefile path. | File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`. |
| `--chicken-density-path` | Set chicken density CSV path. | File path. Default: `livestock_density_10km/chicken livestock density 10km.csv`. |
| `--duck-density-path` | Set duck density CSV path. | File path. Default: `livestock_density_10km/duck livestock density_2015_10km.csv`. |
| `--aiv-2021-path` | Set 2021 AIV fixed data path. | File path. Default: `aiv_fixed_data/EU aiv fixed data 2021.csv`. |
| `--aiv-2022-path` | Set 2022 AIV fixed data path. | File path. Default: `aiv_fixed_data/EU aiv fixed data 2022.csv`. |
| `--bird-abundance-folder` | Set MAD-filtered abundance folder. | Folder path. Default: `validation_prediction_summary/mad_filter_abundance/`. |
| `--output-folder` | Set output root folder. | Folder path. Default: script folder / project root. |
| `--write-date` | Set output date tag. | Format `YYYYMMDD`. Default: today's date. |
| `--help` | Show help message. | No value. |

Simplest usage:

```powershell
Rscript ".\run_aiv_analysis.R" --outbreak-type Wild --species "Anas crecca"
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

| Command | Purpose | Input format / default |
|---|---|---|
| `--eu-shp-path` | Set EU grid shapefile path. | File path. Default: `EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp`. |
| `--checklist-folder` | Set checklist folder. | Folder path. Default: `ebird_filtered_checklist/`. |
| `--abundance-folder` | Set MAD-filtered abundance folder. | Folder path. Default: `validation_prediction_summary/mad_filter_abundance/`. |
| `--species` | Set initial species shown in the dashboard. | Scientific name. Default: first available species. |
| `--host` | Set Dash host. | Host string. Default: `127.0.0.1`. |
| `--port` | Set Dash port. | Integer. Default: `8050`. |
| `--debug` | Enable Dash debug mode. | No value. Default: off. |

Simplest usage:

```powershell
python ".\Plotly_Dashboard.py"
```

Example usage:

```powershell
python ".\Plotly_Dashboard.py" --species "Vanellus vanellus" --host 127.0.0.1 --port 8051
```
