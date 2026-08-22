import os
import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path
import argparse

from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from scipy.stats import spearmanr
import gpboost as gpb

### 1. Load Data
## Define function
def get_checklist_birds(folder):
    
    bird_folder_path = folder/"ebird_filtered_checklist"

    if not bird_folder_path.exists():
        raise FileNotFoundError(f"{bird_folder_path} does not exists !")

    bird_names = []
    for p in bird_folder_path.iterdir():
        bird_name = p.stem.split("_filtered")[0].lstrip("^")
        bird_names.append(bird_name)
    print("number of bird in checklist folder: ", len(bird_names))

    return bird_folder_path, bird_names

## Define cli parameters
parser = argparse.ArgumentParser()

parser.add_argument(
    "--birds",
    nargs="+",
    default=["all"],
    help='Scientific name of bird for example "Anas crecca". Default is "all" to run all birds.'
)

parser.add_argument(
    "--path-folder",
    type=Path,
    default=Path.cwd(),
    help=f"Parent folder containing the ERA5 data, land cover data, and checklist folders. Default: current working directory {Path.cwd()}."
)

parser.add_argument(
    "--output-folder",
    type=Path,
    default=Path.cwd(),
    help=f"Default: current working directory {Path.cwd()}."
)

args = parser.parse_args()

path_folder = args.path_folder
output_folder = args.output_folder

## Check whether the bird names selected by the user all exist in the checklist folder.
bird_folder_path, bird_names = get_checklist_birds(path_folder)

if args.birds == ["all"]:
    choose_bird_names = bird_names
else:
    choose_bird_names = args.birds

for bird_name in choose_bird_names:
    if bird_name not in bird_names:
        raise ValueError(f"{bird_name} is not in checklist bird names!")


## If land cover data and ERA5 data all exist, load them.
lc_path = path_folder/"gee_data"/"land_cover_2016_2022.csv"
era5_2021_path = path_folder/"gee_data"/"era5_2016_2022"/"2021_median_combined_result.csv"
era5_2022_path = path_folder/"gee_data"/"era5_2016_2022"/"2022_median_combined_result.csv"

if not lc_path.exists():
    raise FileNotFoundError(f"{lc_path} does not exist!")
lc = pd.read_csv(lc_path)

if not era5_2021_path.exists():
    raise FileNotFoundError(f"{era5_2021_path} does not exist!")
era5_2021 = pd.read_csv(era5_2021_path)

if not era5_2022_path.exists():
    raise FileNotFoundError(f"{era5_2022_path} does not exist!")
era5_2022 = pd.read_csv(era5_2022_path)


### 2. Data processing
era5_2021 = era5_2021.rename(columns={"Month": "month_number"})
era5_2022 = era5_2022.rename(columns={"Month": "month_number"})
era5_2021["month_number"] = pd.to_datetime(era5_2021["month_number"]).dt.month
era5_2022["month_number"] = pd.to_datetime(era5_2022["month_number"]).dt.month
era5_2022["year_number"] = 2022
era5_2021["year_number"] = 2021

era5_all = pd.concat([era5_2021, era5_2022], ignore_index=True)
print(f"length of era5_2021: {len(era5_2021)}")
print(f"length of era5_2022: {len(era5_2022)}")
print(f"length of era5_all: {len(era5_all)}")

env = pd.merge(lc, era5_all, on=["Id","month_number", "year_number"], how="inner")
print(f"length of lc: {len(lc)}")
print(f"length of env: {len(env)}")


### 3. Model fitting
for bird_name in choose_bird_names:
  ### load checklist data
  print("bird name: ", bird_name)
  path1 = f"{bird_folder_path}/^{bird_name}_filtered_2019to2022.csv"
  path2 = f"{bird_folder_path}/^{bird_name}_filtered2019to2022.csv"

  if os.path.exists(path1):
      checklist_path = path1
  elif os.path.exists(path2):
      checklist_path = path2
  else:
      raise FileNotFoundError(f"Cannot find checklist file for {bird_name}")

  checklist = pd.read_csv(checklist_path)

  checklist['observation_date'] = pd.to_datetime(checklist['observation_date'])
  checklist['year_number'] = checklist['observation_date'].dt.year
  checklist['month_number'] = checklist['observation_date'].dt.month

  cl_include_cols = ['category','common_name','observation_date','all_species_reported','group_identifier','time_started']
  cl_cols = [col for col in checklist.columns if col not in cl_include_cols]
  checklist = checklist[cl_cols].copy()

  checklist["observation_count"] = pd.to_numeric(checklist["observation_count"].replace("X", 0))
  checklist["locality_type"] = checklist["locality_type"].map({"P": 0, "H": 1}).astype(int)

  protocol_cols = ["Traveling", "Stationary"]
  checklist = checklist[checklist["protocol_type"].isin(protocol_cols)].copy()
  checklist["protocol_type"] = checklist["protocol_type"].map({"Stationary": 0, "Traveling": 1}).astype(int)
  checklist = checklist[checklist["protocol_type"] == 1].copy()

  print("lehgth of checklist : ", len(checklist))
  print("lehgth of checklist in 2019 : ", len(checklist[checklist["year_number"]==2019]))
  print("lehgth of checklist in 2020 : ", len(checklist[checklist["year_number"]==2020]))
  print("lehgth of checklist in 2021 : ", len(checklist[checklist["year_number"]==2021]))
  print("lehgth of checklist in 2022 : ", len(checklist[checklist["year_number"]==2022]))

  ### merge checklist and env data
  df = pd.merge(checklist, env, on=["Id","month_number", "year_number"], how="inner")
  print("Before dropna:", len(df))
  df = df.dropna().reset_index(drop=True)
  print("After dropna:", len(df))

  ### filter data
  oc_q99 = df["observation_count"].quantile(0.99)
  df = df[df["observation_count"] <= oc_q99].copy()
  print("99% cutoff of observation count: ", oc_q99)
  print("rows after cut: ", len(df))

  observer_info = pd.DataFrame(df.groupby("observer_id").size().reset_index())
  observer_info.columns = ["observer_id", "count"]
  orc_q70 = observer_info["count"].quantile(0.7)
  observer_list = observer_info["observer_id"][observer_info["count"] > orc_q70].copy()
  df = df[df["observer_id"].isin(observer_list)]
  print(f"70% cutoff of observer count ;{orc_q70}")
  print("number of observers: ", len(observer_list))
  print("rows after cut: ", len(df))

  id_info = pd.DataFrame(df.groupby("Id").size().reset_index())
  id_info.columns = ["Id", "count"]
  idc_q50 = id_info["count"].quantile(0.5)
  id_list = id_info["Id"][id_info["count"] > idc_q50]
  df = df[df["Id"].isin(id_list)]
  print(f"50% cutoff of Id count ;{idc_q50}")
  print("number of Ids: ", len(id_list))
  print("rows after cut: ", len(df))

  df['sample_id'] = range(1,len(df)+1)
  print("all samples:", len(df))

  ### selects model features
  model_include_cols = [
      "scientific_name", "observer_id", "Id", "month_number", "year_number", "sample_id",
      "observation_count", "protocol_type", "number_observers", "longitude", "latitude"
  ]
  model_cols = [col for col in df.columns if col not in model_include_cols]

  x_train = df[model_cols].copy()
  y_train = df["observation_count"].copy()

  group_train = df[["observer_id"]].copy()

  assert len(x_train) == len(y_train) == len(group_train)

  print("Train:", x_train.shape)
  print("Train observers:", group_train["observer_id"].nunique())
  print("X train NA:", x_train.isna().sum().sum())
  print("Y train NA:", pd.isna(y_train).sum())

  X_train_np = x_train.to_numpy(dtype=np.float64)
  y_train_np = np.asarray(y_train, dtype=np.float64)
  group_train_np = group_train.astype(str).to_numpy()

  gp_model = gpb.GPModel(
    group_data=group_train_np,
    likelihood="negative_binomial" #"gaussian" "negative_binomial"
  )

  train_data = gpb.Dataset(
      data=X_train_np,
      label=y_train_np,
      feature_name=list(x_train.columns)
  )

  params = {
      "learning_rate": 0.03,
      "num_leaves": 31,
      "max_depth": -1,
      "min_data_in_leaf": 100,

      "feature_fraction": 0.7,

      "lambda_l1": 0.0,
      "lambda_l2": 1.0,

      "bagging_freq": 0,

      "verbose": 0,
      "seed": 42
  }

  ### fitting model
  model = gpb.train(
      params=params,
      train_set=train_data,
      gp_model=gp_model,
      num_boost_round=700
  )

  cov_pars = gp_model.get_cov_pars()
  print(cov_pars)
  re_train = model.predict_training_data_random_effects()

  re_check = pd.DataFrame({
    "observer_id": df["observer_id"].values,
    "observer_re": re_train["Group_1"].values,
  })

  print(re_check.groupby("observer_id")["observer_re"].nunique().value_counts())
  print(re_check.head())

  ### model performance
  pred_train = model.predict(
    data=X_train_np,
    group_data_pred=group_train_np,
    predict_response=True
  )

  y_pred_train = np.asarray(
      pred_train["response_mean"]
  ).reshape(-1)

  src_train, _ = spearmanr(
      y_train_np,
      y_pred_train
  )

  mae_train = mean_absolute_error(
      y_train_np,
      y_pred_train
  )

  print("Model SRC:", src_train)
  print("Model MAE:", mae_train)


  ### feature importance
  importance_df = pd.DataFrame({
    "feature": model.feature_name(),
    "importance": model.feature_importance(
        importance_type="gain"
    )
  })

  importance_df["importance"] = (
      importance_df["importance"] /
      importance_df["importance"].sum()
  )

  importance_df = importance_df.sort_values(
      "importance",
      ascending=False
  ).reset_index(drop=True)

  print(importance_df)


  ### reconstruction Abundance
  valid_ids = df["Id"].unique()
  recon_df = env[env["Id"].isin(valid_ids)].copy()
  recon_df = recon_df[recon_df["year_number"].isin([2021, 2022])].copy()

  print("df unique Id:", df["Id"].nunique())
  print("recon unique Id:", recon_df["Id"].nunique())
  print("recon rows:", len(recon_df))

  reference_duration = df["duration_mins"].median()
  reference_effort = df["effort_km"].median()
  reference_locality = df["locality_type"].mode()[0]

  recon_df["duration_mins"] = reference_duration
  recon_df["effort_km"] = reference_effort
  recon_df["locality_type"] = reference_locality

  missing_cols = [col for col in model_cols if col not in recon_df.columns]
  if len(missing_cols) > 0:
      raise ValueError(f"recon_df is missing model features: {missing_cols}")

  X_recon = recon_df[model_cols].copy()

  print("Number of reconstruction samples:", len(X_recon))
  print("Number of features:", X_recon.shape[1])
  print("Feature order same as training:", list(X_recon.columns) == list(x_train.columns))

  if list(X_recon.columns) != list(x_train.columns):
      raise ValueError("X_recon feature columns/order does not match x_train")

  print("NA count:", X_recon.isna().sum().sum())

  X_recon_np = X_recon.to_numpy(dtype=np.float64)

  pred_recon = model.predict(
    data=X_recon_np,
    ignore_gp_model=True
  )
  recon_df["birdname"] = bird_name
  recon_df["relative_abundance"] = pred_recon
  recon_df["reference_duration_mins"] = reference_duration
  recon_df["reference_effort_km"] = reference_effort
  recon_df["reference_locality_type"] = reference_locality

  print("===== Reconstruction summary =====")
  print(recon_df["relative_abundance"].describe())


  ### output result
  output_root = output_folder/"gpboost_abundance_outputs"
  os.makedirs(output_root, exist_ok=True)

  bird_slug = re.sub(
      r"[^A-Za-z0-9._-]+",
      "_",
      bird_name
  ).strip("_")

  bird_output_dir = os.path.join(output_root, bird_slug)
  os.makedirs(bird_output_dir, exist_ok=True)

  reconstruction_path = os.path.join(
      bird_output_dir,
      f"{bird_slug}_reconstruction_abundance.csv"
  )

  feature_importance_path = os.path.join(
      bird_output_dir,
      f"{bird_slug}_feature_importance_gain.csv"
  )

  performance_path = os.path.join(
      bird_output_dir,
      f"{bird_slug}_model_performance.csv"
  )

  random_effect_path = os.path.join(
      bird_output_dir,
      f"{bird_slug}_random_effect.csv"
  )

  cov_pars_path = os.path.join(
      bird_output_dir,
      f"{bird_slug}_random_effect_cov_pars.csv"
  )

  recon_df.to_csv(reconstruction_path, index=False)
  importance_df.to_csv(feature_importance_path, index=False)
  re_check.to_csv(random_effect_path, index=False)

  performance_df = pd.DataFrame([{
      "birdname": bird_name,
      "n": len(df),
      "Spearman SRC": src_train,
      "MAE": mae_train,
  }])

  performance_df.to_csv(performance_path, index=False)

  cov_pars_df = pd.DataFrame(cov_pars)
  cov_pars_df.to_csv(cov_pars_path, index=False)

  print("Saved CSV files to:", bird_output_dir) 