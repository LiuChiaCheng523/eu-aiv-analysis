import pandas as pd
import geopandas as gpd
import numpy as np
import os
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import argparse

### Define cli parameters
parser = argparse.ArgumentParser()

parser.add_argument(
    "--bird-folder-path",
    type=Path,
    default=Path.cwd()/"gpboost_abundance_outputs",
    help=f"Parent folder containing the abundance outputs created by gpboost. Default: current working directory {Path.cwd()}/gpboost_abundance_outputs."
)

parser.add_argument(
    "--main-folder",
    type=Path,
    default=Path.cwd(),
    help=f"Parent folder containing the ERA5 data, land cover data, and checklist folders. Default: current working directory {Path.cwd()}."
)

args = parser.parse_args()

bird_folder_path = args.bird_folder_path
main_folder = args.main_folder

if not bird_folder_path.exists():
    raise FileNotFoundError(f"{bird_folder_path} does not exist!")

if not main_folder.exists():
    raise FileNotFoundError(f"{main_folder} does not exist!")


### Load data
lc_path = main_folder/"gee_data"/"land_cover_2016_2022.csv"
era5_2021_path = main_folder/"gee_data"/"era5_2016_2022"/"2021_median_combined_result.csv"
era5_2022_path = main_folder/"gee_data"/"era5_2016_2022"/"2022_median_combined_result.csv"
EU_shp_path = main_folder/"EU_100km_fishnet_simple_by_distance"/"EU_100km_fishnet_simple_by_distance.shp"
bird_types_path = main_folder/"bird_type_lookup.csv"

if not lc_path.exists():
    raise FileNotFoundError(f"{lc_path} does not exist!")
lc = pd.read_csv(lc_path)

if not era5_2021_path.exists():
    raise FileNotFoundError(f"{era5_2021_path} does not exist!")
era5_2021 = pd.read_csv(era5_2021_path)

if not era5_2022_path.exists():
    raise FileNotFoundError(f"{era5_2022_path} does not exist!")
era5_2022 = pd.read_csv(era5_2022_path)

if not EU_shp_path.exists():
    raise FileNotFoundError(f"{EU_shp_path} does not exist!")
EU_shp = gpd.read_file(EU_shp_path)

if not bird_types_path.exists():
   raise FileNotFoundError(f"{bird_types_path} does not exist!")
bird_types = pd.read_csv(bird_types_path)


### make output folder direction
output_plot_folder = Path.cwd()/"plot"
os.makedirs(output_plot_folder, exist_ok=True)

output_table_folder = Path.cwd()/"table"
os.makedirs(output_table_folder, exist_ok=True)

### load bird types
bird_types["bird_name"] = bird_types["birdname"].str.replace(" ", "_", regex=False)

waterbirds = bird_types[["birdname", "bird_name"]][bird_types["waterbird"] == 1]
n_waterbirds = bird_types[["birdname", "bird_name"]][bird_types["waterbird"] == 0]
print("number of waterbirds: ", len(waterbirds))
print("number of non-waterbirds: ", len(n_waterbirds))


### get checklist birds
def get_checklist_birds(folder):
    
    checklist_folder_path = folder/"ebird_filtered_checklist"

    if not checklist_folder_path.exists():
        raise FileNotFoundError(f"{checklist_folder_path} does not exists !")

    checklist_bird_names = []
    for p in checklist_folder_path.iterdir():
        bird_name = p.stem.split("_filtered")[0].lstrip("^")
        checklist_bird_names.append(bird_name)
    print("number of bird in checklist folder: ", len(checklist_bird_names))

    return checklist_folder_path, checklist_bird_names

bird_names = []
for p in bird_folder_path.iterdir():
    if p.is_dir():
        bird_names.append(p.name)

checklist_folder_path, checklist_bird_names = get_checklist_birds(main_folder)

for bird_name in bird_names:
    bird_name_check = bird_name.replace("_", " ")

    if bird_name_check not in checklist_bird_names:
        raise ValueError(f"{bird_name_check} is not in checklist_bird_names!")


### Data processing
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


### Save abundance data
for bird_name in bird_names:
    print("bird name: ", bird_name)

    ### Load Data
    abd_path = bird_folder_path/f"{bird_name}"/f"{bird_name}_reconstruction_abundance.csv"
    abd = pd.read_csv(abd_path)

    ### Save Data
    output_folder = Path.cwd()/"all-birds-abd"
    output_folder.mkdir(parents=True, exist_ok=True)

    saved_abd = abd[['Id', 'year_number', 'month_number','relative_abundance']].copy()
    saved_abd.to_csv(output_folder / f"{bird_name}.csv",index=False)

    print(f"{bird_name} data saved successful !")


### save abundance visualization
def plot_save(bird_name, abd, period):

    output_folder = bird_folder_path/f"{bird_name}"/"plot_python"
    output_folder.mkdir(parents=True, exist_ok=True)

    plot_df = pd.merge(EU_shp, abd, on="Id", how="left")

    if period == 1:
        title_name = "Mean Relative Abundance 2021-05 ~ 2021-10"
        save_name ='mean_relative_abundance_2021_05_2021_10'
    elif period == 2:
        title_name = "Mean Relative Abundance 2021-11 ~ 2022-04"
        save_name ='mean_relative_abundance_2021_11_2022_04'
    else:
        title_name = "Mean Relative Abundance 2022-05 ~ 2022-10"
        save_name ='mean_relative_abundance_2022_05_2022_10'

    #print("all grids: ", len(plot_df))
    #print("NOT NA grids: ", len(plot_df[~plot_df["mean_abundance"].isna()]))
    #print("NA grids: ", len(plot_df[plot_df["mean_abundance"].isna()]))

    fig, ax = plt.subplots(figsize=(8,6))

    plot_df.plot(
        column="mean_abundance",
        cmap="viridis",
        legend=True,
        missing_kwds={"color": "lightgrey"},
        ax=ax
    )

    ax.set_title(title_name)
    plt.savefig(f"{output_folder}/{save_name}.png", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"{bird_name} period{period} image saved successful !")


for bird_name in bird_names:
    print("bird name: ", bird_name)

    ### Load Data
    abd_path = bird_folder_path/f"{bird_name}"/f"{bird_name}_reconstruction_abundance.csv"
    abd = pd.read_csv(abd_path)

    #print(abd.groupby(["Id", "year_number", "month_number"]).size().max())

    vmin = abd["relative_abundance"].min()
    vmax = abd["relative_abundance"].max()
    print("minimum of relative abundance:", round(vmin, 3))
    print("maxmum of relative abundance:", round(vmax,3))

    abd[['Id', 'year_number', 'month_number','relative_abundance']]

    ### Calculate Abundance Mean in period1, period2, and period3
    abd["date"] = pd.to_datetime(dict(year=abd["year_number"], month=abd["month_number"], day=15)) 
    abd_period1 = abd[(abd["date"] >= '2021-05-01') & (abd["date"] <= '2021-10-31')]
    abd_period2 = abd[(abd["date"] >= '2021-11-01') & (abd["date"] <= '2022-04-30')]
    abd_period3 = abd[(abd["date"] >= '2022-05-01') & (abd["date"] <= '2022-10-31')]

    abd_period1_mean = (
        abd_period1
        .groupby("Id")["relative_abundance"]
        .mean()
        .rename("mean_abundance")
        .reset_index()
    )

    abd_period2_mean = (
        abd_period2
        .groupby("Id", as_index=False)
        .agg(mean_abundance=("relative_abundance", "mean"))
    )

    abd_period3_mean = (
        abd_period3
        .groupby("Id", as_index=False)
        .agg(mean_abundance=("relative_abundance", "mean"))
    )

    ### Plot & Save Images
    plot_save(bird_name,abd_period1_mean, period=1)
    plot_save(bird_name,abd_period2_mean, period=2)
    plot_save(bird_name, abd_period3_mean, period=3)


### Featrue importance visualization
def get_FI(bird_cols, FI_cols):
    all_lc_FI = pd.DataFrame(FI_cols)
    all_lc_FI.columns = ['feature']
    all_lc_FI["total_rank_point"] = 0

    for bird_name in bird_cols:
        FI_path = bird_folder_path/f"{bird_name}"/f"{bird_name}_feature_importance_gain.csv"
        FI = pd.read_csv(FI_path)

        lc_FI = FI[FI["feature"].isin(lc_cols)].sort_values("importance", ascending=False).copy()
        lc_FI["rank_point"] = np.arange(len(lc_FI), 0, -1)

        all_lc_FI = pd.merge(all_lc_FI, lc_FI[["feature", "rank_point"]], on="feature", how="left")
        all_lc_FI["total_rank_point"] = all_lc_FI["total_rank_point"] + all_lc_FI["rank_point"]
        all_lc_FI = all_lc_FI[["feature", "total_rank_point"]]

    all_lc_FI = all_lc_FI.sort_values("total_rank_point", ascending=False)
    return all_lc_FI


lc_cols = ['built_MEAN', 'flooded_vegetation_MEAN', 'shrub_and_scrub_MEAN','trees_MEAN',
            'crops_MEAN', 'grass_MEAN', 'water_MEAN', 'snow_and_ice_MEAN', 'bare_MEAN']

all_lc_FI = get_FI(bird_names, lc_cols)

waterbird_names = set(waterbirds["bird_name"])
n_waterbird_names = set(n_waterbirds["bird_name"])

waterbirds_cols = []
n_waterbirds_cols = []

for bird_name in bird_names:
    if bird_name in waterbird_names:
        waterbirds_cols.append(bird_name)

    if bird_name in n_waterbird_names:
        n_waterbirds_cols.append(bird_name)

waterbirds_lc_FI = get_FI(waterbirds_cols, lc_cols)
n_waterbirds_lc_FI = get_FI(n_waterbirds_cols, lc_cols)
#waterbirds_lc_FI = get_FI(waterbirds["bird_name"], lc_cols)
#n_waterbirds_lc_FI = get_FI(n_waterbirds["bird_name"], lc_cols)

features = (
    all_lc_FI
    .sort_values("total_rank_point", ascending=False)["feature"]
    .tolist()
)

def get_values(df, features):
    s = df.set_index("feature")["total_rank_point"]
    return [s.get(f, 0) for f in features]

all_values = get_values(all_lc_FI, features)
water_values = get_values(waterbirds_lc_FI, features)
non_water_values = get_values(n_waterbirds_lc_FI, features)

angles = np.linspace(0, 2 * np.pi, len(features), endpoint=False).tolist()
angles += angles[:1]

all_values += all_values[:1]
water_values += water_values[:1]
non_water_values += non_water_values[:1]

max_value = all_lc_FI["total_rank_point"].max()
yticks = [0, max_value / 2, max_value]

fig, ax = plt.subplots(figsize=(8, 6), subplot_kw=dict(polar=True))

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(1)

# Turn off the default circular grid
ax.grid(False)
ax.spines["polar"].set_visible(False)

# Draw a polygon mesh yourself
for y in yticks[1:]:
    ax.plot(
        angles,
        [y] * len(angles),
        color="lightgrey",
        linewidth=1,
        linestyle="-"
    )

# Manually draw the radial lines for each feature.
for angle in angles[:-1]:
    ax.plot(
        [angle, angle],
        [0, max_value],
        color="lightgrey",
        linewidth=1
    )

ax.plot(angles, water_values, color="blue", linewidth=2, label="water")
ax.fill(angles, water_values, color="blue", alpha=0.08)

ax.plot(angles, non_water_values, color="green", linewidth=2, label="non_water")
ax.fill(angles, non_water_values, color="green", alpha=0.08)

ax.plot(angles, all_values, color="red", linewidth=2, label="all")
ax.fill(angles, all_values, color="red", alpha=0.12)

ax.set_xticks(angles[:-1])
ax.set_xticklabels([])

label_radius = max_value * 1.10

for angle, feature in zip(angles[:-1], features):
    ax.text(
        angle,
        label_radius,
        feature,
        fontsize=10,
        ha="center",
        va="center",
        zorder=20,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1)
    )

ax.set_ylim(0, max_value)
ax.set_yticks(yticks)
ax.set_yticklabels([f"{int(x)}" for x in yticks], fontsize=10)

ax.set_title(
    "Land Cover Feature Importance Ranking",
    fontsize=14,
    fontweight="bold",
    pad=25
)

ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.1))

plt.tight_layout()
plt.savefig(f"{output_plot_folder}/FI_radar.png", dpi=400, bbox_inches="tight")
plt.close()


### Model fitting performance visualization
all_SRC = []

for bird_name in bird_names:
    SRC_path = bird_folder_path/f"{bird_name}"/f"{bird_name}_model_performance.csv"
    SRC = pd.read_csv(SRC_path)
    all_SRC.append(SRC)

all_SRC = pd.concat(all_SRC, axis=0, ignore_index=True)

all_SRC["bird_type"] =np.where(all_SRC["birdname"].isin(waterbirds["birdname"]), "waterbird", "non_waterbird")
all_SRC = all_SRC[['birdname', 'bird_type', 'n', 'Spearman SRC']].copy()
all_SRC.columns = ['birdname', 'bird_type', 'sample_size', 'SRC']
print(all_SRC.head())


for bird_type, color in {
    "waterbird": "blue",
    "non_waterbird": "green"
}.items():
    df = all_SRC[all_SRC["bird_type"] == bird_type]

    plt.scatter(
        df["sample_size"],
        df["SRC"],
        color=color,
        label=bird_type,
        alpha=0.8
    )

plt.grid(True, linestyle="--", alpha=0.4)
plt.xlabel("Sample Size")
plt.ylabel("Spearman Rank Correlation")
plt.title("Model Fit Evaluation")
plt.legend()
plt.savefig(f"{output_plot_folder}/SRC_scatter.png", dpi=400, bbox_inches="tight")
plt.close()


### Save checklist and model information
def checklist_process(checklist):

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

    observer_info = pd.DataFrame(df.groupby("observer_id").size().reset_index())
    observer_info.columns = ["observer_id", "count"]
    orc_q70 = observer_info["count"].quantile(0.7)
    observer_list = observer_info["observer_id"][observer_info["count"] > orc_q70].copy()
    df = df[df["observer_id"].isin(observer_list)]

    id_info = pd.DataFrame(df.groupby("Id").size().reset_index())
    id_info.columns = ["Id", "count"]
    idc_q50 = id_info["count"].quantile(0.5)
    id_list = id_info["Id"][id_info["count"] > idc_q50]
    df = df[df["Id"].isin(id_list)]

    return df

all_ob_info_rows = []
for bird_name in bird_names:
    birdname = bird_name.replace("_", " ")
    print("bird name: ", birdname)

    checklist_path1 = checklist_folder_path/f"^{birdname}_filtered_2019to2022.csv"
    checklist_path2 = checklist_folder_path/f"^{birdname}_filtered2019to2022.csv"
    if checklist_path1.exists():
        checklist = pd.read_csv(checklist_path1)
    elif checklist_path2.exists():
        checklist = pd.read_csv(checklist_path2)
    else:
        raise FileNotFoundError(f"Cannot find checklist file for {birdname}")
        
    ob_re_path = bird_folder_path/f"{bird_name}"/f"{bird_name}_random_effect.csv"
    ob_re_cov_path = bird_folder_path/f"{bird_name}"/f"{bird_name}_random_effect_cov_pars.csv"
    ob_re = pd.read_csv(ob_re_path)
    ob_re_cov = pd.read_csv(ob_re_cov_path)

    df = checklist_process(checklist)

    print(f"number of observer in checklist: {len(df["observer_id"].unique())}")
    print(f"number of observer in random effect table: {len(ob_re["observer_id"].unique())}")

    ob_count = df["observer_id"].value_counts().reset_index()
    ob_re_unique = ob_re.drop_duplicates(subset="observer_id")
    ob_count = pd.merge(ob_count, ob_re_unique, on="observer_id", how="left")
    corr = ob_count[["count", "observer_re"]].corr()
    #print(corr.loc["observer_re", "count"])

    all_ob_info_rows.append({
        "birdname": birdname,
        "sample_size": len(df),
        "re_var": ob_re_cov["Group_1"].iloc[0],
        "cl_std": df["observer_id"].value_counts().std()
    })

all_ob_info = pd.DataFrame(all_ob_info_rows)
all_info = pd.merge(all_SRC[["birdname", "bird_type", "SRC"]], all_ob_info, on="birdname")
print(all_info.head())
all_info.to_csv(f"{output_table_folder}/checklist_and_model_info.csv", index=False)
print(f"checklist and model information are saved sucessfully to {output_table_folder}/checklist_and_model_info.csv")