import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

### Define cli parameters
parser = argparse.ArgumentParser()

parser.add_argument(
    "--csv-date",
    required=True,
    help=f"The date preceding .csv."
)

parser.add_argument(
    "--main-folder",
    type=Path,
    default=Path.cwd(),
    help=f"Parent folder containing aiv analysis results. Default: current working directory {Path.cwd()}."
)

args = parser.parse_args()

csv_date = args.csv_date
main_folder = args.main_folder

if not main_folder.exists():
    raise FileNotFoundError(f"{main_folder} does not exist!")

### make output folder direction
output_plot_folder = Path.cwd()/"plot"
os.makedirs(output_plot_folder, exist_ok=True)

output_table_folder = Path.cwd()/"table"
os.makedirs(output_table_folder, exist_ok=True)

### Load data
bird_types_path = main_folder/"bird_type_lookup.csv"
bird_types = pd.read_csv(bird_types_path)
bird_types["bird_name"] = bird_types["birdname"].str.replace(" ", "_", regex=False)

waterbirds = bird_types[["birdname", "bird_name"]][bird_types["waterbird"] == 1]

result_folder_path = main_folder/"aiv_analysis"
W_abd_result_path = result_folder_path/f"Wild_all_birds_single_stage_abundance_{csv_date}.csv"
D_abd_result_path = result_folder_path/f"Domestic_all_birds_single_stage_abundance_{csv_date}.csv"

if not W_abd_result_path.exists():
    raise FileNotFoundError(f"{W_abd_result_path} does not exist!")

if not D_abd_result_path.exists():
    raise FileNotFoundError(f"{D_abd_result_path} does not exist!")

W_abd_result = pd.read_csv(W_abd_result_path)
W_abd_result.columns = ["bird_name", "C1", "P1", "C2", "P2", "C3", "P3", "num_grid"]
W_abd_result["outbreak_type"] = 'wild'
W_abd_result["bird_type"] = np.where(W_abd_result["bird_name"].isin(waterbirds["bird_name"]), "waterbird", "non_waterbird")

D_abd_result = pd.read_csv(D_abd_result_path)
D_abd_result.columns = ["bird_name", "C1", "P1", "C2", "P2", "C3", "P3", "num_grid"]
D_abd_result["outbreak_type"] = 'domestic'
D_abd_result["bird_type"] = np.where(D_abd_result["bird_name"].isin(waterbirds["bird_name"]), "waterbird", "non_waterbird")


for col in ["label_1", "label_2", "label_3"]:
    W_abd_result[col] = False
    D_abd_result[col] = False

def cat_label(df):
    for i in range(len(df)):
        if df["C1"].iloc[i] > 0 and df["P1"].iloc[i] < 0.05:
            df.loc[df.index[i], "label_1"] = True
        else:
            df.loc[df.index[i], "label_1"] = False

        if df["C2"].iloc[i] > 0 and df["P2"].iloc[i] < 0.05:
            df.loc[df.index[i], "label_2"] = True
        else:
            df.loc[df.index[i], "label_2"] = False

        if df["C3"].iloc[i] > 0 and df["P3"].iloc[i] < 0.05:
            df.loc[df.index[i], "label_3"] = True
        else:
            df.loc[df.index[i], "label_3"] = False

    return df

W_abd_result_label = cat_label(W_abd_result)
D_abd_result_label = cat_label(D_abd_result)
print(W_abd_result_label.head())
print(D_abd_result_label.head())



def count_labels(df, type_name):
    rows = []

    for bird_type in ["waterbird", "non_waterbird"]:
        row = {
            "type": type_name,
            "bird_type": bird_type,
        }

        for i in [1, 2, 3]:
            row[f"period{i}"] = len(
                df[
                    (df[f"label_{i}"] == True) &
                    (df["bird_type"] == bird_type)
                ]
            )

        rows.append(row)

    return rows


rows = []
rows.extend(count_labels(W_abd_result_label, "Wild"))
rows.extend(count_labels(D_abd_result_label, "Domestic"))

df = pd.DataFrame(rows)

df["bird_type"] = df["bird_type"].replace({
    "non_waterbird": "non-waterbird"
})

df.to_csv(f"{output_table_folder}/abundance_positive_signi_info.csv", index=False)
print(f"positive and significant effct of abundance information are saved sucessfully to {output_table_folder}/abundance_positive_signi_info.csv")

# wide -> long
df_long = df.melt(
    id_vars=["type", "bird_type"],
    value_vars=["period1", "period2", "period3"],
    var_name="period",
    value_name="count"
)

period_labels = {
    "period1": "2021/5–2021/10",
    "period2": "2021/11–2022/4",
    "period3": "2022/5–2022/10",
}

df_long["period"] = df_long["period"].map(period_labels)

df_long["bird_type"] = pd.Categorical(
    df_long["bird_type"],
    categories=["non-waterbird", "waterbird"],
    ordered=True
)

for bird_origin in ["Wild", "Domestic"]:
    sub = df_long[df_long["type"] == bird_origin]

    plt.figure(figsize=(8, 6))

    ax = sns.barplot(
        data=sub,
        x="period",
        y="count",
        hue="bird_type",
        width=0.55,
        palette={
            "non-waterbird": "#2ca02c",
            "waterbird": "#1f77b4",
        }
    )

    ax.set_title(
        f"{bird_origin}: positive & significant effects",
        fontsize=13,
        fontweight="bold"
    )
    ax.set_xlabel("period")
    ax.set_ylabel("count of significant positive effects")

    ax.legend(title="Group", frameon=False)

    plt.tight_layout()
    plt.savefig(f"{output_plot_folder}/{bird_origin}_abundance_effect_hist.png", dpi=400, bbox_inches="tight")
    plt.close()
    print(f"{bird_origin} abundance effect histogram is saved sucessfully to {output_plot_folder}/{bird_origin}_abundance_effect_hist.png")


### Livestock Density Visualization (Domestic)
ld_result_path = result_folder_path/f"Domestic_all_birds_single_stage_density_{csv_date}.csv"

if not ld_result_path.exists():
    raise FileNotFoundError(f"{ld_result_path} does not exist!")

ld_result = pd.read_csv(ld_result_path)

ld_result.columns = [
    'bird_name',
    'CK_C1', 'CK_P1', 'DK_C1', 'DK_P1',
    'CK_C2', 'CK_P2', 'DK_C2', 'DK_P2',
    'CK_C3', 'CK_P3', 'DK_C3', 'DK_P3', 
    'num_grid'
]


ld_result1 = ld_result[['bird_name', 'CK_C1', 'CK_P1', 'DK_C1', 'DK_P1']].copy()
ld_result1.columns = ['bird_name', 'CK_C', 'CK_P', 'DK_C', 'DK_P']
ld_result1['label'] = 'None'

ld_result2 = ld_result[['bird_name', 'CK_C2', 'CK_P2', 'DK_C2', 'DK_P2']].copy()
ld_result2.columns = ['bird_name', 'CK_C', 'CK_P', 'DK_C', 'DK_P']
ld_result2['label'] = 'None'

ld_result3 = ld_result[['bird_name', 'CK_C3', 'CK_P3', 'DK_C3', 'DK_P3']].copy()
ld_result3.columns = ['bird_name', 'CK_C', 'CK_P', 'DK_C', 'DK_P']
ld_result3['label'] = 'None'


def cat_LD_label(df):
    for i in range(len(df)):
        if df["DK_C"].iloc[i] > 0 and df["DK_P"].iloc[i] < 0.05:
            if df["CK_C"].iloc[i] > 0 and df["CK_P"].iloc[i] < 0.05:
                df.loc[df.index[i], "label"] = "all"
            else:
                df.loc[df.index[i], "label"] = "only duck"
        elif df["CK_C"].iloc[i] > 0 and df["CK_P"].iloc[i] < 0.05:
            df.loc[df.index[i], "label"] = "only chicken"
        else:
            df.loc[df.index[i], "label"] = "Non"

    return df

ld_result1_label = cat_LD_label(ld_result1)
ld_result2_label = cat_LD_label(ld_result2)
ld_result3_label = cat_LD_label(ld_result3)


#print(ld_result1_label["label"].value_counts())
#print(ld_result2_label["label"].value_counts())
#print(ld_result3_label["label"].value_counts())

c_max = ld_result2_label["CK_C"].max()
ld_result2_label = ld_result2_label[ld_result2_label["CK_C"] < c_max]


plt.figure(figsize=(6, 5))

ax = sns.scatterplot(
    data=ld_result2_label,
    x="CK_C",
    y="DK_C",
    hue="label",
    s=70,
    alpha=0.8,
    palette={
        "all": "#1f77b4",
        "only duck": "#ff7f0e",
        "only chicken": "#2ca02c",
        "Non": "#7f7f7f",
    }
)

ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
ax.axvline(0, color="black", linewidth=0.8, linestyle="--")

ax.set_xlabel("Chicken effect")
ax.set_ylabel("Duck effect")
ax.set_title("Livestock Density Effect & Significance")

ax.legend(title="label", frameon=False)

plt.tight_layout()
plt.savefig(f"{output_plot_folder}/Domestic_Livestock_effect_scatter.png", dpi=400, bbox_inches="tight")
plt.close()
print(f"Domestic Livestock effect scatter plot is saved sucessfully to {output_plot_folder}/Domestic_Livestock_effect_scatter.png")

