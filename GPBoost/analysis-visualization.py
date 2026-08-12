# %%
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# %%
bird_types = pd.read_csv("C:/Users/bruce/R/aiv-gpboost/bird_type_lookup.csv")
bird_types["bird_name"] = bird_types["birdname"].str.replace(" ", "_", regex=False)

waterbirds = bird_types[["birdname", "bird_name"]][bird_types["waterbird"] == 1]
waterbirds.head()

# %%
result_folder = Path(r"C:/Users/bruce/R/aiv-gpboost/aiv_analysis_new")
for p in result_folder.iterdir():
    print(p)

# %%
W_abd_result = pd.read_csv(result_folder/"Wild_all_birds_single_stage_abundance_20260811.csv")
W_abd_result.columns = ["bird_name", "C1", "P1", "C2", "P2", "C3", "P3", "num_grid"]
W_abd_result["outbreak_type"] = 'wild'
W_abd_result["bird_type"] = np.where(W_abd_result["bird_name"].isin(waterbirds["bird_name"]), "waterbird", "non_waterbird")


D_abd_result = pd.read_csv(result_folder/"Domestic_all_birds_single_stage_abundance_20260811.csv")
D_abd_result.columns = ["bird_name", "C1", "P1", "C2", "P2", "C3", "P3", "num_grid"]
D_abd_result["outbreak_type"] = 'domestic'
D_abd_result["bird_type"] = np.where(D_abd_result["bird_name"].isin(waterbirds["bird_name"]), "waterbird", "non_waterbird")


W_abd_result["label_1"] = "None"
W_abd_result["label_2"] = "None"
W_abd_result["label_3"] = "None"

D_abd_result["label_1"] = "None"
D_abd_result["label_2"] = "None"
D_abd_result["label_3"] = "None"

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


# %%
print("D1 waterbird:", len(D_abd_result_label[(D_abd_result_label["label_1"]==True) & (D_abd_result_label["bird_type"]=='waterbird')]))
print("D2 waterbird:", len(D_abd_result_label[(D_abd_result_label["label_2"]==True) & (D_abd_result_label["bird_type"]=='waterbird')]))
print("D3 waterbird:", len(D_abd_result_label[(D_abd_result_label["label_3"]==True) & (D_abd_result_label["bird_type"]=='waterbird')]))

print("D1 non-waterbird:", len(D_abd_result_label[(D_abd_result_label["label_1"]==True) & (D_abd_result_label["bird_type"]=='non_waterbird')]))
print("D2 non-waterbird:", len(D_abd_result_label[(D_abd_result_label["label_2"]==True) & (D_abd_result_label["bird_type"]=='non_waterbird')]))
print("D3 non-waterbird:", len(D_abd_result_label[(D_abd_result_label["label_3"]==True) & (D_abd_result_label["bird_type"]=='non_waterbird')]))

print("W1 waterbird:", len(W_abd_result_label[(W_abd_result_label["label_1"]==True) & (W_abd_result_label["bird_type"]=='waterbird')]))
print("W2 waterbird:", len(W_abd_result_label[(W_abd_result_label["label_2"]==True) & (W_abd_result_label["bird_type"]=='waterbird')]))
print("W3 waterbird:", len(W_abd_result_label[(W_abd_result_label["label_3"]==True) & (W_abd_result_label["bird_type"]=='waterbird')]))

print("W1 non-waterbird:", len(W_abd_result_label[(W_abd_result_label["label_1"]==True) & (W_abd_result_label["bird_type"]=='non_waterbird')]))
print("W2 non-waterbird:", len(W_abd_result_label[(W_abd_result_label["label_2"]==True) & (W_abd_result_label["bird_type"]=='non_waterbird')]))
print("W3 non-waterbird:", len(W_abd_result_label[(W_abd_result_label["label_3"]==True) & (W_abd_result_label["bird_type"]=='non_waterbird')]))


df = pd.DataFrame({
    "type": ["Wild", "Wild", "Domestis", "Domestis"],
    "bird_type": ["waterbird", "non-waterbird", "waterbird", "non-waterbird"],
    "period1": [16, 4, 3, 0],
    "period2": [14, 7, 17, 7],
    "period3": [14, 4, 12 ,0],
})

print(df)

# %%
df = pd.DataFrame({
    "type": ["Wild", "Wild", "Domestis", "Domestis"],
    "bird_type": ["waterbird", "non-waterbird", "waterbird", "non-waterbird"],
    "period1": [16, 4, 3, 0],
    "period2": [14, 7, 17, 7],
    "period3": [14, 4, 12 ,0],
})

print(df)

# %%
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

for bird_origin in ["Wild", "Domestis"]:
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
    plt.show()

# %% [markdown]
# ### Livestock Density Visualization

# %%
#ld_result = pd.read_csv(result_folder/"Domestic_all_birds_single_stage_density_20260811.csv")
ld_result = pd.read_csv(result_folder/"Wild_all_birds_single_stage_density_20260811.csv")
ld_result.columns = [
    'bird_name',
    'CK_C1', 'CK_P1', 'DK_C1', 'DK_P1',
    'CK_C2', 'CK_P2', 'DK_C2', 'DK_P2',
    'CK_C3', 'CK_P3', 'DK_C3', 'DK_P3', 
    'num_grid'
]

# %%
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

# %%
print(ld_result1_label["label"].value_counts())
print(ld_result2_label["label"].value_counts())
print(ld_result3_label["label"].value_counts())

# %%
print(ld_result2_label["label"].value_counts())
#c_max = ld_result2_label["CK_C"].max()
#ld_result2_label = ld_result2_label[ld_result2_label["CK_C"] < c_max]

# %%
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
plt.show()

# %%



