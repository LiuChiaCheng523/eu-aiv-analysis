import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns

from dash import Dash, dcc, html, Input, Output
import plotly.express as px

EU_shp_path = "D:/aiv_project_2025/EU_100km_fishnet_simple_by_distance/EU_100km_fishnet_simple_by_distance.shp"
Vanellus_checkist_path = "D:/aiv_project_2025/ebird filtered checklist2/^Vanellus vanellus_filtered_2019to2022.csv"
Vanellus_abundance_path = "D:/aiv_project_table/abundance_spatiotemporal_sampling_method/MAD_filtered_abundance/Anas crecca.csv"

EU_shp = gpd.read_file(EU_shp_path)
checklist_df = pd.read_csv(Vanellus_checkist_path)
abd_df = pd.read_csv(Vanellus_abundance_path)


# 建立 time_index
abd_df["time_index"] = abd_df["year_number"] * 12 + abd_df["month_number"]

start_year, start_month = 2021, 1
end_year, end_month = 2022, 12

start = start_year * 12 + start_month
end = end_year * 12 + end_month

abd_sub = abd_df[
    (abd_df["time_index"] >= start) &
    (abd_df["time_index"] <= end)
]

# mean abundance
mean_df = abd_sub.groupby("Id")["abundance_filtered_mean"].mean().reset_index()
mean_df.rename(columns={"abundance_filtered_mean": "mean_abundance"}, inplace=True)

# trend
def calc_slope(group):
    x = group["time_index"].values
    y = group["abundance_filtered_mean"].values
    
    if len(x) < 2:
        return np.nan
    
    return np.polyfit(x, y, 1)[0]

trend_df = abd_sub.groupby("Id").apply(calc_slope).reset_index(name="trend")

# max - min range
diff_df = abd_sub.groupby("Id")["abundance_filtered_mean"].agg(
    abd_max="max",
    abd_min="min"
).reset_index()

diff_df["abd_range"] = diff_df["abd_max"] - diff_df["abd_min"]

# checklist sample size
count_df = checklist_df.groupby("Id").size().reset_index(name="n_samples")

# merge
result = EU_shp.merge(mean_df, on="Id", how="left")
result = result.merge(trend_df, on="Id", how="left")
result = result.merge(diff_df[["Id", "abd_range"]], on="Id", how="left")
result = result.merge(count_df, on="Id", how="left")
result.columns

# feature engineering
result["log_samples"] = np.log1p(result["n_samples"])
result["log_abd_range"] = np.log1p(result["abd_range"])
result["trend_clip"] = result["trend"].clip(-5, 5)

# centroid（給 dashboard 用）
result["lon"] = result.geometry.centroid.x
result["lat"] = result.geometry.centroid.y

# 檢查
print(result.head())
print(result.columns)


result["log_samples"] = result["log_samples"].fillna(0)
result["mean_abundance"] = result["mean_abundance"].fillna(0)

map_fig = px.scatter_map(
    result,
    lat="lat",
    lon="lon",
    color="mean_abundance",
    size="log_samples",
    color_continuous_scale="viridis_r",
    size_max=20,
    zoom=3,
    height=600,
    hover_data=["Id", "mean_abundance", "log_samples"],
    custom_data=["Id"]   
)

map_fig.update_layout(mapbox_style="carto-positron")
#map_fig.show()

# Dash App
app = Dash(__name__)

app.layout = html.Div([
    
    html.H2("Spatiotemporal Bird Analysis Dashboard"),
    
    html.Div([
        
        # 左：地圖
        html.Div([
            dcc.Graph(id="map", figure=map_fig)
        ], style={"width": "50%", "display": "inline-block"}),
        
        # 右：時間序列 + 統計
        html.Div([
            dcc.Graph(id="trend_plot"),
            html.Div(id="stats_output")
        ], style={"width": "50%", "display": "inline-block"})
        
    ])
])



@app.callback(
    Output("trend_plot", "figure"),
    Output("stats_output", "children"),
    Input("map", "clickData")
)
def update_output(clickData):

    print(clickData)

    if clickData is None:
        return {}, "Click a grid"

    # ===== 取得點擊位置 =====
    point = clickData["points"][0]
    lon = point["lon"]
    lat = point["lat"]

    # ===== 找最近 grid =====
    idx = ((result["lon"] - lon)**2 + (result["lat"] - lat)**2).idxmin()
    selected = result.loc[idx]
    selected_id = selected["Id"]

    # ===== 取 time series =====
    ts = abd_df[abd_df["Id"] == selected_id].copy()

    # 防呆（避免空資料爆炸）
    if ts.empty:
        return {}, f"No data for Id {selected_id}"

    # ===== 建立日期 + 排序 =====
    ts["date"] = pd.to_datetime(
        ts["year_number"].astype(str) + "-" + ts["month_number"].astype(str)
    )
    ts = ts.sort_values("date")

    # ===== 建立時間 index =====
    ts["time_index"] = ts["year_number"] * 12 + ts["month_number"]

    # ===== 計算 trend（線性回歸）=====
    if len(ts) > 1:
        coef = np.polyfit(ts["time_index"], ts["abundance_filtered_mean"], 1)
        ts["trend_line"] = coef[0] * ts["time_index"] + coef[1]
        slope = coef[0]
    else:
        ts["trend_line"] = np.nan
        slope = np.nan

    # ===== rolling mean（平滑）=====
    ts["rolling_mean"] = ts["abundance_filtered_mean"].rolling(3).mean()

    # ===== 畫圖 =====
    fig = px.line(
        ts,
        x="date",
        y="abundance_filtered_mean",
        title=f"Trend for Id {selected_id}",
        markers=True
    )

    # 加 trend 線
    fig.add_scatter(
        x=ts["date"],
        y=ts["trend_line"],
        mode="lines",
        name="Trend",
        line=dict(dash="dash", color="red")
    )

    # 加平滑線
    fig.add_scatter(
        x=ts["date"],
        y=ts["rolling_mean"],
        mode="lines",
        name="3-month avg",
        line=dict(color="green")
    )

    # ===== 統計 =====
    abd_range = ts["abundance_filtered_mean"].max() - ts["abundance_filtered_mean"].min()

    n_samples = checklist_df[checklist_df["Id"] == selected_id].shape[0]

    stats = html.Div([
        html.H4("Statistics"),
        html.P(f"Id: {selected_id}"),
        html.P(f"Mean: {ts['abundance_filtered_mean'].mean():.2f}"),
        html.P(f"Slope: {slope:.3f}"),
        html.P(f"Range: {abd_range:.2f}"),
        html.P(f"Samples: {n_samples}")
    ])

    return fig, stats

# 啟動
if __name__ == "__main__":
    app.run(debug=True)