from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import plotly.graph_objects as go


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Interactive bird abundance dashboard for project outputs."
    )
    parser.add_argument(
        "--eu-shp-path",
        default=str(
            script_dir
            / "EU_100km_fishnet_simple_by_distance"
            / "EU_100km_fishnet_simple_by_distance.shp"
        ),
        help="Path to the EU grid shapefile.",
    )
    parser.add_argument(
        "--checklist-folder",
        default=str(script_dir / "ebird_filtered_checklist"),
        help="Folder containing checklist CSV files.",
    )
    parser.add_argument(
        "--abundance-folder",
        default=str(script_dir / "validation_prediction_summary" / "mad_filter_abundance"),
        help="Folder containing MAD-filter abundance CSV files.",
    )
    parser.add_argument(
        "--species",
        default=None,
        help="Initial bird species scientific name. Defaults to first available species.",
    )
    parser.add_argument("--host", default="127.0.0.1", help="Dash host.")
    parser.add_argument("--port", type=int, default=8050, help="Dash port.")
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run Dash in debug mode.",
    )
    return parser.parse_args()


def extract_species_name_from_file(filename: str) -> str:
    name = filename.removesuffix(".csv")
    name = name.removeprefix("^")
    if name.endswith("_filtered_2019to2022"):
        name = name[: -len("_filtered_2019to2022")]
    if name.endswith("_filtered2019to2022"):
        name = name[: -len("_filtered2019to2022")]
    return name.strip()


def build_species_catalog(checklist_folder: Path, abundance_folder: Path) -> pd.DataFrame:
    checklist_paths = sorted(checklist_folder.glob("*.csv"))
    abundance_paths = sorted(abundance_folder.glob("*.csv"))

    checklist_df = pd.DataFrame(
        {
            "species": [extract_species_name_from_file(path.name) for path in checklist_paths],
            "checklist_path": [str(path) for path in checklist_paths],
        }
    )
    abundance_df = pd.DataFrame(
        {
            "species": [path.stem for path in abundance_paths],
            "abundance_path": [str(path) for path in abundance_paths],
        }
    )

    catalog = checklist_df.merge(abundance_df, on="species", how="inner")
    catalog = catalog.sort_values("species").reset_index(drop=True)
    if catalog.empty:
      raise ValueError(
          "No overlapping species found between checklist_folder and abundance_folder."
      )
    return catalog


def load_base_map(eu_shp_path: Path) -> gpd.GeoDataFrame:
    eu_shp = gpd.read_file(eu_shp_path)
    eu_shp["lon"] = eu_shp.geometry.centroid.x
    eu_shp["lat"] = eu_shp.geometry.centroid.y
    return eu_shp


def build_species_result(eu_shp: gpd.GeoDataFrame, checklist_path: Path, abundance_path: Path) -> tuple[gpd.GeoDataFrame, pd.DataFrame, pd.DataFrame]:
    checklist_df = pd.read_csv(checklist_path)
    abd_df = pd.read_csv(abundance_path)

    abd_df["time_index"] = abd_df["year_number"] * 12 + abd_df["month_number"]

    start = 2021 * 12 + 1
    end = 2022 * 12 + 12
    abd_sub = abd_df[(abd_df["time_index"] >= start) & (abd_df["time_index"] <= end)].copy()

    mean_df = abd_sub.groupby("Id")["abundance_filtered_mean"].mean().reset_index()
    mean_df.rename(columns={"abundance_filtered_mean": "mean_abundance"}, inplace=True)

    def calc_slope(group: pd.DataFrame) -> float:
        x = group["time_index"].to_numpy()
        y = group["abundance_filtered_mean"].to_numpy()
        if len(x) < 2:
            return np.nan
        return np.polyfit(x, y, 1)[0]

    trend_df = abd_sub.groupby("Id").apply(calc_slope).reset_index(name="trend")
    diff_df = abd_sub.groupby("Id")["abundance_filtered_mean"].agg(abd_max="max", abd_min="min").reset_index()
    diff_df["abd_range"] = diff_df["abd_max"] - diff_df["abd_min"]
    count_df = checklist_df.groupby("Id").size().reset_index(name="n_samples")

    result = eu_shp.merge(mean_df, on="Id", how="left")
    result = result.merge(trend_df, on="Id", how="left")
    result = result.merge(diff_df[["Id", "abd_range"]], on="Id", how="left")
    result = result.merge(count_df, on="Id", how="left")

    result["log_samples"] = np.log1p(result["n_samples"].fillna(0))
    result["log_abd_range"] = np.log1p(result["abd_range"].fillna(0))
    result["trend_clip"] = result["trend"].clip(-5, 5)
    result["mean_abundance"] = result["mean_abundance"].fillna(0)

    return result, checklist_df, abd_df


def build_map_figure(result: gpd.GeoDataFrame, species: str) -> go.Figure:
    map_fig = px.scatter_map(
        result,
        lat="lat",
        lon="lon",
        color="mean_abundance",
        size="log_samples",
        color_continuous_scale="viridis_r",
        size_max=20,
        zoom=3,
        height=650,
        hover_data=["Id", "mean_abundance", "log_samples", "trend_clip"],
        custom_data=["Id"],
        title=f"{species} mean abundance by grid",
    )
    map_fig.update_layout(mapbox_style="carto-positron", margin=dict(l=10, r=10, t=50, b=10))
    return map_fig


def build_empty_trend_figure() -> go.Figure:
    fig = go.Figure()
    fig.update_layout(
        title="Click a grid to inspect the time series",
        xaxis_title="Date",
        yaxis_title="Abundance",
        template="plotly_white",
    )
    return fig


def create_dashboard(args: argparse.Namespace) -> Dash:
    eu_shp_path = Path(args.eu_shp_path)
    checklist_folder = Path(args.checklist_folder)
    abundance_folder = Path(args.abundance_folder)

    catalog = build_species_catalog(checklist_folder, abundance_folder)
    eu_shp = load_base_map(eu_shp_path)

    initial_species = args.species or catalog.iloc[0]["species"]
    if initial_species not in set(catalog["species"]):
        raise ValueError(f"Species not found in dashboard catalog: {initial_species}")

    species_cache: dict[str, tuple[gpd.GeoDataFrame, pd.DataFrame, pd.DataFrame]] = {}

    def get_species_payload(species: str) -> tuple[gpd.GeoDataFrame, pd.DataFrame, pd.DataFrame]:
        if species not in species_cache:
            row = catalog.loc[catalog["species"] == species].iloc[0]
            species_cache[species] = build_species_result(
                eu_shp.copy(),
                Path(row["checklist_path"]),
                Path(row["abundance_path"]),
            )
        return species_cache[species]

    initial_result, _, _ = get_species_payload(initial_species)
    initial_map = build_map_figure(initial_result, initial_species)

    app = Dash(__name__)
    app.title = "Bird Abundance Dashboard"

    app.layout = html.Div(
        [
            html.H2("Spatiotemporal Bird Analysis Dashboard"),
            html.Div(
                [
                    html.Label("Bird species"),
                    dcc.Dropdown(
                        id="species_dropdown",
                        options=[{"label": sp, "value": sp} for sp in catalog["species"]],
                        value=initial_species,
                        clearable=False,
                    ),
                ],
                style={"width": "36%", "marginBottom": "16px"},
            ),
            html.Div(
                [
                    html.Div(
                        [dcc.Graph(id="map", figure=initial_map)],
                        style={"width": "52%", "display": "inline-block", "verticalAlign": "top"},
                    ),
                    html.Div(
                        [
                            dcc.Graph(id="trend_plot", figure=build_empty_trend_figure()),
                            html.Div(id="stats_output", children="Click a grid"),
                        ],
                        style={"width": "48%", "display": "inline-block", "verticalAlign": "top"},
                    ),
                ]
            ),
        ],
        style={"padding": "12px"},
    )

    @app.callback(Output("map", "figure"), Input("species_dropdown", "value"))
    def update_map(species: str) -> go.Figure:
        result, _, _ = get_species_payload(species)
        return build_map_figure(result, species)

    @app.callback(
        Output("trend_plot", "figure"),
        Output("stats_output", "children"),
        Input("map", "clickData"),
        Input("species_dropdown", "value"),
    )
    def update_output(click_data: dict | None, species: str):
        result, checklist_df, abd_df = get_species_payload(species)

        if click_data is None:
            return build_empty_trend_figure(), "Click a grid"

        point = click_data["points"][0]
        selected_id = point["customdata"][0] if "customdata" in point else point.get("hovertext")
        if pd.isna(selected_id):
            lon = point["lon"]
            lat = point["lat"]
            idx = ((result["lon"] - lon) ** 2 + (result["lat"] - lat) ** 2).idxmin()
            selected_id = result.loc[idx, "Id"]

        ts = abd_df[abd_df["Id"] == selected_id].copy()
        if ts.empty:
            return build_empty_trend_figure(), f"No data for Id {selected_id}"

        ts["date"] = pd.to_datetime(
            ts["year_number"].astype(str) + "-" + ts["month_number"].astype(str) + "-01"
        )
        ts = ts.sort_values("date")
        ts["time_index"] = ts["year_number"] * 12 + ts["month_number"]

        if len(ts) > 1:
            coef = np.polyfit(ts["time_index"], ts["abundance_filtered_mean"], 1)
            ts["trend_line"] = coef[0] * ts["time_index"] + coef[1]
            slope = coef[0]
        else:
            ts["trend_line"] = np.nan
            slope = np.nan

        ts["rolling_mean"] = ts["abundance_filtered_mean"].rolling(3).mean()

        fig = px.line(
            ts,
            x="date",
            y="abundance_filtered_mean",
            title=f"{species} trend for Id {selected_id}",
            markers=True,
        )
        fig.add_scatter(
            x=ts["date"],
            y=ts["trend_line"],
            mode="lines",
            name="Trend",
            line=dict(dash="dash", color="red"),
        )
        fig.add_scatter(
            x=ts["date"],
            y=ts["rolling_mean"],
            mode="lines",
            name="3-month avg",
            line=dict(color="green"),
        )
        fig.update_layout(template="plotly_white", margin=dict(l=20, r=20, t=50, b=20))

        abd_range = ts["abundance_filtered_mean"].max() - ts["abundance_filtered_mean"].min()
        n_samples = checklist_df[checklist_df["Id"] == selected_id].shape[0]
        selected_row = result[result["Id"] == selected_id].iloc[0]

        stats = html.Div(
            [
                html.H4("Statistics"),
                html.P(f"Species: {species}"),
                html.P(f"Id: {selected_id}"),
                html.P(f"Mean abundance: {ts['abundance_filtered_mean'].mean():.2f}"),
                html.P(f"Slope: {slope:.3f}" if not np.isnan(slope) else "Slope: NA"),
                html.P(f"Range: {abd_range:.2f}"),
                html.P(f"Checklist samples: {n_samples}"),
                html.P(f"Map mean abundance: {selected_row['mean_abundance']:.2f}"),
            ]
        )

        return fig, stats

    return app


if __name__ == "__main__":
    args = parse_args()
    app = create_dashboard(args)
    app.run(host=args.host, port=args.port, debug=args.debug)
