from datetime import datetime

import param
import panel as pn
import numpy as np
import pandas as pd
import hvplot.pandas
import geoviews as gv
import holoviews as hv
from holoviews.streams import Tap
from bokeh.themes import Theme

pn.extension(throttled=True, notifications=True)

VAR_OPTIONS = {
    "Maximum Air Temperature [F]": "max_temp_f",
    "Minimum Air Temperature [F]": "min_temp_f",
    "Maximum Dew Point [F]": "max_dewpoint_f",
    "Minimum Dew Point [F]": "min_dewpoint_f",
    "Daily Precipitation [inch]": "precip_in",
    "Average Wind Speed [knots]": "avg_wind_speed_kts",
    "Average Wind Direction [deg]": "avg_wind_drct",
    "Minimum Relative Humidity [%]": "min_rh",
    "Average Relative Humidity [%]": "avg_rh",
    "Maximum Relative Humidity [%]": "max_rh",
    "NCEI 1991-2020 Daily High Temperature Climatology [F]": "climo_high_f",
    "NCEI 1991-2020 Daily Low Temperature Climatology [F]": "climo_low_f",
    "NCEI 1991-2020 Daily Precipitation Climatology [inch]": "climo_precip_in",
    "Reported Snowfall [inch]": "snow_in",
    "Reported Snow Depth [inch]": "snowd_in",
    "Minimum 'Feels Like' Temperature [F]": "min_feel",
    "Average 'Feels Like' Temperature [F]": "avg_feel",
    "Maximum 'Feels Like' Temperature [F]": "max_feel",
    "Maximum sustained wind speed [knots]": "max_wind_speed_kts",
    "Maximum wind gust [knots]": "max_wind_gust_kts",
    "Daily Solar Radiation MJ/m2": "srad_mj",
}
WELCOME_MESSAGE = """
### Welcome to the Climaviewer!

This app allows you to compare a single year of data from a weather station to the average of a range of years.

1. Select a network station from the dropdowns; alternatively you can click on the map to select the nearest station.
2. Choose the desired variable and year to plot, and change the average range if desired.
3. Hover over the plot to see the value for each day; use the mouse wheel to zoom in/out.
4. Tap on the top plot to see the value for a specific day, compared to all other years in the average range.
"""

FOOTER_MESSAGE = """
Made entirely with OSS packages: [Panel](https://panel.holoviz.org/), [Holoviews](https://holoviews.org/), [GeoViews](https://geoviews.org/), [Bokeh](https://bokeh.org/), [Pandas](https://pandas.pydata.org/), [Numpy](https://numpy.org/)

Data sourced from the [Iowa Environmental Mesonet](https://mesonet.agron.iastate.edu/).
"""

VAR_OPTIONS_R = {v: k for k, v in VAR_OPTIONS.items()}
ONI_URL = "https://raw.githubusercontent.com/ahuang11/oni/master/oni.csv"
NETWORKS_URL = "https://mesonet.agron.iastate.edu/sites/networks.php?network=_ALL_&format=csv&nohtml=on"
STATION_URL_FMT = (
    "https://mesonet.agron.iastate.edu/cgi-bin/request/daily.py?network={network}&stations={station}"
    "&year1=1928&month1=1&day1=1&year2=2024&month2=12&day2=31&var={var}&na=blank&format=csv"
)
DARK_RED = "#FF5555"
DARK_BLUE = "#5588FF"
SEASON_TO_MONTH = {
    "DJF": "JAN",
    "JFM": "FEB",
    "FMA": "MAR",
    "MAM": "APR",
    "AMJ": "MAY",
    "MJJ": "JUN",
    "JJA": "JUL",
    "JAS": "AUG",
    "ASO": "SEP",
    "SON": "OCT",
    "OND": "NOV",
    "NDJ": "DEC",
}
XTICKS = [
    (1, "JAN"),
    (31, "FEB"),
    (59, "MAR"),
    (90, "APR"),
    (120, "MAY"),
    (151, "JUN"),
    (181, "JUL"),
    (212, "AUG"),
    (243, "SEP"),
    (273, "OCT"),
    (304, "NOV"),
    (334, "DEC"),
]
MONTH_TO_JULIAN_DAY = {month: day for day, month in XTICKS}
ONI_COLORS = {
    "El Nino": DARK_RED,
    "Neutral": "grey",
    "La Nina": DARK_BLUE,
}

THEME_JSON = {
    "attrs": {
        "figure": {
            "background_fill_color": "#1b1e23",
            "border_fill_color": "#1b1e23",
            "outline_line_alpha": 0,
        },
        "Grid": {
            "grid_line_color": "#808080",
            "grid_line_alpha": 0.1,
        },
        "Axis": {
            # tick color and alpha
            "major_tick_line_color": "#4d4f51",
            "minor_tick_line_alpha": 0,
            # tick labels
            "major_label_text_font": "Courier New",
            "major_label_text_color": "#808080",
            "major_label_text_align": "left",
            "major_label_text_font_size": "0.95em",
            "major_label_text_font_style": "normal",
            # axis labels
            "axis_label_text_font": "Courier New",
            "axis_label_text_font_style": "normal",
            "axis_label_text_font_size": "1.15em",
            "axis_label_text_color": "lightgrey",
            "axis_line_color": "#4d4f51",
        },
        "Legend": {
            "spacing": 8,
            "glyph_width": 15,
            "label_standoff": 8,
            "label_text_color": "#808080",
            "label_text_font": "Courier New",
            "label_text_font_size": "0.95em",
            "label_text_font_style": "bold",
            "border_line_alpha": 0,
            "background_fill_alpha": 0.25,
            "background_fill_color": "#1b1e23",
        },
        "BaseColorBar": {
            # axis labels
            "title_text_color": "lightgrey",
            "title_text_font": "Courier New",
            "title_text_font_size": "0.95em",
            "title_text_font_style": "normal",
            # tick labels
            "major_label_text_color": "#808080",
            "major_label_text_font": "Courier New",
            "major_label_text_font_size": "0.95em",
            "major_label_text_font_style": "normal",
            "background_fill_color": "#1b1e23",
            "major_tick_line_alpha": 0,
            "bar_line_alpha": 0,
        },
        "Title": {
            "text_font": "Courier New",
            "text_font_style": "normal",
            "text_color": "lightgrey",
        },
    }
}
theme = Theme(json=THEME_JSON)
this_year = datetime.now().year

hv.renderer("bokeh").theme = theme


class ClimateApp(pn.viewable.Viewer):
    network = param.Selector(
        default="WA_ASOS", label="Network (delete & type to search)"
    )
    station = param.Selector(default="SEA", label="Station (delete & type to search)")
    year = param.Integer(default=this_year - 1, bounds=(1928, this_year))
    year_range = param.Range(
        default=(1990, 2020), bounds=(1928, this_year), label="Average Range"
    )
    var = param.Selector(default="max_temp_f", objects=sorted(VAR_OPTIONS.values()))
    stat = param.Selector(default="Mean", objects=["Mean", "Median"])

    _title = param.String()
    _ylabel = param.String()

    def __init__(self, **params):
        super().__init__(**params)
        self._sidebar = pn.Column(sizing_mode="stretch_both")
        self._main = pn.Column(
            pn.indicators.LoadingSpinner(
                value=True, width=25, height=25, name="Loading, please wait a moment..."
            ),
            sizing_mode="stretch_both",
        )
        self._modal = pn.Column(width=850, height=500, align="center")
        self._template = pn.template.FastListTemplate(
            sidebar=[self._sidebar],
            main=[self._main],
            modal=[self._modal],
            theme="dark",
            theme_toggle=False,
            main_layout=None,
            title="Year vs Climatology",
            accent="grey",
        )
        pn.state.onload(self._onload)

    def _onload(self):
        try:
            self._sidebar.loading = True
            self._populate_sidebar()
            self._populate_main()
            self._populate_modal()
        finally:
            self._sidebar.loading = False

    def _populate_sidebar(self):
        self._network_df = self._get_network_df()
        self._oni_df = self._get_oni_df()
        networks = sorted(self._network_df["iem_network"].unique())
        self.param["network"].objects = networks

        open_button = pn.widgets.Button(
            name="Select station from map",
            button_type="primary",
            sizing_mode="stretch_width",
        )
        open_button.on_click(self._open_modal)
        network_select = pn.widgets.AutocompleteInput.from_param(
            self.param.network, min_characters=0, case_sensitive=False
        )
        station_select = pn.widgets.AutocompleteInput.from_param(
            self.param.station, min_characters=0, case_sensitive=False
        )
        var_select = pn.widgets.Select.from_param(self.param.var, options=VAR_OPTIONS)
        year_slider = pn.widgets.IntSlider.from_param(self.param.year)
        year_range_slider = pn.widgets.RangeSlider.from_param(self.param.year_range)
        stat_select = pn.widgets.RadioButtonGroup.from_param(
            self.param.stat, sizing_mode="stretch_width"
        )
        self._sidebar.objects = [
            pn.pane.Markdown(WELCOME_MESSAGE),
            open_button,
            network_select,
            station_select,
            var_select,
            year_slider,
            year_range_slider,
            stat_select,
            pn.pane.Markdown(FOOTER_MESSAGE),
        ]

    def _populate_main(self):
        self._tap_x = Tap()
        self._station_pane = pn.pane.HoloViews(
            sizing_mode="stretch_both",
            min_width=800,
            min_height=350,
        )
        self._day_pane = pn.pane.HoloViews(
            sizing_mode="stretch_both",
            min_width=800,
            min_height=350,
        )
        self._update_stations()
        self._update_var_station_dependents()
        self._day_pane.object = pn.bind(
            self._update_day_plot,
            self.param.var,
            self.param.station,
            self.param.year,
            self.param.year_range,
            self._tap_x.param.x,
        )
        pointer_vline = hv.DynamicMap(self._update_vline, streams=[self._tap_x])
        oni_plot = hv.DynamicMap(
            self._update_oni_plot
        )
        station_plot = hv.DynamicMap(
            pn.bind(
                self._update_station_plot,
                self.param.var,
                self.param.station,
                self.param.year,
                self.param.year_range,
                self.param.stat,
            ),
        )
        self._station_pane.object = (
            (oni_plot * station_plot * pointer_vline)
            .opts(
                xlabel="Time of Year",
                gridstyle={"ygrid_line_alpha": 0},
                xticks=XTICKS,
                show_grid=True,
                padding=(0, (0, 0.45)),
                responsive=True,
                shared_axes=False,
                legend_position="top_right"
            )
            .apply.opts(title=self.param._title, ylabel=self.param._ylabel)
        )
        self._update_ylabel()
        self._update_title()

        self._main.objects = [self._station_pane, self._day_pane]

    def _populate_modal(self):
        network_points = self._network_df.hvplot.points(
            "lon",
            "lat",
            legend=False,
            cmap="category10",
            color="iem_network",
            hover_cols=["stid", "station_name", "iem_network"],
            size=10,
            geo=True,
            responsive=True,
            xlabel="Longitude",
            ylabel="Latitude",
        ).opts(
            "Points",
            fill_alpha=0,
            tools=["tap", "hover"],
            active_tools=["wheel_zoom"],
        )

        tap = Tap(source=network_points)
        pn.bind(self._update_station, x=tap.param.x, y=tap.param.y, watch=True)
        instructions = pn.pane.Markdown(
            "#### The nearest station will be selected when you click on the map."
        )
        network_pane = pn.pane.HoloViews(
            network_points * gv.tile_sources.CartoDark(),
        )
        self._modal.objects = [instructions, network_pane]

    def _open_modal(self, event):
        self._template.open_modal()

    @pn.cache
    def _get_oni_df(self):
        df = pd.read_csv(ONI_URL)
        df["month"] = df["season"].map(SEASON_TO_MONTH)
        df["julian_day"] = df["month"].map(MONTH_TO_JULIAN_DAY)
        df["julian_day_end"] = df["julian_day"].shift(-1).fillna(365)
        df["oni"] = df["oni"].str.replace("_", " ").str.title()
        return df

    @pn.cache
    def _get_network_df(self):
        network_df = pd.read_csv(NETWORKS_URL)
        return network_df.loc[network_df["iem_network"].str.contains("ASOS")]

    @pn.depends("network", watch=True)
    def _update_stations(self):
        self._template.close_modal()
        network_df_subset = self._network_df.loc[
            self._network_df["iem_network"] == self.network,
            ["stid", "station_name"],
        ]
        names = sorted(network_df_subset["station_name"].unique())
        stids = sorted(network_df_subset["stid"].unique())
        self.param["station"].objects = names + stids

    def _update_station(self, x, y):
        if x is None or y is None:
            return

        def haversine_vectorized(lon1, lat1, lon2, lat2):
            R = 6371  # Radius of the Earth in kilometers
            dlat = np.radians(lat2 - lat1)
            dlon = np.radians(lon2 - lon1)
            a = (
                np.sin(dlat / 2.0) ** 2
                + np.cos(np.radians(lat1))
                * np.cos(np.radians(lat2))
                * np.sin(dlon / 2.0) ** 2
            )
            c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
            return R * c

        distances = haversine_vectorized(
            self._network_df["lon"].values, self._network_df["lat"].values, x, y
        )

        min_distance_index = np.argmin(distances)

        closest_row = self._network_df.iloc[min_distance_index]
        with param.parameterized.batch_call_watchers(self):
            self.network = closest_row["iem_network"]
            self.station = closest_row["stid"]

    @pn.cache
    def _get_station_df(self, station, var):
        if station in self._network_df["station_name"].unique():
            station = self._network_df.loc[
                self._network_df["station_name"] == station, "stid"
            ].iloc[0]
            if station.startswith("K"):
                station = station.lstrip("K")
        station_url = STATION_URL_FMT.format(
            network=self.network, station=station, var=var
        )
        station_df = (
            pd.read_csv(
                station_url,
                parse_dates=True,
                index_col="day",
            )
            .drop(columns=["station"])
            .astype("float16")
            .assign(
                dayofyear=lambda df: df.index.dayofyear,
                year=lambda df: df.index.year,
            )
            .dropna()
        )
        return station_df

    @pn.depends("var", "station", watch=True)
    def _update_var_station_dependents(self):
        try:
            self._main.loading = True
            self._station_df = self._get_station_df(self.station, self.var).dropna()
            if len(self._station_df) == 0:
                return

            year_range_min = self._station_df["year"].min()
            year_range_max = self._station_df["year"].max()
            self.param["year"].bounds = (year_range_min, year_range_max)
            if self.year < year_range_min:
                self.year = year_range_min
            if self.year > year_range_max:
                self.year = year_range_max
        finally:
            self._main.loading = False

    def _update_vline(self, x, y):
        if x is None:
            x = 0
        if y is None:
            y = 0
        vline = hv.VLine(x).opts(line_width=0.9, color="lightgrey")
        text = hv.Text(
            x,
            y,
            f"Julian Day {int(x)}",
        ).opts(
            text_color="lightgrey",
            text_align="left",
            text_baseline="bottom",
            text_alpha=0.8,
        )
        return vline * text

    @param.depends("year", watch=True)
    def _update_oni_plot(self):
        df = self._oni_df
        df_year = df.loc[df["year"] == 1998]
        df_year.iloc[-1, -1] = 365

        overlay = hv.Overlay([])
        for oni in df_year["oni"].unique():
            df_subset = df_year.loc[df_year["oni"] == oni, ["julian_day", "julian_day_end"]]
            overlay *= hv.VSpans(
                df_subset,
                ["julian_day", "julian_day_end"],
                label=oni,
            ).opts(color=ONI_COLORS[oni], alpha=0.18, line_alpha=0)
        return overlay

    def _update_day_plot(self, var, station, year, year_range, x):
        df = self._station_df
        if not x:
            x = 1
        x = int(x)
        df_day_year = df.query(f"year == {year}")
        if x > df_day_year["dayofyear"].max():
            x = df_day_year["dayofyear"].max()
        day_year = df_day_year.loc[df_day_year["dayofyear"] == x, self.var].iloc[0]
        df_subset = df.loc[df["year"].between(*year_range)]
        df_day_climo = df_subset.loc[df_subset["dayofyear"] == x]
        df_day_climo = df_day_climo.assign(
            above_or_below=df_day_climo[self.var] >= day_year
        )
        title = (
            f"{VAR_OPTIONS_R[self.var]} across {year_range[0]}-{year_range[1]} on "
            + df_day_climo.index.strftime("%B %d")[0]
            + f" (Julian Day {x}) "
        )

        days_above = df_day_climo.loc[df_day_climo["above_or_below"] == True].shape[0]
        days_below = df_day_climo.loc[df_day_climo["above_or_below"] == False].shape[0]

        min_x = df[self.var].min()

        plot = hv.Overlay([])
        plot *= df_day_climo.hvplot.hist(
            self.var,
            responsive=True,
            by="above_or_below",
            bins=11,
            legend=False,
            color=hv.Cycle([DARK_BLUE, DARK_RED]),
        ).opts("Histogram", fill_alpha=0.7, line_alpha=0)
        plot *= hv.VLine(day_year).opts(line_width=0.9, color="lightgrey")
        plot *= hv.Text(
            day_year,
            0.1,
            f"{year}",
        ).opts(
            text_color="lightgrey",
            text_align="left",
            text_baseline="bottom",
            text_alpha=0.8,
        )
        plot *= self._create_text_days_labels(
            df,
            days_above,
            days_below,
            text_x=min_x + 5,
            text_y=4,
            spacing=1,
            suffix=f"YEARS",
        )

        return plot.opts(
            xlabel=VAR_OPTIONS_R[self.var],
            ylabel="Number of Days",
            title=title,
            shared_axes=False,
            show_grid=True,
            gridstyle={"xgrid_line_alpha": 0},
            xlim=(min_x, df[self.var].max()),
        )

    def _update_station_plot(self, var, station, year, year_range, stat):
        if len(self._station_df) == 0:
            return

        # base dataframes
        df = self._station_df
        df_subset = df.loc[df["year"].between(*year_range)]
        df_avg = df_subset.groupby("dayofyear").mean()
        df_year = df[df.year == year]

        # above/below
        df_year = df_year[["dayofyear", self.var]].merge(
            df_avg.reset_index()[["dayofyear", self.var]],
            on="dayofyear",
            suffixes=("", "_avg"),
        )
        df_year["above_or_below"] = df_year[self.var] >= df_year[f"{self.var}_avg"]
        days_above = df_year.loc[df_year["above_or_below"] == True].shape[0]
        days_below = df_year.loc[df_year["above_or_below"] == False].shape[0]

        # stats
        if stat == "Mean":
            year_avg = df_year[self.var].mean()
        else:
            year_avg = df_year[self.var].median()
        year_max = df_year[self.var].max()
        year_min = df_year[self.var].min()

        plots = self._create_line_plots(df, df_year, df_avg)
        lines = self._create_hlines(year_avg, year_max, year_min)
        texts = self._create_text_labels(year_avg, year_max, year_min)
        text_days = self._create_text_days_labels(df, days_above, days_below)

        # Overlay all elements
        station_overlay = plots * lines * texts * text_days
        return station_overlay

    @pn.depends("var", watch=True)
    def _update_ylabel(self):
        self._ylabel = VAR_OPTIONS_R[self.var]

    @pn.depends("station", "year", "year_range", watch=True)
    def _update_title(self):
        df = self._station_df
        df_subset = df.loc[df["year"].between(*self.year_range)]
        # hack to get the title and ylabel to update
        year_min = df_subset["year"].min()
        if self.year_range[0] > year_min:
            year_min = self.year_range[0]
        year_max = df_subset["year"].max()
        if self.year_range[1] < year_max:
            year_max = self.year_range[1]
        year_range_label = f"{year_min}-{year_max}"
        self._title = f"{self._get_station_label()} - {self.year} vs Average ({year_range_label})"

    def _create_line_plots(self, df, df_year, df_avg):
        plot_kwargs = {
            "x": "dayofyear",
            "y": self.var,
            "legend": False,
            "responsive": True,
        }
        plot = df.hvplot(
            by="year",
            color="grey",
            alpha=0.02,
            hover=False,
            **plot_kwargs,
        )
        df_above = df_year.copy()
        df_above.loc[df_above["above_or_below"]] = np.nan
        df_below = df_year.copy()
        df_below.loc[~df_below["above_or_below"]] = np.nan
        plot_year = df_year.hvplot(
            color="lightgrey", hover="vline", alpha=0.5, **plot_kwargs
        ).redim.label(**{"dayofyear": "Julian Day", self.var: str(self.year)})
        plot_above = df_above.hvplot(
            hover="vline", color=DARK_BLUE, **plot_kwargs
        ).redim.label(**{"dayofyear": "Julian Day", self.var: str(self.year)})
        plot_below = df_below.hvplot(
            hover="vline", color=DARK_RED, **plot_kwargs
        ).redim.label(**{"dayofyear": "Julian Day", self.var: str(self.year)})
        plot_avg = df_avg.hvplot(
            color="lightgrey", hover="vline", **plot_kwargs
        ).redim.label(**{"dayofyear": "Julian Day", self.var: "Average"})
        return plot * plot_year * plot_above * plot_below * plot_avg

    def _create_hlines(self, year_avg, year_max, year_min):
        # Create horizontal lines
        plot_year_avg = hv.HLine(year_avg).opts(
            line_color="lightgrey", line_dash="dashed", line_width=0.5
        )
        plot_year_max = hv.HLine(year_max).opts(
            line_color=DARK_RED, line_dash="dashed", line_width=0.5
        )
        plot_year_min = hv.HLine(year_min).opts(
            line_color=DARK_BLUE, line_dash="dashed", line_width=0.5
        )
        return plot_year_avg * plot_year_max * plot_year_min

    def _create_text_labels(self, year_avg, year_max, year_min):
        text_year_opts = {
            "text_align": "right",
            "text_baseline": "bottom",
            "text_alpha": 0.8,
        }
        text_year_label = "AVERAGE" if self.stat == "Mean" else "MEDIAN"
        text_year_avg = hv.Text(
            360, year_avg + 3, f"{text_year_label} {year_avg:.1f}", fontsize=8
        ).opts(
            text_color="lightgrey",
            **text_year_opts,
        )
        text_year_max = hv.Text(
            360, year_max + 3, f"MAX {year_max:.1f}", fontsize=8
        ).opts(
            text_color=DARK_RED,
            **text_year_opts,
        )
        text_year_min = hv.Text(
            360, year_min + 3, f"MIN {year_min:.1f}", fontsize=8
        ).opts(
            text_color=DARK_BLUE,
            **text_year_opts,
        )
        return text_year_avg * text_year_max * text_year_min

    def _create_areas(self, df_above, df_below):
        area_kwargs = {
            "x": "dayofyear",
            "y": f"{self.var}_avg",
            "y2": self.var,
            "hover": False,
            "responsive": True,
        }
        area_opts = {"fill_alpha": 0.2, "line_alpha": 0.8}
        plot_above = df_above.hvplot.area(**area_kwargs).opts(
            line_color=DARK_RED, fill_color=DARK_RED, **area_opts
        )
        plot_below = df_below.hvplot.area(**area_kwargs).opts(
            line_color=DARK_BLUE, fill_color=DARK_BLUE, **area_opts
        )
        return plot_above * plot_below

    def _create_text_days_labels(
        self,
        df,
        days_above,
        days_below,
        text_x=None,
        text_y=None,
        spacing=None,
        suffix=None,
    ):
        text_x = text_x or 30
        text_y = text_y or df[self.var].max() + 3
        spacing = spacing or 2
        suffix = suffix or "DAYS"
        text_days_above = hv.Text(text_x, text_y, f"{days_above}", fontsize=14).opts(
            text_align="right",
            text_baseline="bottom",
            text_color=DARK_RED,
            text_alpha=0.8,
        )
        text_days_below = hv.Text(text_x, text_y, f"{days_below}", fontsize=14).opts(
            text_align="right",
            text_baseline="top",
            text_color=DARK_BLUE,
            text_alpha=0.8,
        )
        text_above = hv.Text(
            text_x + spacing, text_y, f"{suffix} ABOVE", fontsize=7
        ).opts(
            text_align="left",
            text_baseline="bottom",
            text_color="lightgrey",
            text_alpha=0.8,
        )
        text_below = hv.Text(
            text_x + spacing, text_y, f"{suffix} BELOW", fontsize=7
        ).opts(
            text_align="left",
            text_baseline="top",
            text_color="lightgrey",
            text_alpha=0.8,
        )

        return text_days_above * text_days_below * text_above * text_below

    def _get_station_label(self):
        if self.station not in self._network_df["station_name"].unique():
            stid = self.station
            station_name = self._network_df.loc[
                self._network_df["stid"] == self.station, "station_name"
            ].iloc[0]
        else:
            stid = self._network_df.loc[
                self._network_df["station_name"] == self.station, "stid"
            ].iloc[0]
            station_name = self.station
        station_label = f"{station_name.title()} ({stid})"
        return station_label

    def __panel__(self):
        return self._template


ClimateApp().servable()
