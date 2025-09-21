import os
import warnings
from mpl_toolkits.basemap import Basemap
import matplotlib.cbook
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import conda

CondaFileDir = conda.__file__
CondaDir = CondaFileDir.split("lib")[0]
ProjLib = os.path.join(os.path.join(CondaDir, "share"), "proj")
os.environ["PROJ_LIB"] = ProjLib


warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


def createFigure(PlotConf):
    try:
        fig, ax = plt.subplots(1, 1, figsize=PlotConf["FigSize"])
    except:
        fig, ax = plt.subplots(1, 1)

    return fig, ax


def saveFigure(fig, Path):
    Dir = os.path.dirname(Path)
    try:
        os.makedirs(Dir)
    except:
        pass
    fig.savefig(Path, dpi=150.0, bbox_inches="tight")


def prepareAxis(PlotConf, ax):
    ax.set_axisbelow(True)

    # for key in PlotConf:
    if "Title" in PlotConf:
        ax.set_title(PlotConf["Title"])

    # set x coordinate label
    if "xLabel" in PlotConf:
        ax.set_xlabel(PlotConf["xLabel"])

    if "xTicks" in PlotConf:
        ax.set_xticks(PlotConf["xTicks"])

    if "xTicksLabels" in PlotConf:
        ax.set_xticklabels(PlotConf["xTicksLabels"])

    if "xLim" in PlotConf:
        ax.set_xlim(PlotConf["xLim"])

    # set y coordinate label
    if "yLabel" in PlotConf:
        ax.set_ylabel(PlotConf["yLabel"])

    if "yFormat" in PlotConf:
        ax.yaxis.set_major_formatter(
            mpl.ticker.StrMethodFormatter(PlotConf["yFormat"]))

    if "yTicks" in PlotConf:
        ax.set_yticks(PlotConf["yTicks"])

    if "yTicksLabels" in PlotConf:
        ax.set_yticklabels(PlotConf["yTicksLabels"])

    if "yLim" in PlotConf:
        ax.set_ylim(PlotConf["yLim"])

    # set grid
    if "Grid" in PlotConf and PlotConf["Grid"] == True:
        ax.grid(linestyle="--", linewidth=0.5, which="both")


def prepareColorBar(PlotConf, ax, Values):
    try:
        Min = PlotConf["ColorBarMin"]
    except:
        Mins = []
        for v in Values.values():
            Mins.append(min(v))
        Min = min(Mins)
    try:
        Max = PlotConf["ColorBarMax"]
    except:
        Maxs = []
        for v in Values.values():
            Maxs.append(max(v))
        Max = max(Maxs)
    normalize = mpl.cm.colors.Normalize(vmin=Min, vmax=Max)

    divider = make_axes_locatable(ax)
    # size size% of the plot and gap of pad% from the plot
    color_ax = divider.append_axes("right", size="3%", pad="2%")

    cmap = mpl.cm.get_cmap(PlotConf["ColorBar"])
    mpl.colorbar.ColorbarBase(
        color_ax,
        cmap=cmap,
        norm=mpl.colors.Normalize(vmin=Min, vmax=Max),
        label=PlotConf["ColorBarLabel"],
        ticks=PlotConf["ColorBarTicks"] if "ColorBarTicks" in PlotConf else None,
    )

    return normalize, cmap


def drawMap(
    PlotConf,
    ax,
):
    Map = Basemap(
        projection="cyl",
        llcrnrlat=PlotConf["LatMin"] - 0,
        urcrnrlat=PlotConf["LatMax"] + 0,
        llcrnrlon=PlotConf["LonMin"] - 0,
        urcrnrlon=PlotConf["LonMax"] + 0,
        lat_ts=10,
        resolution="l",
        ax=ax,
    )

    # Draw map meridians
    Map.drawmeridians(
        np.arange(PlotConf["LonMin"], PlotConf["LonMax"] +
                  1, PlotConf["LonStep"]),
        labels=[0, 0, 0, 1],
        fontsize=6,
        linewidth=0.2,
        labelstyle="+/-",
    )

    # Draw map parallels
    Map.drawparallels(
        np.arange(PlotConf["LatMin"], PlotConf["LatMax"] +
                  1, PlotConf["LatStep"]),
        labels=[1, 0, 0, 0],
        fontsize=6,
        linewidth=0.2,
        labelstyle="+/-",
        fmt="%i",
    )

    # Draw coastlines
    Map.drawcoastlines(linewidth=0.5)

    # Draw countries
    Map.drawcountries(linewidth=0.25)


def prepareThirdAxis(PlotConf, ax):
    z_ax = ax.twinx()
    z_ax.set_ylabel(PlotConf["zAxisLabel"])

    if "zTicks" in PlotConf:
        z_ax.set_yticks(PlotConf["zTicks"])

    if "zLim" in PlotConf:
        z_ax.set_ylim(PlotConf["zLim"])

    ax_label = PlotConf["zLabel"]

    z_ax.plot(
        PlotConf["xData"],
        PlotConf["zData"][ax_label],
        color=PlotConf["Colors"][ax_label] if "Colors" in PlotConf else None,
        marker=PlotConf["Marker"],
        label=ax_label,
        linewidth=PlotConf["LineWidth"] if "LineWidth" in PlotConf else 1.5,
    )

    return z_ax


def generateLinesPlot(PlotConf, output):
    LineWidth = 1.5

    fig, ax = createFigure(PlotConf)

    prepareAxis(PlotConf, ax)
    z_ax = None

    if "LineWidth" in PlotConf:
        LineWidth = PlotConf["LineWidth"]

    if "ColorBar" in PlotConf:
        normalize, cmap = prepareColorBar(PlotConf, ax, PlotConf["zData"])
    elif "zLabel" in PlotConf:
        # If there is a zLabel, we prepare a third axis that is not a color bar
        z_ax = prepareThirdAxis(PlotConf, ax)

    if "Map" in PlotConf and PlotConf["Map"] == True:
        drawMap(PlotConf, ax)

    for Label in PlotConf["yData"].keys():
        if "ColorBar" in PlotConf:
            ax.scatter(
                PlotConf["xData"],
                PlotConf["yData"][Label],
                marker=PlotConf["Marker"],
                alpha=PlotConf["Alpha"] if "Alpha" in PlotConf else 1.0,
                linewidth=LineWidth,
                c=cmap(normalize(np.array(PlotConf["zData"][Label]))),
                label=Label,
            )

        else:
            ax.plot(
                PlotConf["xData"],
                PlotConf["yData"][Label],
                color=PlotConf["Colors"][Label] if "Colors" in PlotConf else None,
                marker=PlotConf["Marker"],
                linewidth=LineWidth,
                label=Label,
            )

    # Only add legend if position is specified
    if "LegendPosition" in PlotConf:
        lines1, labels1 = ax.get_legend_handles_labels()

        lines2, labels2 = [], []
        if z_ax is not None:
            lines2, labels2 = z_ax.get_legend_handles_labels()

        ax.legend(
            lines1 + lines2,
            labels1 + labels2,
            loc=(
                PlotConf["LegendPosition"] if "LegendPosition" in PlotConf else "best"
            ),
            frameon=True,
            fontsize=9,
            markerscale=2.0,
        )

    saveFigure(fig, output)

    plt.close(fig)


def generatePlot(PlotConf, output):
    if PlotConf["Type"] == "Lines":
        generateLinesPlot(PlotConf, output)
