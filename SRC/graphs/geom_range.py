from COMMON.gnss_constants import S_IN_H
from COMMON.constants import LOS_IDX


def geom_range(nav_solution, los_data):
    if nav_solution == 'GPS':
        constellation_id = 'G'
        constellation = 'GPS'
    else:
        constellation_id = 'E'
        constellation = 'Galileo'

    PlotConf = {
        "Title": "Simulated GNSS " + constellation + " Geometrical Range",
        "Type": "Lines",
        "Grid": True,
        "FigSize": (8.4, 7.6),
        "Marker": ".",
        "LineWidth": 0.2,
    }

    PlotConf["xLabel"] = "Hour of Day "
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["yLabel"] = "Geometrical Range [km]"

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elev"
    PlotConf["ColorBarMin"] = 0.0
    PlotConf["ColorBarMax"] = 90.0

    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    filter_condition = los_data[LOS_IDX["CONST"]] == constellation_id

    Label = 0
    PlotConf["xData"] = los_data[LOS_IDX["SOD"]][filter_condition] / S_IN_H
    PlotConf["yData"][Label] = los_data[LOS_IDX["GEOM-RANGE"]
                                        ][filter_condition]
    PlotConf["zData"][Label] = los_data[LOS_IDX["ELEV"]][filter_condition]

    return PlotConf
