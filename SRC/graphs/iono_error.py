import numpy as np
from COMMON.constants import LOS_IDX


def iono_error(nav_solution, los_data):
    if nav_solution == 'GPS':
        constellation_id = 'G'
        constellation = 'GPS'
    else:
        constellation_id = 'E'
        constellation = 'Galileo'
    PlotConf = {
        "Title": "Simulated GNSS " + constellation + " IONO Error",
        "Type": "Lines",
        "FigSize": (8.4, 7.6),
        "Grid": True,
        "Marker": ".",
        "lineWidth": 0.001,
        "Alpha": 0.3,
    }

    filter_condition = los_data[LOS_IDX["CONST"]] == constellation_id

    PlotConf["xLabel"] = "Elevation [deg]"
    PlotConf["xTicks"] = range(0, 91, 10)
    PlotConf["xLim"] = [0, 90]

    PlotConf["yLabel"] = "Iono Error [m]"

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "PRN"
    PlotConf["ColorBarTicks"] = sorted(
        np.unique(los_data[LOS_IDX["PRN"]][filter_condition]))

    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    PlotConf["xData"] = los_data[LOS_IDX["ELEV"]][filter_condition]
    PlotConf["yData"][Label] = los_data[LOS_IDX["IONO-ERR"]][filter_condition]
    PlotConf["zData"][Label] = los_data[LOS_IDX["PRN"]][filter_condition]

    return PlotConf
