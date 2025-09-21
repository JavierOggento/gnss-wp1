import numpy as np
from COMMON.coordinates import xyz2llh
from COMMON.constants import LOS_IDX


def sat_tracks(nav_solution, los_data):
    if nav_solution == 'GPS':
        constellation_id = 'G'
        constellation = 'GPS'
    else:
        constellation_id = 'E'
        constellation = 'Galileo'

    plot_conf = {
        "Title": "Simulated GNSS " + constellation + " Satellite Tracks",
        "Type": "Lines",
        "Grid": True,
        "FigSize": (16.8, 15.2),
        "Marker": ".",
        "LineWidth": 0.1,
    }
    plot_conf["Map"] = True

    plot_conf["LonMin"] = -180
    plot_conf["LonMax"] = 180
    plot_conf["LatMin"] = -60
    plot_conf["LatMax"] = 60

    lonStep = 15
    latStep = 10
    plot_conf["LonStep"] = lonStep
    plot_conf["LatStep"] = latStep

    # plot_conf["yLabel"] = "Latitude [deg]"
    plot_conf["yTicks"] = range(
        plot_conf["LatMin"], plot_conf["LatMax"] + 1, lonStep)
    plot_conf["yLim"] = [plot_conf["LatMin"], plot_conf["LatMax"]]

    # plot_conf["xLabel"] = "Longitude [deg]"
    plot_conf["xTicks"] = range(
        plot_conf["LonMin"], plot_conf["LonMax"] + 1, latStep)
    plot_conf["xLim"] = [plot_conf["LonMin"], plot_conf["LonMax"]]

    plot_conf["ColorBar"] = "gnuplot"
    plot_conf["ColorBarLabel"] = "Elevation [deg]"
    plot_conf["ColorBarMin"] = 0.0
    plot_conf["ColorBarMax"] = 90.0

    # Transform ECEF to Geodetic

    filter_condition = los_data[LOS_IDX["CONST"]] == constellation_id

    sat_x = los_data[LOS_IDX["SAT-X"]][filter_condition]

    DataLen = len(sat_x)

    count = 0

    longitude = np.zeros(DataLen)
    latitude = np.zeros(DataLen)
    elev = np.zeros(DataLen)

    for index in sat_x.keys():
        x = sat_x[index]
        y = los_data[LOS_IDX["SAT-Y"]][index]
        z = los_data[LOS_IDX["SAT-Z"]][index]

        longitude[count], latitude[count], _ = xyz2llh(x, y, z)
        elev[count] = los_data[LOS_IDX["ELEV"]][index]

        count += 1

    plot_conf["xData"] = {}
    plot_conf["yData"] = {}
    plot_conf["zData"] = {}

    Label = 0
    plot_conf["xData"] = longitude
    plot_conf["yData"][Label] = latitude
    plot_conf["zData"][Label] = elev

    return plot_conf
