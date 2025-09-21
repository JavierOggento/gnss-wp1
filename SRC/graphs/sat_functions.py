# Copyright (C) GNSS ACADEMY
##
# Name          : SatFunctions.py
# Purpose       : Satellite Analyses functions
# Project       : WP0-JSNP
# Component     :
# Author        : GNSS Academy
# Creation date : 2021
# File Version  : 1.0
# Version date  :
##


from pandas import read_csv
from COMMON.plots import generatePlot
from COMMON.constants import LOS_IDX
from .sat_tracks import sat_tracks
from .geom_range import geom_range
from .range_rate import range_rate
from .sis_error import sis_error
from .iono_error import iono_error
from .tropo_error import tropo_error
from .local_error import local_error
from .range_error import range_error

BASE_PATH = 'SCN/SCEN-SENFUS/OUT/LOS/'


def getData(los_file, usecols):
    """
    Read LOS data from a file and return the specified columns.
    """
    los_data = read_csv(
        los_file,
        delim_whitespace=True,
        skiprows=1,
        header=None,
        usecols=usecols,
    )
    return los_data


def plot_satellite_tracks(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["CONST"],
            LOS_IDX["SAT-X"],
            LOS_IDX["SAT-Y"],
            LOS_IDX["SAT-Z"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot satellite Tracks ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(sat_tracks('GPS', los_data),
                     BASE_PATH + "GPSSatTracks.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(sat_tracks('GAL', los_data),
                     BASE_PATH + "GALSatTracks.png")


def plot_geometrical_range(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["CONST"],
            LOS_IDX["GEOM-RANGE"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot Geometrical Range ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(geom_range('GPS', los_data),
                     BASE_PATH + "GPSGeomRange.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(geom_range('GAL', los_data),
                     BASE_PATH + "GALGeomRange.png")


def plot_range_rate(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["CONST"],
            LOS_IDX["RANGE-RATE"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot Geometrical Range Rate ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(range_rate('GPS', los_data),
                     BASE_PATH + "GPSRangeRate.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(range_rate('GAL', los_data),
                     BASE_PATH + "GALRangeRate.png")


def plot_sis_errors(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["PRN"],
            LOS_IDX["CONST"],
            LOS_IDX["SIS-ERR"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot Sis Error ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(sis_error('GPS', los_data),
                     BASE_PATH + "GPSSisError.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(sis_error('GAL', los_data),
                     BASE_PATH + "GALSisError.png")


def plot_iono_errors(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["PRN"],
            LOS_IDX["CONST"],
            LOS_IDX["IONO-ERR"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot Iono Error ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(iono_error('GPS', los_data),
                     BASE_PATH + "GPSIonoError.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(iono_error('GAL', los_data),
                     BASE_PATH + "GALIonoError.png")


def plot_tropo_errors(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["PRN"],
            LOS_IDX["CONST"],
            LOS_IDX["TROPO-ERR"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot Tropo Error ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(tropo_error('GPS', los_data),
                     BASE_PATH + "GPSTropoError.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(tropo_error('GAL', los_data),
                     BASE_PATH + "GALTropoError.png")


def plot_local_errors(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["PRN"],
            LOS_IDX["CONST"],
            LOS_IDX["LOCAL-ERR"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot Local Error ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(local_error('GPS', los_data),
                     BASE_PATH + "GPSLocalError.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(local_error('GAL', los_data),
                     BASE_PATH + "GALLocalError.png")


def plot_range_errors(nav_solution, los_file):
    los_data = getData(
        los_file,
        [
            LOS_IDX["SOD"],
            LOS_IDX["PRN"],
            LOS_IDX["CONST"],
            LOS_IDX["RANGE-ERR"],
            LOS_IDX["ELEV"],
        ],
    )

    print("Plot Range Error ...")

    if (nav_solution == 'GPSGAL') or (nav_solution == 'GPS'):
        generatePlot(range_error('GPS', los_data),
                     BASE_PATH + "GPSRangeError.png")
    if (nav_solution == 'GPSGAL') or (nav_solution == 'GAL'):
        generatePlot(range_error('GAL', los_data),
                     BASE_PATH + "GALRangeError.png")
