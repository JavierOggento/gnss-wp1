
########################################################################
# input_output.py:
# This is the I/O Module of SENFUS tool
#
#  Project:        SENFUS
#  File:           input_output.py
#
#   Author: GNSS Academy
#   Copyright 2025 GNSS Academy
########################################################################

# Import External functions/variables
# ----------------------------------------------------------------------
import sys
import os
from collections import OrderedDict
import pandas as pd
import numpy as np

from COMMON.gnss_constants import EARTH_RADIUS, ECCENTRICITY

# Import Internal functions/variables
# ----------------------------------------------------------------------


def check_conf_param(key, fields, min_fields, max_fields, low_lim, upp_lim):
    values = []
    len_fields = len(fields) - 1

    if len_fields < min_fields:
        sys.stderr.write(
            f"ERROR: Too few fields ({len_fields}) for configuration parameter {key}. Minimum = {min_fields}\n")
        sys.exit(-1)

    if len_fields > max_fields:
        sys.stderr.write(
            f"ERROR: Too many fields ({len_fields}) for configuration parameter {key}. Maximum = {max_fields}\n")
        sys.exit(-1)

    for field in fields[1:]:
        try:
            values.append(float(field))
        except:
            try:
                check = field.isnumeric()
            except:
                check = field.isnumeric()

            if check:
                values.append(int(field))
            else:
                values.append(field)

    for i, field in enumerate(values):
        if isinstance(low_lim[i], (int, float)):
            try:
                if field < low_lim[i] or field > upp_lim[i]:
                    sys.stderr.write(
                        f"ERROR: Configuration parameter {key} {field} is out of range [{low_lim[i]}, {upp_lim[i]}]\n")
            except:
                sys.stderr.write(
                    f"ERROR: Wrong type for configuration parameter {key}\n")
                sys.exit(-1)

    return values[0] if len(values) == 1 else values


def read_conf(cfg_file):

    conf = OrderedDict()

    n_read_params = 0

    with open(cfg_file, 'r') as f:
        lines = f.readlines()

        for line in lines:
            if line[0] != '#':
                fields = line.rstrip('\n').split(' ')
                fields = list(filter(None, fields)) if '' in fields else fields

                if fields:
                    if len(fields) == 1:
                        sys.stderr.write(
                            f"ERROR: Configuration file contains a parameter with no value: {line}")
                        sys.exit(-1)

                    key = fields[0]

                    if key == 'SAMPLING_RATE':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [0], [86400])
                        n_read_params += 1

                    elif key == 'NAV_SOLUTION':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [None], [None])
                        n_read_params += 1

                    elif key == 'GNSS_OUT':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [0], [1])
                        n_read_params += 1

                    elif key == 'SEED':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [0.0], [np.inf])
                        n_read_params += 1

                    elif key == 'RCVR_TRAJ':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [None], [None])
                        n_read_params += 1

                    elif key == 'RCVR_MASK':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [0], [90])
                        n_read_params += 1

                    elif key == 'GNSS_RATE':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [0], [86400])
                        n_read_params += 1

                    elif key == 'SP3_FILE':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [None], [None])
                        n_read_params += 1

                    elif key == 'APO_FILE':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [None], [None])
                        n_read_params += 1

                    elif key == 'RCVR_CLK':
                        conf[key] = check_conf_param(
                            key, fields, 2, 2, [None, None], [None, None])
                        n_read_params += 1

                    elif key == 'RANGE_RATE_ERR_SIGMA':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [None], [None])
                        n_read_params += 1

                    elif key == 'ORBCLK_ERR_SIGMA':
                        conf[key] = check_conf_param(
                            key, fields, 2, 2, [None, None], [None, None])
                        n_read_params += 1

                    elif key == 'IONO_ERR_SIGMA':
                        conf[key] = check_conf_param(
                            key, fields, 2, 2, [None, None], [None, None])
                        n_read_params += 1

                    elif key == 'TROPO_ERR_SIGMA':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [None], [None])
                        n_read_params += 1

                    elif key == 'LOCAL_ERR_SIGMA':
                        conf[key] = check_conf_param(
                            key, fields, 2, 2, [None, None], [None, None])
                        n_read_params += 1

                    if key == 'GENERATE_LOS_POS_FILES':
                        conf[key] = check_conf_param(
                            key, fields, 1, 1, [0], [1])
                        n_read_params += 1

    conf['RCVR_MASK_RAD'] = np.deg2rad(conf['RCVR_MASK'])

    np.random.seed(int(conf['SEED']))

    return conf


def read_traj_data(file_path):
    """
    Reads a GNSS data file with fixed-width columns into a Pandas DataFrame.

    Args:
        file_path (str): Path to the input data file.

    Returns:
        pd.DataFrame: DataFrame containing the parsed GNSS data.
    """
    column_names = [
        "Time(s)", "Longitude(deg)", "Latitude(deg)", "Height(m)",
        "NorthVel(m/s)", "EastVel(m/s)", "DownVel(m/s)",
        "Roll(deg)", "Pitch(deg)", "Yaw(deg)"
    ]

    # Define column widths based on fixed format (adjust as needed)
    col_widths = [10, 20, 20, 13, 13, 13, 13, 11, 11, 11]

    # Read the file, skipping the first header row
    df = pd.read_fwf(file_path, widths=col_widths,
                     names=column_names, skiprows=1)

    # Add delta time column
    df["T-T0(s)"] = df["Time(s)"] - df.iloc[0]["Time(s)"]

    # Add ECEF coordinates
    # WGS-84 ellipsoid parameters
    a = 6378137.0  # Semi-major axis (m)
    e2 = 0.00669437999014  # Square of first eccentricity

    # Convert degrees to radians
    lat_rad = np.radians(df["Latitude(deg)"])
    lon_rad = np.radians(df["Longitude(deg)"])
    h = df["Height(m)"]

    # Compute N (prime vertical radius of curvature)
    N = a / np.sqrt(1 - e2 * np.sin(lat_rad) ** 2)

    # Compute ECEF coordinates
    df["X(m)"] = (N + h) * np.cos(lat_rad) * np.cos(lon_rad)
    df["Y(m)"] = (N + h) * np.cos(lat_rad) * np.sin(lon_rad)
    df["Z(m)"] = (N * (1 - e2) + h) * np.sin(lat_rad)

    return df


def read_sp3(file_path):
    """
    Reads an SP3 file and extracts the header, satellite list, and ephemeris data.

    Parameters:
        file_path (str): Path to the SP3 file.

    Returns:
        dict: Contains 'header', 'satellites', and 'ephemeris' as DataFrame.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header = []
    satellites = []
    ephemeris_data = []
    reading_ephemeris = False

    for line in lines:
        if line.startswith('#'):  # Header lines
            header.append(line.strip())
        elif line.startswith('+'):  # Satellite list
            satellites.extend(line[3:].strip().split())
        elif line.startswith('*'):  # Epoch time marker
            reading_ephemeris = True
            time_parts = list(map(float, line[2:].split()))
            # Convert to seconds of the day
            current_epoch = time_parts[3] * 3600 + \
                time_parts[4] * 60 + time_parts[5]
        elif reading_ephemeris and (line.startswith('P') or line.startswith('V')):
            values = line.split()
            prn = values[0][1:]
            # Only include GPS (G) and Galileo (E)
            if prn.startswith('G') or prn.startswith('E'):
                x, y, z, _ = map(float, values[1:])
                ephemeris_data.append(
                    [current_epoch, prn, x, y, z])

    positions_df = pd.DataFrame(ephemeris_data, columns=[
                                'SOD', 'PRN', 'X', 'Y', 'Z',]).groupby('PRN')

    return positions_df


def read_apo(file_path):
    """
    Reads the APO file.

    Parameters:
        file_path (str): Path to the APO file.

    Returns:
        df: Contains APO data.
    """
    column_names = [
        "Const", "PRN", "x_f1", "y_f1", "z_f1",
        "x_f2", "y_f2", "z_f2",
    ]

    # Define column widths based on fixed format (adjust as needed)
    col_widths = [14, 3, 11, 11, 11, 11, 11, 11]

    # Read the file, skipping the first header row
    df = pd.read_fwf(file_path, widths=col_widths,
                     names=column_names, skiprows=1)

    apo_data = {}
    for i, row in df.iterrows():
        apo_data[row["Const"] + "%02d" % row["PRN"]] = row["z_f1"] / 1000

    return apo_data


def create_output_file(path, hdr):
    """
    Purpose: Open output file and write its header

    Parameters
    ==========
    path: str
        Path to file
    hdr: str
        File header

    Returns
    =======
    f: File descriptor
        Descriptor of output file
    """

    # Display message
    print(f"INFO: Creating file: {path}...")

    # Create output directory, if needed
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))

    # Open output file
    f = open(path, 'w')

    # Write header
    f.write(hdr)

    return f

# End of create_output_file()
