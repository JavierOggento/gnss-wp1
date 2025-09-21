#!/usr/bin/env python

########################################################################
# senfus.py:
# This is the Main Module of SENFUS tool
#
#  Project:        SENFUS
#  File:           senfus.py
#
#   Author: GNSS Academy
#   Copyright 2025 GNSS Academy
#
# Usage:
#   senfus.py $SCEN_PATH
########################################################################


# Import External functions/variables
# ----------------------------------------------------------------------
import sys
import os
import numpy as np

# Import Internal functions/variables
# ----------------------------------------------------------------------
from input_output import read_conf, read_traj_data, read_sp3, read_apo
from input_output import create_output_file
from gnss_simulator import simulate_gnss_corrected_meas
from geometry import euler_to_rotmat, ned_to_ecef
# from gnss_positioning import estimate_gnss_position_lsq
# from geometry import compute_ned_errors
from gnss_out import pos_hdr, los_hdr
# from gnss_out import write_gnss_pos, plot_gnss_pos
from gnss_out import write_gnss_los, plot_gnss_los

# ----------------------------------------------------------------------
# INTERNAL FUNCTIONS
# ----------------------------------------------------------------------


def display_usage():
    sys.stderr.write(
        "ERROR: Please provide path to SCENARIO as a unique argument\n")


#######################################################
# MAIN BODY
#######################################################

# Print header
print('------------------------------------')
print('--> RUNNING SENFUS:')
print('------------------------------------')

# Check InputOutput Arguments
if len(sys.argv) != 2:
    display_usage()
    sys.exit()

# Extract the arguments
scen = sys.argv[1]

# Read conf
###############
# Select the Configuration file name
cfg_file = scen + '/CFG/senfus.cfg'

# Read conf file
# start here..., remove the following line and add a line to read the conf file:
conf = read_conf(cfg_file)

# Read input files
#####################
# Read true trajectory
# ------------------------
# Get trajectory data file path

traj_path = scen + '/INP/RCVR/' + conf["RCVR_TRAJ"]
true_traj = read_traj_data(traj_path)


# Read satellite APOs
# ------------------------
# Get APO data file path
apo_path = scen + '/INP/ATX/' + conf["APO_FILE"]
apo_data = read_apo(apo_path)

# Read satellite orbits
# ------------------------
# Get SP3 data file path
sp3_path = scen + '/INP/SP3/' + conf["SP3_FILE"]
sp3_data = read_sp3(sp3_path)


# sys.exit()


gnss_los_file_path = scen + '/OUT/GNSS/los_file.dat'

# Loop over trajectory records
if conf['GENERATE_LOS_POS_FILES']:
    # Open output files
    #######################
    gnss_los_file = create_output_file(gnss_los_file_path, los_hdr)

    for index, receiver_state in true_traj.iterrows():
        sod = receiver_state["Time(s)"]
        roll = np.deg2rad(receiver_state["Roll(deg)"])
        pitch = np.deg2rad(receiver_state["Pitch(deg)"])
        yaw = np.deg2rad(receiver_state["Yaw(deg)"])

        # Transform Euler angles to Rotation Matrix
        true_R = euler_to_rotmat((roll, pitch, yaw))

        # Get ECEF position, velocity and attitude
        receiver_ecef = ned_to_ecef(receiver_state, true_R)

        if sod % conf["GNSS_RATE"] == 0.0:
            # Simulate GNSS corrected Measurements
            # ---------------------------------------
            gnss_corrected_meas, los_info = simulate_gnss_corrected_meas(
                conf, receiver_state, receiver_ecef, sp3_data, apo_data)

            # [CHALLENGE] Estimate the GNSS-only position
            # ---------------------------------------
            # gnss_pvt = estimate_gnss_position_lsq(gnss_corrected_meas)

            if conf['GNSS_OUT']:
                # output line of sight information
                write_gnss_los(gnss_los_file, len(
                    gnss_corrected_meas), los_info)

                # [CHALLENGE] Estimate the Position Errors
                # ---------------------------------------
                # gnss_position_errors, gnss_velocity_errors, gnss_pos_geodetics = compute_ned_errors(gnss_pvt, receiver_state)

                # [CHALLENGE] output GNSS Position and Position Errors
                # write_gnss_pos(gnss_pos_file, sod, len(gnss_corrected_meas), gnss_pos_geodetics, gnss_pvt, gnss_position_errors)

    gnss_los_file.close()

if conf['GNSS_OUT']:
    # plot LOS outputs
    plot_gnss_los(conf['NAV_SOLUTION'], gnss_los_file_path)

    # [CHALLENGE] plot GNSS Position and Position Errors
    # plot_gnss_pos(gnss_pos_file_path)


# Print header
print('------------------------------------')
print('--> RUNNING SENFUS:')
print('------------------------------------')
