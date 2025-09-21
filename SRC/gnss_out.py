

########################################################################
# gnss_out.py:
# This is the GNSS Outputs Module of SENFUS tool
#
#  Project:        SENFUS
#  File:           gnss_out.py
#
#   Author: GNSS Academy
#   Copyright 2025 GNSS Academy
########################################################################

# Import External functions/variables
# ----------------------------------------------------------------------


# Import Internal functions/variables
# ----------------------------------------------------------------------
from graphs.sat_functions import plot_satellite_tracks, plot_geometrical_range, plot_range_rate, plot_sis_errors, plot_iono_errors, plot_tropo_errors, plot_local_errors, plot_range_errors

# LOS file
# ----------------------------------------------------------------------
# Header
los_hdr = "\
#SOD C PRN   ELEV    AZIM    SAT-X          SAT-Y          SAT-Z         SAT-VX       SAT-VY       SAT-VZ    MEASUREMENT  GEOM-RANGE  RANGE-RATE  SIS-ERR  IONO-ERR TROPO-ERR LOCAL-ERR RANGE-ERR\n"

# Line format
los_fmt = "%05d %1s %2s %8.3f %8.3f %14.3f %14.3f %14.3f %10.3f %10.3f %10.3f\
 %14.3f %14.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f".split()


def write_gnss_los(gnss_los_file, length_gnss_corrected_meas, los_info):
    for i in range(length_gnss_corrected_meas):
        line_data = los_info[i]
        line_str = los_fmt[0] % line_data[0]
        for j in range(1, len(los_fmt)):
            line_str += ' ' + los_fmt[j] % line_data[j]
        gnss_los_file.write(line_str + '\n')

    # gnss_los_file.flush()


def plot_gnss_los(nav_solution, gnss_los_file_path):
    plot_satellite_tracks(nav_solution, gnss_los_file_path)
    plot_geometrical_range(nav_solution, gnss_los_file_path)
    plot_range_rate(nav_solution, gnss_los_file_path)
    plot_sis_errors(nav_solution, gnss_los_file_path)
    plot_iono_errors(nav_solution, gnss_los_file_path)
    plot_tropo_errors(nav_solution, gnss_los_file_path)
    plot_local_errors(nav_solution, gnss_los_file_path)
    plot_range_errors(nav_solution, gnss_los_file_path)


# POS file
# ----------------------------------------------------------------------
# Header
pos_hdr = "\
#SOD    LONEST    LATEST     ALTEST       CLKEST      NSV   HPE     VPE      NPE     EPE     DPE     GDOP    PDOP    TDOP\n"

# Line format
pos_fmt = "%07.2f %10.5f %10.5f %10.3f %12.3f %3d \
            %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f".split()


def write_gnss_pos(gnss_pos_file, sod, length_gnss_corrected_meas, pos_info):
    # TODO [OUTPUT]: implement write_gnss_pos()
    pass


def plot_gnss_pos(gnss_pos_file_path):
    # TODO [GRAPHICS]: implement plot_gnss_pos()
    pass
