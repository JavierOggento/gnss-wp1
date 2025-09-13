

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
from collections import OrderedDict


# Import Internal functions/variables
# ----------------------------------------------------------------------


# LOS file
# ----------------------------------------------------------------------
# Header
los_hdr = "\
#SOD C PRN   ELEV    AZIM    SAT-X          SAT-Y          SAT-Z         SAT-VX       SAT-VY       SAT-VZ    MEASUREMENT  GEOM-RANGE  RANGE-RATE  SIS-ERR  IONO-ERR TROPO-ERR LOCAL-ERR RANGE-ERR\n"

# Line format
los_fmt = "%05d %1s %2s %8.3f %8.3f %14.3f %14.3f %14.3f %10.3f %10.3f %10.3f\
 %14.3f %14.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f".split()

# File columns
los_idx = OrderedDict({})
los_idx["SOD"] = 0
los_idx["CONST"] = 1
los_idx["PRN"] = 2
los_idx["ELEV"] = 3
los_idx["AZIM"] = 4
los_idx["SAT-X"] = 5
los_idx["SAT-Y"] = 6
los_idx["SAT-Z"] = 7
los_idx["SAT-VX"] = 8
los_idx["SAT-VY"] = 9
los_idx["SAT-VZ"] = 10
los_idx["MEASUREMENT"] = 11
los_idx["GEOM-RANGE"] = 12
los_idx["RANGE-RATE"] = 13
los_idx["SIS-ERR"] = 14
los_idx["IONO-ERR"] = 15
los_idx["TROPO-ERR"] = 16
los_idx["LOCAL-ERR"] = 17
los_idx["RANGE-ERR"] = 18


def write_gnss_los(gnss_los_file, length_gnss_corrected_meas, los_info):
    for i in range(length_gnss_corrected_meas):
        line_data = los_info[i]
        line_str = los_fmt[0] % line_data[0]
        for j in range(1, len(los_fmt)):
            line_str += ' ' + los_fmt[j] % line_data[j]
        gnss_los_file.write(line_str + '\n')

    # gnss_los_file.flush()


def plot_gnss_los(gnss_los_file_path):
    # TODO [GRAPHICS]: implement plot_gnss_los()
    pass


# POS file
# ----------------------------------------------------------------------
# Header
pos_hdr = "\
#SOD    LONEST    LATEST     ALTEST       CLKEST      NSV   HPE     VPE      NPE     EPE     DPE     GDOP    PDOP    TDOP\n"

# Line format
pos_fmt = "%07.2f %10.5f %10.5f %10.3f %12.3f %3d \
            %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f".split()

# File columns
pos_idx = OrderedDict({})
pos_idx["SOD"] = 0
pos_idx["LONEST"] = 1
pos_idx["LATEST"] = 2
pos_idx["ALTEST"] = 3
pos_idx["CLKEST"] = 4
pos_idx["NSV"] = 5
pos_idx["HPE"] = 6
pos_idx["VPE"] = 7
pos_idx["NPE"] = 8
pos_idx["EPE"] = 9
pos_idx["DPE"] = 10
pos_idx["GDOP"] = 11
pos_idx["PDOP"] = 12
pos_idx["TDOP"] = 13


def write_gnss_pos(gnss_pos_file, sod, length_gnss_corrected_meas, pos_info):
    # TODO [OUTPUT]: implement write_gnss_pos()
    pass


def plot_gnss_pos(gnss_pos_file_path):
    # TODO [GRAPHICS]: implement plot_gnss_pos()
    pass
