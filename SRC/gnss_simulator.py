
########################################################################
# gnss_simulator.py:
# This is the GNSS Simulation Module of SENFUS tool
#
#  Project:        SENFUS
#  File:           gnss_simulator.py
#
#   Author: GNSS Academy
#   Copyright 2025 GNSS Academy
########################################################################

# Import External functions/variables
# ----------------------------------------------------------------------


# Import Internal functions/variables
# ----------------------------------------------------------------------

# TODO: implement propagate_orbits()
def propagate_orbits(receiverInfo, sp3Data, apoData):
    #  Extract receiverTime from receiverState

    #  Define tolerance (1e-09) and maxIterations (10)

    #  Initialize satList, satsPos, and satsVel

    #  For each satelliteLabel in sp3Data:
    #   If satelliteLabel corresponds to a configured constellation:
    #
    #     Interpolate the position of the satellite from SP3 data at Reception Time
    #     Compute initial estimate of transmissionTime using Geometrical Range

    #     For a limited number of iterations:
    #         Interpolate satellite position at estimated transmissionTime
    #         Compute geometric range to receiver
    #         Update transmissionTime estimate
    #         If convergence condition is met:
    #             Exit loop
    #
    #     If satelliteLabel exists in apoData:
    #         Retrieve apoCorrection
    #     Else:
    #         Continue to next satellite
    #
    #     Apply apoCorrection to satellite position
    #
    #     Store corrected satellite position in satsPos
    #
    #     Estimate satellite velocity at transmissionTime
    #
    #     Store estimated velocity in satsVel

    #     Add satelliteLabel to satList
    #  Return satList, satsPos, satsVel
    pass


# TODO: implement simulate_gnss_corrected_meas()
def simulate_gnss_corrected_meas():
    # Extract timeDelta from receiverState as t-t0
    #
    # Extract receiver position in geodetics
    # Extract receiver position and velocity in ECEF
    #
    #  (satList, satsPosTt, satsVelTt) = propagateOrbits(receiverInfo, sp3Data, apoData)
    #  (gnssCorrectedMeas, losInfo) = generateGnssMeasurements(conf,
    #                                                         timeDelta,
    #                                                         satInfo,
    #                                                         receiverInfo)
    #  Return (gnssCorrectedMeas, losInfo)
    pass
