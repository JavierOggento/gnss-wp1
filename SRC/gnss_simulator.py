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

import math
import sys
import numpy as np
import pandas as pd
from COMMON.gnss_constants import EARTH_RADIUS, SPEED_OF_LIGHT, EARTH_ROTATION_RATE
from geometry import ecef_to_ned, skew_symmetric


def estimate_velocity(current_sat_sp3_data, sod):
    # TODO: implement estimate_velocity()
    # Define time_step = 1s
    time_step = 1
    #  Compute posMinus by interpolating sp3Data at (sod - timeStep/2)
    pos_minus = interpolate_sp3(current_sat_sp3_data, sod - time_step/2)
    #  Compute posPlus by interpolating sp3Data at (sod + timeStep/2)
    pos_plus = interpolate_sp3(current_sat_sp3_data, sod + time_step/2)
    #  Estimate satellite_velocity as (pos_plus - pos_minus) / time_step
    satellite_velocity = (pos_plus - pos_minus) / time_step
    #  Return satellite_velocity
    return satellite_velocity


def interpolate_sp3(current_sat_ephemeris, receiver_time):
    SP3_RATE = 300
    # TODO: implement interpolate_sp3()
    #   Extract satellite positions and times from SP3 data

    if (current_sat_ephemeris.loc[current_sat_ephemeris['SOD'] == receiver_time].shape[0] > 0):
        return current_sat_ephemeris.loc[current_sat_ephemeris['SOD'] == receiver_time][['X', 'Y', 'Z']].values[0]

    sp3_rows_lowest = current_sat_ephemeris.loc[current_sat_ephemeris['SOD'] < receiver_time].nlargest(
        5, 'SOD')

    sp3_rows_highest = current_sat_ephemeris.loc[current_sat_ephemeris['SOD'] > receiver_time].nsmallest(
        5, 'SOD')

    sp3_rows = pd.concat([sp3_rows_lowest, sp3_rows_highest])
    # print("sp3_rows:", sp3_rows)

    #   Perform linear interpolation to find satellite position at receiver_time

    # implement 10-point Lagrange interpolation between closest 10 SP3 positions
    if len(sp3_rows) < 10:
        sys.stderr.write(
            f"ERROR: Not enough SP3 data points for interpolation at {receiver_time}\n")
        sys.exit(-1)

    # Perform Lagrange interpolation
    # Prepare arrays for interpolation
    sod_points = sp3_rows['SOD'].values
    x_points = sp3_rows['X'].values
    y_points = sp3_rows['Y'].values
    z_points = sp3_rows['Z'].values

    def lagrange_interpolate(points, sod_points, receiver_time):
        result = 0.0
        n = len(points)
        for i in range(n):
            L = 1.0
            for j in range(n):
                if i != j:
                    L *= (receiver_time -
                          sod_points[j]) / (sod_points[i] - sod_points[j])
            result += L * points[i]
        return result

    interpolated_x = lagrange_interpolate(x_points, sod_points, receiver_time)
    interpolated_y = lagrange_interpolate(y_points, sod_points, receiver_time)
    interpolated_z = lagrange_interpolate(z_points, sod_points, receiver_time)

    interpolated_position = np.array(
        [interpolated_x, interpolated_y, interpolated_z])

    return interpolated_position


def propagate_orbits(conf, receiver_info, sp3_data, apo_data):
    receiver_ecef, receiver_state = receiver_info

    #  Extract receiverTime from receiverState
    receiver_time = receiver_state["Time(s)"]

    #  Define tolerance (1e-09) and maxIterations (10)
    tolerance = 1e-09
    max_iterations = 10

    #  Initialize satList, satsPos, and satsVel
    sats_list = []
    sats_pos = []
    sats_vel = []

    receiver_p = receiver_ecef[0]
    # print("receiver_p:", receiver_p)

    # satellites = sp3_data["satellites"]
    # ephemeris = sp3_data["ephemeris"]
    satellites = sp3_data.groups.keys()
    # ephemeris = sp3_data
    if conf["NAV_SOLUTION"] == "GPS":
        satellites = [k for k in satellites if 'G' in k]
    elif conf["NAV_SOLUTION"] == "GAL":
        satellites = [k for k in satellites if 'E' in k]

    for satellite_label in satellites:
        current_sat_ephemeris = sp3_data.get_group(satellite_label)
        # print("")
        # print("satellite_label:", satellite_label)
        estimated_transmission_time = receiver_time

        # Interpolate the position of the satellite from SP3 data at Reception Time
        # Compute initial estimate of transmissionTime using Geometrical Range

        for index in range(max_iterations):
            # print("Iteration:", index)
            # print("Estimated transmission_time:", estimated_transmission_time)

            # Interpolate satellite position at estimated transmissionTime
            current_interpolated_position = interpolate_sp3(
                current_sat_ephemeris, estimated_transmission_time) * 1000
            # current_sat_ephemeris, estimated_transmission_time)

            # print("current_position:", current_interpolated_position)

            current_geometric_range = np.linalg.norm(
                current_interpolated_position - receiver_p)

            # print("current_geometric_range:", current_geometric_range)

            current_transmission_time = receiver_time - \
                (current_geometric_range / SPEED_OF_LIGHT)

            # print("interpolated transmission_time:", current_transmission_time)

            # If convergence condition is met:
            if abs(current_transmission_time - estimated_transmission_time) < tolerance:
                estimated_transmission_time = current_transmission_time
                break

            estimated_transmission_time = current_transmission_time

        if not satellite_label in apo_data:
            continue

        #     If satelliteLabel exists in apoData:
        #         Retrieve apoCorrection
        #     Apply apoCorrection to satellite position
        k = current_interpolated_position / \
            np.linalg.norm(current_interpolated_position)
        # print("increment: ", apo_data[satellite_label] * k)

        apo_correction = apo_data[satellite_label] * k

        # final_estimated_position = current_interpolated_position + apo_correction
        final_estimated_position = current_interpolated_position - apo_correction

        # print("final_estimated_position:", final_estimated_position)

        #     Store corrected satellite position in satsPos
        sats_pos.append(final_estimated_position)

        #     Estimate satellite velocity at transmissionTime
        estimated_velocity = estimate_velocity(
            current_sat_ephemeris, estimated_transmission_time)

        # print("estimated_velocity:", estimated_velocity)

        sats_vel.append(estimated_velocity)

        sats_list.append(satellite_label)
    return sats_list, sats_pos, sats_vel


def iono_mapping(vertical_error_meters, elevation_degrees, iono_height_m=350000.0):
    """
    Transforms a vertical ionospheric error into a slant ionospheric error
    using the thin-shell ionospheric model.

    This function is based on the assumption that all free electrons in the
    ionosphere are concentrated in a single, thin shell at a fixed height.

    Args:
        vertical_error_meters (float): The vertical ionospheric error in meters.
        elevation_degrees (float): The satellite elevation angle in degrees.
        iono_height_m (float, optional): The effective height of the
                                          ionospheric shell in meters.
                                          Defaults to 350000.0 m, a commonly
                                          used value.

    Returns:
        float: The calculated slant ionospheric error in meters. Returns 0 if
               the elevation angle is too low for a valid calculation.
    """

    # Convert elevation angle from degrees to radians for trigonometric functions
    elevation_rad = math.radians(elevation_degrees)

    # Calculate the argument for the mapping function based on the provided
    # reference image's formula: m(elev) = [1 - (RE / (RE + h) * cos(elev))^2]^(-1/2)
    cos_elevation_term = math.cos(
        elevation_rad) * (EARTH_RADIUS / (EARTH_RADIUS + iono_height_m))

    # Check for valid mathematical domain to prevent math domain errors.
    if cos_elevation_term > 1.0:
        cos_elevation_term = 1.0
    elif cos_elevation_term < -1.0:
        cos_elevation_term = -1.0

    # Calculate the mapping function (MF) or obliquity factor
    mapping_function = 1.0 / math.sqrt(1.0 - cos_elevation_term**2)

    # Calculate the slant ionospheric error
    slant_error_meters = vertical_error_meters * mapping_function

    return slant_error_meters


def tropo_mapping(zenith_delay_meters, elevation_degrees):
    """
    Transforms a zenith tropospheric delay to a slant tropospheric delay
    using a simple geometric mapping function for pseudorange.

    This model is a more refined geometric model compared to the simple
    cosecant function.

    Args:
        zenith_delay_meters (float): The zenith tropospheric delay in meters.
        elevation_degrees (float): The satellite elevation angle in degrees.

    Returns:
        float: The calculated slant tropospheric delay in meters. Returns 0 if
               the elevation angle is not valid.
    """
    # Convert elevation angle from degrees to radians
    elevation_rad = math.radians(elevation_degrees)

    # Handle cases where the satellite is below the horizon
    if elevation_degrees <= 0:
        print("Warning: Elevation angle must be greater than 0.")
        return 0.0

    # Calculate the mapping function (m_pp) based on the image's formula
    # m_pp = 1.001 / sqrt(0.002001 + sin^2(E))
    sin_sq_elevation = math.sin(elevation_rad)**2
    mapping_function = 1.001 / math.sqrt(0.002001 + sin_sq_elevation)

    # Calculate the slant tropospheric delay
    slant_delay_meters = zenith_delay_meters * mapping_function

    return slant_delay_meters


def mops_mp_sigma(elevation_degrees):
    """
    Calculates the standard deviation of the multipath error (sigma) based on
    the MOPS multipath model.


    Args:
        elevation_degrees (float): The satellite elevation angle in degrees ($\theta$).

    Returns:
        float: The standard deviation of the multipath error in meters.
    """
    # Handle cases where the elevation is not valid for the model
    if elevation_degrees < 2:
        return 0.0

    # Calculate the sigma multipath error using the MOPS formula
    sigma = 0.13 + 0.53 * math.exp(-elevation_degrees / 10.0)

    return sigma


def simulate_ranging_errors(conf, satellite, elevation):
    #  Determine Satellite constellation

    if (satellite[0] == 'G'):
        sat_sigma_index = 0
    else:
        sat_sigma_index = 1

    # Initialize rangingError array
    ranging_errors = []

    # Draw sample for SIS error using orbital clock error model
    orbclk_err_sigma = conf['ORBCLK_ERR_SIGMA'][sat_sigma_index]
    sis_error = float(orbclk_err_sigma) * np.random.normal(0, 1)
    # print("sis_error:", sis_error)
    ranging_errors.append(sis_error)

    # Draw sample for ionospheric error and transform to slant using iono mapping function
    iono_error_sigma = conf['IONO_ERR_SIGMA'][sat_sigma_index]
    iono_error = float(iono_error_sigma) * np.random.normal(0, 1)
    slant_iono_error = iono_mapping(iono_error, elevation)
    # print("slant_iono_error:", slant_iono_error)
    ranging_errors.append(slant_iono_error)

    # Draw sample for tropospheric error and transform to slant using tropo mapping function
    tropo_error_sigma = conf['TROPO_ERR_SIGMA']
    tropo_error = float(tropo_error_sigma) * np.random.normal(0, 1)
    slant_tropo_error = tropo_mapping(tropo_error, elevation)
    # print("slant_tropo_error:", slant_tropo_error)
    ranging_errors.append(slant_tropo_error)

    # Draw sample for local error and add multipath contribution with MOPS model
    local_err_sigma = conf['LOCAL_ERR_SIGMA'][sat_sigma_index]
    local_err = (float(local_err_sigma) + mops_mp_sigma(elevation)
                 ) * np.random.normal(0, 1)
    # print("local_err:", local_err)

    ranging_errors.append(local_err)
    # Compute total range error as sum of individual errors
    # NOTE1: to draw samples from N(0, s), get the s first and multiply it by a sample
    # drawn from N(0,1)
    # NOTE2: make sure the configured seed is taken into account to ensure repeatability in
    # the results
    total_ranging_error = sum(ranging_errors)
    # Return rangingError
    # print("total_ranging_error:", total_ranging_error)
    return sis_error, slant_iono_error, slant_tropo_error, local_err, total_ranging_error


def calculate_elevation(ned_pos):
    n = ned_pos[0]
    e = ned_pos[1]
    d = ned_pos[2]  # Assuming 'd' is positive downwards

    # Calculate the magnitude of the horizontal components
    horizontal_magnitude = np.sqrt(n**2 + e**2)

    if horizontal_magnitude == 0:
        # Satellite is directly overhead (or directly below if d > 0)
        if d < 0:  # Above horizon
            elevation_rad = np.pi / 2
        else:  # Below horizon
            elevation_rad = -np.pi / 2
    else:
        # Use atan2 for robust calculation. The argument order is atan2(y, x)
        # In our case, y is the vertical component (-d), and x is the horizontal magnitude.
        elevation_rad = np.arctan2(-d, horizontal_magnitude)

    # Convert radians to degrees
    elevation_deg = np.degrees(elevation_rad)
    return elevation_deg


def get_rotation_matrix(r_sat_ecef, r_rec_ecef):

    # Compute time of signal flight (approximate)
    delta_r = r_sat_ecef - r_rec_ecef
    rho_geom = np.linalg.norm(delta_r)
    tau = rho_geom / SPEED_OF_LIGHT

    # Earth's rotation vector and skew-symmetric matrix
    e_r = skew_symmetric(np.array([0, 0, EARTH_ROTATION_RATE]))

    # Earth rotation matrix using matrix exponential
    e_r_tau = e_r * tau
    return e_r_tau


def generateGnssMeasurements(conf, time_delta, sats_info, receiver_info):
    receiver_ecef, receiver_state = receiver_info
    # TODO: implement generateGnssMeasurements()
    #   Get constants: speedOfLight and earthRotationRate

    #   Initialize number of satellites and measurement arrays
    los_info = []
    gnss_measurements = []

    #   Compute ECEF to NED (North, East, Down) transformation matrix
    C_e_n = ecef_to_ned(receiver_state)

    #   For each satellite in satList:
    sats_list, sats_pos, sats_vel = sats_info

    for i, satellite in enumerate(sats_list):

        #       Compute line-of-sight vector and geometric range
        pos = sats_pos[i]
        vel = sats_vel[i] * 1000

        los_vector = pos - receiver_ecef[0]
        ned_pos = C_e_n @ (los_vector)

        elevation = calculate_elevation(ned_pos)
        azimuth = np.degrees(np.arctan2(ned_pos[1], ned_pos[0]))

        if elevation > conf['RCVR_MASK']:

            pos_with_sagnac_correction = pos - \
                get_rotation_matrix(pos, receiver_ecef[0]) @ pos

            vel_with_sagnac_correction = vel - \
                get_rotation_matrix(pos, receiver_ecef[0]) @ vel

            geometric_range = np.linalg.norm(
                pos_with_sagnac_correction - receiver_ecef[0])

            sis_error, slant_iono_error, slant_tropo_error, local_err, total_ranging_error = simulate_ranging_errors(
                conf, satellite, elevation)

            pseudo_measurement = geometric_range + \
                conf["RCVR_CLK"][0] + total_ranging_error

            range_rate = np.dot((vel_with_sagnac_correction - receiver_ecef[1]),
                                (pos_with_sagnac_correction - receiver_ecef[0])) / geometric_range


            drift = conf["RCVR_CLK"][1]

            # Draw sample for range rate error using range rate error model
            range_rate_sigma_error = conf['RANGE_RATE_ERR_SIGMA'] * \
                np.random.normal(0, 1)

            range_rate_error = drift + range_rate_sigma_error

            pseudo_range_rate = range_rate + range_rate_error

            gnss_measurements.append([
                pseudo_measurement,
                pseudo_range_rate,
                pos_with_sagnac_correction,
                vel_with_sagnac_correction,
            ])

            #           Store line-of-sight information (LOS) in output dataset
            los_info.append([
                time_delta,
                satellite[0],
                satellite[1:],
                elevation,
                azimuth,
                pos_with_sagnac_correction[0],
                pos_with_sagnac_correction[1],
                pos_with_sagnac_correction[2],
                vel_with_sagnac_correction[0],
                vel_with_sagnac_correction[1],
                vel_with_sagnac_correction[2],
                pseudo_measurement,
                geometric_range / 1000,
                pseudo_range_rate,
                sis_error,
                slant_iono_error,
                slant_tropo_error,
                local_err,
                total_ranging_error
            ])
            #           Store GNSS measurement

    #  Return gnssCorrectedMeas, losInfo

    return gnss_measurements, los_info


def simulate_gnss_corrected_meas(conf, receiver_state, receiver_ecef, sp3_data, apo_data):

    # TODO: implement simulate_gnss_corrected_meas()
    # Extract timeDelta from receiverState as t-t0
    receiver_time = receiver_state["Time(s)"]

    # Extract receiver position in geodetics
    # latitude, longitude, height = receiver_state["Longitude(deg)"], receiver_state[
    #     "Latitude(deg)"], receiver_state["Height(m)"]

    receiver_info = [receiver_ecef, receiver_state]

    sat_info = propagate_orbits(conf, receiver_info, sp3_data, apo_data)

    # time_#delta = 0  # receiver_state["TimeDelta(s)"]
    gnss_corrected_meas, los_info = generateGnssMeasurements(conf,
                                                             receiver_time,
                                                             sat_info,
                                                             receiver_info)

    return gnss_corrected_meas, los_info
