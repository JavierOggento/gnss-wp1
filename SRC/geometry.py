
########################################################################
# geometry.py:
# This is the Geometry Module of SENFUS tool
#
#  Project:        SENFUS
#  File:           geometry.py
#
#   Author: GNSS Academy
#   Copyright 2025 GNSS Academy
########################################################################

# Import External functions/variables
# ----------------------------------------------------------------------
import numpy as np

# Import Internal functions/variables
# ----------------------------------------------------------------------
from COMMON.GnssConstants import ECCENTRICITY, EARTH_RADIUS, IONO_HEIGHT


def skew_symmetric(vec):
    return np.array([[0, -vec[2], vec[1]],
                     [vec[2], 0, -vec[0]],
                     [-vec[1], vec[0], 0]])


def euler_to_rotmat(eul):
    # Precalculate sines and cosines of the Euler angles
    sin_phi = np.sin(eul[0])
    cos_phi = np.cos(eul[0])
    sin_theta = np.sin(eul[1])
    cos_theta = np.cos(eul[1])
    sin_psi = np.sin(eul[2])
    cos_psi = np.cos(eul[2])

    # Calculate coordinate transformation matrix
    C = np.zeros((3, 3))
    C[0, 0] = cos_theta * cos_psi
    C[0, 1] = cos_theta * sin_psi
    C[0, 2] = -sin_theta
    C[1, 0] = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi
    C[1, 1] = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi
    C[1, 2] = sin_phi * cos_theta
    C[2, 0] = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi
    C[2, 1] = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi
    C[2, 2] = cos_phi * cos_theta

    return C


def ecef_to_ned(receiver_state):
    """
    Compute the ECEF to NED (North, East, Down) transformation matrix.
    Args:
        receiver_state (dict): Receiver state containing position and orientation information
    Returns:
        C_e_n (np.ndarray): 3x3 transformation matrix from ECEF to NED
    """
    longitude = np.deg2rad(receiver_state["Longitude(deg)"])
    latitude = np.deg2rad(receiver_state["Latitude(deg)"])

    sin_lat = np.sin(latitude)
    cos_lat = np.cos(latitude)
    sin_long = np.sin(longitude)
    cos_long = np.cos(longitude)

    # ECEF to NED transformation matrix
    C_e_n = np.array([
        [-sin_lat * cos_long, -sin_lat * sin_long, cos_lat],
        [-sin_long,           cos_long,           0],
        [-cos_lat * cos_long, -cos_lat * sin_long, -sin_lat]
    ])

    return C_e_n


def ned_to_ecef(receiver_state, true_R):

    # Parameters

    longitude = np.deg2rad(receiver_state["Longitude(deg)"])
    latitude = np.deg2rad(receiver_state["Latitude(deg)"])
    height = receiver_state["Height(m)"]
    north_vel = receiver_state["NorthVel(m/s)"]
    east_vel = receiver_state["EastVel(m/s)"]
    down_vel = receiver_state["DownVel(m/s)"]

    # v: indica que es un vector de velocidad.
    # eb: significa "Earth to Body" (de la Tierra al cuerpo), pero en este contexto específico, se refiere al movimiento del cuerpo (por ejemplo, un receptor GPS) respecto a la Tierra.
    # n: indica que el vector está expresado en el sistema de coordenadas NED (North-East-Down).
    v_eb_n = np.array((north_vel, east_vel, down_vel))

    # Calculate transverse radius of curvature
    R_E = EARTH_RADIUS / np.sqrt(1 - (ECCENTRICITY * np.sin(latitude))**2)

    # Convert position
    cos_lat = np.cos(latitude)
    sin_lat = np.sin(latitude)
    cos_long = np.cos(longitude)
    sin_long = np.sin(longitude)

    r_eb_e = np.array([
        (R_E + height) * cos_lat * cos_long,
        (R_E + height) * cos_lat * sin_long,
        ((1 - ECCENTRICITY**2) * R_E + height) * sin_lat
    ])

    # Calculate ECEF to NED coordinate transformation matrix
    C_e_n = np.array([
        [-sin_lat * cos_long, -sin_lat * sin_long, cos_lat],
        [-sin_long, cos_long, 0],
        [-cos_lat * cos_long, -cos_lat * sin_long, -sin_lat]
    ])

    # Transform velocity
    v_eb_e = np.dot(C_e_n.T, v_eb_n)

    # Transform attitude
    C_b_e = np.dot(C_e_n.T, true_R)

    return r_eb_e, v_eb_e, C_b_e


def pv_ecef_to_geodetics(r_eb_e, v_eb_e):
    # Parameters

    # Convert position using Borkowski closed-form exact solution
    longitude = np.arctan2(r_eb_e[1], r_eb_e[0])

    # Calculate auxiliary values for Borkowski's solution
    k1 = np.sqrt(1 - ECCENTRICITY**2) * abs(r_eb_e[2])
    k2 = ECCENTRICITY**2 * EARTH_RADIUS
    beta = np.sqrt(r_eb_e[0]**2 + r_eb_e[1]**2)
    E = (k1 - k2) / beta
    F = (k1 + k2) / beta

    # Calculate Borkowski parameters
    P = 4/3 * (E * F + 1)
    Q = 2 * (E**2 - F**2)
    D = P**3 + Q**2

    V = (np.sqrt(D) - Q)**(1/3) - (np.sqrt(D) + Q)**(1/3)

    G = 0.5 * (np.sqrt(E**2 + V) + E)
    T = np.sqrt(G**2 + (F - V * G) / (2 * G - E)) - G

    latitude = np.sign(r_eb_e[2]) * \
        np.arctan((1 - T**2) / (2 * T * np.sqrt(1 - ECCENTRICITY**2)))

    height = (beta - EARTH_RADIUS * T) * np.cos(latitude) + \
        (r_eb_e[2] - np.sign(r_eb_e[2]) * EARTH_RADIUS *
         np.sqrt(1 - ECCENTRICITY**2)) * np.sin(latitude)

    # Calculate ECEF to NED coordinate transformation matrix
    cos_lat = np.cos(latitude)
    sin_lat = np.sin(latitude)
    cos_long = np.cos(longitude)
    sin_long = np.sin(longitude)

    C_e_n = np.array([
        [-sin_lat * cos_long, -sin_lat * sin_long, cos_lat],
        [-sin_long, cos_long, 0],
        [-cos_lat * cos_long, -cos_lat * sin_long, -sin_lat]
    ])

    # Transform velocity from ECEF to NED
    v_eb_n = np.dot(C_e_n, v_eb_e)

    return latitude, longitude, height, v_eb_n


# def geodetics_to_pv_ecef(latitude, longitude, height, v_eb_n):
    # TODO: check if it's correct
    # Convert geodetic coordinates to ECEF
    latitude = np.deg2rad(latitude)
    longitude = np.deg2rad(longitude)

    # Compute ECEF position
    N = EARTH_RADIUS / np.sqrt(1 - (ECCENTRICITY * np.sin(latitude))**2)
    x = (N + height) * np.cos(latitude) * np.cos(longitude)
    y = (N + height) * np.cos(latitude) * np.sin(longitude)
    z = (N * (1 - ECCENTRICITY**2) + height) * np.sin(latitude)

    r_eb_e = np.array([x, y, z])

    # Compute ECEF velocity
    R_N, R_E = radio_of_curvature(latitude)
    v_eb_e = np.array([
        v_eb_n[0] * (R_N + height),
        v_eb_n[1] * (R_E + height) * np.cos(latitude),
        -v_eb_n[2]
    ])

    return r_eb_e, v_eb_e


def radio_of_curvature(latitude):

    R_N = EARTH_RADIUS * (1 - ECCENTRICITY**2) / (1 - (ECCENTRICITY * np.sin(latitude))
                                                  ** 2)**1.5  # Meridian radius of curvature
    # Transverse radius of curvature
    R_E = EARTH_RADIUS / np.sqrt(1 - (ECCENTRICITY * np.sin(latitude))**2)

    return R_N, R_E


def compute_ned_errors(gnss_pvt, true_traj):
    true_longitude = np.deg2rad(true_traj["Longitude(deg)"])
    true_latitude = np.deg2rad(true_traj["Latitude(deg)"])
    true_height = true_traj["Height(m)"]

    north_vel = true_traj["NorthVel(m/s)"]
    east_vel = true_traj["EastVel(m/s)"]
    down_vel = true_traj["DownVel(m/s)"]

    receiver_vel_ned = np.array((north_vel, east_vel, down_vel))

    gnss_p = gnss_pvt[0]
    gnss_v = gnss_pvt[1]

    est_latitude, est_longitude, est_height, est_v_eb_ned = pv_ecef_to_geodetics(
        gnss_p, gnss_v)

    gnss_pos_geodetics = np.array(
        (np.rad2deg(est_latitude), np.rad2deg(est_longitude), est_height))

    # Position error calculation
    R_N, R_E = radio_of_curvature(true_latitude)
    delta_r_eb_n = np.zeros(3)
    delta_r_eb_n[0] = (est_latitude - true_latitude) * (R_N + true_height)
    delta_r_eb_n[1] = (est_longitude - true_longitude) * \
        (R_E + true_height) * np.cos(true_latitude)
    delta_r_eb_n[2] = -(est_height - true_height)

    # Velocity error calculation
    delta_v_eb_n = est_v_eb_ned - receiver_vel_ned

    return delta_r_eb_n, delta_v_eb_n, gnss_pos_geodetics


def compute_tropo_mpp_rad(Elev):
    if Elev >= np.deg2rad(4):
        TropoMpp = (1.001/(np.sqrt(0.002001+(np.sin(Elev))**2)))
    elif Elev >= np.deg2rad(2):
        TropoMpp = (1.001/(np.sqrt(0.002001+(np.sin(Elev))**2))) *\
            (1+0.015*(max(0, 4-Elev))**2)
    else:
        TropoMpp = np.nan

    return TropoMpp


def compute_iono_mpp_rad(ElevRad):

    Fpp = (1.0-((EARTH_RADIUS * np.cos(ElevRad)) /
                (EARTH_RADIUS + IONO_HEIGHT))**2)**(-0.5)

    return Fpp
