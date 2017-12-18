from com.nuaa.gps.nrc.common_func.geo_ecef import geo_ccef
from com.nuaa.gps.nrc.common_func.sate_posivelo import sate_posivelo
import numpy as np
import math

from com.nuaa.gps.nrc.config.configuration import global_gpsConfig


def satellite_positions(gps_config, t, orbit):
    nsat = len(gps_config.sign_set_prnmat);
    time_m = np.zeros(nsat, 1);
    time_r = np.zeros(nsat, 1);
    xs_tx = np.zeros(nsat, 3);
    vs_tx = np.zeros(nsat, 3);
    # position in GEOREF : long (deg)  lati (deg)  alti (m)
    posi_m = gps_config.sign_set_init_posi_m
    posi_r = gps_config.sign_set_init_posi_r

    # Master
    trace_m_posi__geo = posi_m;
    trace_m_posi__ecef = geo_ccef(posi_m);

    # Rover
    tracer_posi_geo = posi_r;
    tracer_posi_ecef = geo_ccef(posi_r);

    Xr1 = tracer_posi_ecef;
    t_orbit = gps_config.sign_set_t_sate + t;
    time_r = ''

    for j in range(nsat):
        # 解算卫星位置速度
        Xs, Vs = sate_posivelo(t_orbit, orbit[j, :]);  # 信号接收时刻卫星位置
        Tp1 = np.sqrt((Xs - Xr1).transpose() * (Xs - Xr1)) / global_gpsConfig.sign_set_c;
        Tp_old1 = 0;
        while ((Tp1 - Tp_old1) > 1e-12):
            Xs, Vs = sate_posivelo(t_orbit - Tp1, orbit[j, :]);
            C = np.matrix(
                [math.cos(global_gpsConfig.sign_set_wiee * Tp1), math.sin(global_gpsConfig.sign_set_wie * Tp1), 0],
                [-math.sin(global_gpsConfig.sign_set_wie * Tp1), math.cos(global_gpsConfig.sign_set_wie * Tp1), 0],
                [0, 0, 1])
            Xs = C * Xs;
            Vs = C * Vs;
            Tp_old1 = Tp1;
            Tp1 = np.sqrt((Xs - Xr1).transpose() * (Xs - Xr1)) / global_gpsConfig.sign_set_c;
        time_m[j, 1] = Tp1;
        xs_tx[j, :], vs_tx[j, :] = sate_posivelo(t_orbit - time_m[j, 1], orbit[j, :])

    Xr2 = tracer_posi_ecef;
    for j in range(nsat):
        Xs, Vs = sate_posivelo(t_orbit, orbit[j, :]);
        Tp2 = np.sqrt((Xs - Xr2).transpose() * (Xs - Xr2)) / global_gpsConfig.sign_set_c;
        Tp_old2 = 0;
        while ((Tp2 - Tp_old2) > 1e-12):
            Xs, Vs = sate_posivelo(t_orbit - Tp2, orbit[j, :]);
            C = np.matrix(
                [math.cos(global_gpsConfig.sign_set_wie * Tp2), math.sin(global_gpsConfig.sign_set_wie * Tp2), 0],
                [-math.sin(global_gpsConfig.sign_set_wie * Tp2), math.cos(global_gpsConfig.sign_set_wie * Tp2), 0],
                [0, 0, 1])
            Xs = C * Xs;
            Vs = C * Vs;
            Tp_old2 = Tp2;
            Tp2 = np.sqrt((Xs - Xr2).transpose() * (Xs - Xr2)) / global_gpsConfig.sign_set_c;
    time_r[j, 1] = Tp2;
    xs_tx[j, :], vs_tx[j, :] = sate_posivelo(t_orbit - time_r(j, 1), orbit[j, :]);
    return xs_tx, vs_tx, time_m, time_r
