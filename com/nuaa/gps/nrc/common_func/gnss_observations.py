from com.nuaa.gps.nrc.config.configuration import global_gpsConfig
import numpy as np


def gnss_observations(t, XS_tx, VS_tx, Xr_M, Vr_M, Xr_R, Vr_R):
    Xb0 = Xr_M - Xr_R;
    nsat = len(global_gpsConfig.sign_set_prnmat);
    R1_M = np.zeros((nsat, global_gpsConfig.sign_set_data_length))
    e1_M = np.zeros((nsat, global_gpsConfig.sign_set_data_length))
    pr1_M = R1_M
    for i in range(nsat):
        for j in range(global_gpsConfig.sign_set_data_length):
            R1_M[i, j] = np.sqrt((XS_tx[i, :].transfer() - Xr_M).transfer() * (XS_tx[i, :].transfer() - Xr_M))
            e1_M[i, :, j] = (XS_tx[i, :] - Xr_M.transfer()) / R1_M(i, j);
            pr1_M[i, j] = R1_M[ i, j] + global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_code * np.random.randint(
                1, 5);

            pass
