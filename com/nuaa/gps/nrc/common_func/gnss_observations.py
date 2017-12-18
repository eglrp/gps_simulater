from com.nuaa.gps.nrc.config.configuration import global_gpsConfig
import numpy as np
import math
import math.floor as floor

from com.nuaa.gps.nrc.entity.gnss_observations import GNSSObservations


def gnss_observations(t, XS_tx, VS_tx, Xr_M, Vr_M, Xr_R, Vr_R):
    Xb0 = Xr_M - Xr_R
    nsat = len(global_gpsConfig.sign_set_prnmat)
    gnss_obs = GNSSObservations(np.zeros((nsat, global_gpsConfig.sign_set_data_length)))
    for i in range(nsat):
        for j in range(global_gpsConfig.sign_set_data_length):
            gnss_obs.R1_M[i, j] = np.sqrt(
                (XS_tx[i, :].T - Xr_M).T * (XS_tx[i, :].T - Xr_M))
            gnss_obs.e1_M[i, :, j] = (XS_tx[i, :] - Xr_M.T) / gnss_obs.R1_M(i, j)
            gnss_obs.pr1_M[i, j] = gnss_obs.R1_M[
                                                i, j] + global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_code * np.random.randint()
            gnss_obs.pr2_M[i, j] = gnss_obs.pr1_M(i, j)

            gnss_obs.ph1_M[i, j] = gnss_obs.R1_M[i, j] / global_gpsConfig.sign_set_g_bo1 - math.floor[
                gnss_obs.R1_M[
                    i, j] / global_gpsConfig.sign_set_g_bo1] + global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_g_bo1 * np.randint()

            gnss_obs.N1_M[i, j] = floor[
                gnss_obs.R1_M[i, j] / global_gpsConfig.sign_set_g_bo1]  # L1频点 整周模糊度
            gnss_obs.N2_M[i, j] = floor[
                gnss_obs.R1_M[i, j] / global_gpsConfig.sign_set_g_bo2]  # L2频点 整周模糊度

            gnss_obs.dop1_M[i, j] = np.matrix([Vr_M.T, -VS_tx[i, :]]) * gnss_obs.e1_M[i, :,
                                                                                          j].T / global_gpsConfig.sign_set_c * global_gpsConfig.sign_set_l1  # 当卫星与接收机相对远离时，多普勒频移为负，载波相位测量值变大
            gnss_obs.dop2_M[i, j] = np.matrix([Vr_M.T - VS_tx[i, :]]) * gnss_obs.e1_M[i, :,
                                                                                          j].T / global_gpsConfig.sign_set_c * global_gpsConfig.sign_set_l2
            gnss_obs.e1_Mx[i, j] = [XS_tx[i, 1] - Xr_M[1]] / gnss_obs.R1_M[i, j]
            gnss_obs.e1_My[i, j] = [XS_tx[i, 2] - Xr_M[2]] / gnss_obs.R1_M[i, j]
            gnss_obs.e1_Mz[i, j] = [XS_tx[i, 3] - Xr_M[3]] / gnss_obs.R1_M[i, j]
            # Rover 的L1、L2、L5频点伪距、载波相位观测值和整周模糊度、多普勒频移观测量
            gnss_obs.R1_R[i, j] = np.sqrt[
                [XS_tx[i, :].T - Xr_R].T * [XS_tx[i, :].T - Xr_R]]  # 卫星到流动站的距离
            gnss_obs.e1_R[i, :, j] = [XS_tx[i, :] - Xr_R.T] / gnss_obs.R1_R[
                i, j]  # 卫星到流动站的方向余弦
            gnss_obs.pr1_R[i, j] = gnss_obs.R1_R[
                                                i, j] + global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_code * np.randint()  # 流动站伪距观测量    ???
            gnss_obs.pr2_R[i, j] = gnss_obs.pr1_R[i, j]

            gnss_obs.ph1_R[i, j] = gnss_obs.R1_R[i, j] / global_gpsConfig.sign_set_g_bo1 - floor[
                gnss_obs.R1_R[
                    i, j] / global_gpsConfig.sign_set_g_bo1] + global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_g_bo1 * np.randint()  # L1频点 载波相位观测值
            gnss_obs.ph2_R[i, j] = gnss_obs.R1_R[i, j] / global_gpsConfig.sign_set_l2 - floor[
                gnss_obs.R1_R[
                    i, j] / global_gpsConfig.sign_set_l2] + global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_l2 * np.randint()  # L2频点 载波相位观测值
            # ph5_R[i,j]  = R1_R[i,j]/sign_set.G_bo5 - floor[R1_R[i,j]/sign_set.G_bo5] + global_gpsConfig.sign_set_phase*sign_set.G_bo5*np.randint()    #L5频点 载波相位观测值

            gnss_obs.N1_R[i, j] = floor[
                gnss_obs.R1_R[i, j] / global_gpsConfig.sign_set_g_bo1]  # L1频点 整周模糊度
            gnss_obs.N2_R[i, j] = floor[
                gnss_obs.R1_R[i, j] / global_gpsConfig.sign_set_l2]  # L2频点 整周模糊度
            # N5_R[i,j]  =floor[R1_R[i,j]/sign_set.G_bo5]         #L5频点 整周模糊度

            gnss_obs.dop1_R[i, j] = np.matrix([Vr_R.T, -VS_tx[i, :]]) * gnss_obs.e1_R[i, :,
                                                                                          j].T / global_gpsConfig.sign_set_c * global_gpsConfig.sign_set_l1  # 当卫星与接收机相对远离时，多普勒频移为负，载波相位测量值变大
            gnss_obs.dop2_R[i, j] = np.matrix([Vr_R.T, -VS_tx[i, :]]) * gnss_obs.e1_R[i, :,
                                                                                          j].T / global_gpsConfig.sign_set_c * global_gpsConfig.sign_set_l2
            # dop5_R[i,j] = [Vr' -VS_tx[i,:]]*e1_R[i,j]/global_gpsConfig.sign_set_c*sign_set.L5
            gnss_obs.e1_Rx[i, j] = [XS_tx[i, 1] - Xr_R[1]] / gnss_obs.R1_R[i, j]
            gnss_obs.e1_Ry[i, j] = [XS_tx[i, 2] - Xr_R[2]] / gnss_obs.R1_R[i, j]
            gnss_obs.e1_Rz[i, j] = [XS_tx[i, 3] - Xr_R[3]] / gnss_obs.R1_R[i, j]
    return gnss_obs
