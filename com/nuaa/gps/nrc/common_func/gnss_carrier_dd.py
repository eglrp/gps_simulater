from com.nuaa.gps.nrc.config.configuration import global_gpsConfig
from com.nuaa.gps.nrc.entity.gnss_carrier_dd import GNSS_Carrier_DD
import numpy as np


def gnss_carrier_dd(t, gnss_obs):
    nsat = len[global_gpsConfig.sign_set_prnmat];
    gnss_dd = GNSS_Carrier_DD(np.zeros((nsat, global_gpsConfig.sign_set_data_length)));

    #    产生双差观测量   ############
    for i in range(nsat):  # i代表哪颗卫星
        for j in range(global_gpsConfig.sign_set_data_length):  # j代表第几个采样点
            # 单差观测值     单差为站间作差
            gnss_dd.dph1[i, j] = gnss_obs.ph1_R[i, j] - gnss_obs.ph1_M[i, j];  # L1频点  载波相位单差观测值
            gnss_dd.dph2[i, j] = gnss_obs.ph2_R[i, j] - gnss_obs.ph2_M[i, j];  # L2频点  载波相位单差观测值
            # dph5[i,j]  =ph5_R[i,j]  - ph5_M[i,j];              #L5频点  载波相位单差观测值

            gnss_dd.dpr1[i, j] = gnss_obs.pr1_R[i, j] - gnss_obs.pr1_M[i, j];  # 站间单差伪距观测量
            gnss_dd.de1x[i, j] = gnss_obs.e1_Mx[i, j] + gnss_obs.e1_Rx[i, j];  # 计算方向余弦
            gnss_dd.de1y[i, j] = gnss_obs.e1_My[i, j] + gnss_obs.e1_Ry[i, j];
            gnss_dd.de1z[i, j] = gnss_obs.e1_Mz[i, j] + gnss_obs.e1_Rz[i, j];

            # 单差整周模糊度
            gnss_dd.dN1[i, j] = gnss_obs.N1_R[i, j] - gnss_obs.N1_M[i, j];  # L1频点 站间单差整周模糊度
            gnss_dd.dN2[i, j] = gnss_obs.N2_R[i, j] - gnss_obs.N2_M[i, j];  # L2频点 站间单差整周模糊度
            # dN5[i,j]   = N5_R[i,j]  - N5_M[i,j];              # L5频点 站间单差整周模糊度

    ##  求双差观测值       星间作差
    for k in range(nsat - 1):  # 双差为卫星作差, 以第1颗卫星为参考星
        ##  双差观测值
        gnss_dd.ddph1[k, :] = gnss_dd.dph1[k + 1, :] - gnss_dd.dph1[1, :];  # L1双差载波相位
        gnss_dd.ddph2[k, :] = gnss_dd.dph2[k + 1, :] - gnss_dd.dph2[1, :];  # L2
        # ddph5[k,:]  = dph5[k+1,:]  - dph5[1,:];            #L5

        gnss_dd.ddN1[k, :] = gnss_dd.dN1[k + 1, :] - gnss_dd.dN1[1, :];  # L1双差模糊度
        gnss_dd.ddN2[k, :] = gnss_dd.dN2[k + 1, :] - gnss_dd.dN2[1, :];  # L2双差模糊度
        # ddN5[k,:]   = dN5[k+1,:]   - dN5[1,:];  #L5双差模糊度

        gnss_dd.ddpr1[k, :] = gnss_dd.dpr1[k + 1, :] - gnss_dd.dpr1[1, :];  # 双差伪距观测量
        gnss_dd.dde1x[k, :] = -0.5 * [gnss_dd.de1x[k + 1, :] - gnss_dd.de1x[1, :]];  # 双差方向余弦
        gnss_dd.dde1y[k, :] = -0.5 * [gnss_dd.de1y[k + 1, :] - gnss_dd.de1y[1, :]];
        gnss_dd.dde1z[k, :] = -0.5 * [gnss_dd.de1z[k + 1, :] - gnss_dd.de1z[1, :]];

    return gnss_dd