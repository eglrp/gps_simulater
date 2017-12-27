import numpy as np

from com.nuaa.gps.nrc.config.configuration import global_gpsConfig


def LSbaseline_ErrorModel(var_gnss_carrier_dd, global_transfer_Config, var_gnss_observations):

    ddph1 = var_gnss_carrier_dd.ddph1
    ddN1 = var_gnss_carrier_dd.ddN1
    ddph2 = var_gnss_carrier_dd.ddph2
    ddN2 = var_gnss_carrier_dd.ddN2
    Bk3 = global_transfer_Config.Bk3
    Xb0 = var_gnss_observations.Xb0

    delta_xb = np.zeros((3, 1));  # Xb=[delta_x  delta_y  delta_z]
    Xb = np.zeros([3, global_gpsConfig.sign_set_data_length]);
    delta_ddN = np.zeros((len(global_gpsConfig.sign_set_prnmat), 1));

    ##  最小二乘求解基线矢量
    for i in range(global_gpsConfig.sign_set_data_length):
        # 最小二乘法求解基线矢量
        ks = len(global_gpsConfig.sign_set_prnmat)

        D1 = 2 * [np.ones((ks - 1, ks - 1)) + np.eye[ks - 1]] * np.square([
            global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_g_bo1]);  # 求权矩阵  GPS-王惠南 P161页  公式（6.104）
        P1 = np.inv[D1];  # 求协方差矩阵 最小二乘中的权重值

        Af = Bk3[:, :, i];  # 载波的矩阵方程基线改正数系数A      [A B]*[delta_Xb  delta_ddN]'=L
        Bf = global_gpsConfig.sign_set_g_bo1 * np.eye([ks - 1]);  # 双差模糊度改正数系数B
        Lf = global_gpsConfig.sign_set_g_bo1 * [ddph1[:, i] + ddN1[:, i] + Af * Xb[:, i]];

        N11 = Af.transfer() * P1 * Af;
        N12 = Af.transfer() * P1 * Bf;
        N21 = Bf.transfer() * P1 * Af;
        N22 = Bf.transfer() * P1 * Bf;
        U1 = Af.transfer() * P1 * Lf;
        U2 = Bf.transfer() * P1 * Lf;
        N_0 = N21 * [N11].I * N12;
        P_Q = N22 - N_0;
        Segma_A = [P_Q].I;
        delta_ddN = Segma_A * [U2 - N21 * [N11].I * U1];  # 模糊度改正数
        Q_ddn = Segma_A;  # 模糊度协方差矩阵
        delta_xb = [N11].I * [U1 - N12 * delta_ddN];  # 基线改正数
        Q_xb = [N11].I + [N11].I * N12 * Segma_A * N21 * [N11].I;  # 基线协方差矩阵
        Q_xn = -[N11].I * N12 * Segma_A;  # 基线-模糊度协方差矩阵
        Q_nx = -Segma_A * N21 * [N11].I;
        Xb[:, i] = Xb0 + delta_xb;  # 基线向量
        ddN1[:, i] = ddN1[:, i] + delta_ddN;  # 模糊度浮点解

        ## 利用LAMBDA算法固定模糊度


        [afixed, sqnorm, Ps, Qzhat, Z, nfixed, mu] = LAMBDA[delta_ddN, Q_ddn, 4];  # ？? 单历元固定

        # 再求基线向量改正数
        delta_xbf = delta_xb[:, i] - Q_xn * Q_ddn * [delta_ddN - afixed];
        Q_xf = Q_xb - Q_xn * P_Q * Q_nx;  # 协方差
        Xb[:, i] = Xb0 + delta_xbf;
        ddN1[:, i] = ddN1[:, i] + afixed;
    return Xb, ddN1