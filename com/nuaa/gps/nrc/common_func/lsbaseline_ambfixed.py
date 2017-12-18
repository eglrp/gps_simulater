from com.nuaa.gps.nrc.config.configuration import global_gpsConfig
import numpy as np


def LSbaseline_Ambfixed(ddph1, ddN1, Bk3, Xb):
    for i in range(global_gpsConfig.sign_set_data_length):
        # 最小二乘法求解基线矢量
        ks = len(global_gpsConfig.sign_set_prnmat);

        D1 = 2 * (np.ones(ks - 1, ks - 1) + np.eye(ks - 1)) * np.square(
            global_gpsConfig.sign_set_phase * global_gpsConfig.sign_set_g_bo1);  # 求权矩阵  GPS-王惠南 P161页  公式（6.104）
        P1 = (D1).I;  # 求协方差矩阵 最小二乘中的权重值

        Af = Bk3[:, :, i];  # 载波的矩阵方程基线改正数系数A      (A B)*[delta_Xb  delta_ddN]'=L
        Bf = global_gpsConfig.sign_set_g_bo1 * np.eye(ks - 1);  # 双差模糊度改正数系数B
        Lf = global_gpsConfig.sign_set_g_bo1 * (ddph1[:, i] + ddN1[:, i] + Af * Xb[:, i]);
        N11 = Af.T * P1 * Af;
        N12 = Af.T * P1 * Bf;
        N21 = Bf.T * P1 * Af;
        N22 = Bf.T * P1 * Bf;
        U1 = Af.T * P1 * Lf;
        U2 = Bf.T * P1 * Lf;
        N_0 = N21 * (N11).I * N12;
        P_Q = N22 - N_0;
        Segma_A = (P_Q).I;
        delta_ddN = Segma_A * (U2 - N21 * (N11).I * U1);  # 模糊度改正数
        Q_ddn = Segma_A;  # 模糊度协方差矩阵
        delta_xb = (N11).I * (U1 - N12 * delta_ddN);  # 基线改正数
        Q_xb = (N11).I + (N11).I * N12 * Segma_A * N21 * (N11).I;  # 基线协方差矩阵
        Q_xn = -(N11).I * N12 * Segma_A;  # 基线-模糊度协方差矩阵
        Q_nx = -Segma_A * N21 * (N11).I;
        Xb[:, i] = Xb[:, i] + delta_xb;  # 基线向量
    return Xb, Q_xb
