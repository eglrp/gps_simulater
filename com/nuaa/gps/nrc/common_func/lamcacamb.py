from com.nuaa.gps.nrc.config.configuration import global_gpsConfig
import numpy as np


def LamCacAmb(ddph1,ddN1,Bk3,ddN1_dop,bo1,Xb):
    i = 5;  # 观测历元数   2历元需同时观测可见卫星数至少4颗
    flag = 1;  # 标志位 表明模糊度多历元解算过程是否结束
    ks = len[global_gpsConfig.sign_set_prnmat]  # ks为卫星个数

    Af = Bk3;  # 载波的矩阵方程基线改正数系数A
    # while[flag]
    #  构造载波观测方程 L=A*x+B*y  利用LAMBDA法求解模糊度
    L = []
    A = []
    B = []
    N1_dop = [];
    # De = 2*[ones[ks-1,ks-1] + eye[ks-1]];
    De = 4 * np.square[global_gpsConfig.sign_set_phase * bo1] * np.eye[ks - 1];  # 求权矩阵   ？？
    Qe = [De].I  # 求协方差矩阵
    P = [];
    # k 根据需要选取多个时刻的数据进行模糊度求解
    for k in range[1, i]:
        L = np.mat([L, bo1 * [ddph1[:, 2 * k] - Af[:, :, 2 * k] * Xb[:, 2 * k]]]);  # 存储k个历元载波观测值   5*k  ???
        B = np.mat([B, np.eye[ks - 1]]);  # 存储k个历元模糊度系数矩阵
        N1_dop = np.mat([N1_dop, ddN1_dop[:, 2 * k]]);  # 存储k个历元模糊度值    ???

        P1 = np.mat([P, np.zeros((ks - 1) * (k - 1), ks - 1)]);  # 对应基线协方差矩阵
        P2 = np.mat([np.zeros([ks - 1, [ks - 1] * [k - 1]]), Qe]);  # 对应模糊度协方差矩阵
        P = np.mat([P1, P2]);  # 存储k个历元观测量协方差矩阵

        B = B * bo1;
        L = L - N1_dop * bo1;
        af1 = [B.T * P * B].I * B.T * P * L;
        Q1 = [B.T * P * B].I;  # Q1为模糊度的方差阵
        # ncands = 4;                             # 模糊度候选值个数
        ## LAMBDA法求解模糊度
        # [afixed,sqnorms] = lambda2[af1,Qx1,ncands];
        [afixed, sqnorm, Ps, Qzhat, Z, nfixed, mu] = LAMBDA[af1, Q1, 4];
    return afixed,i,Ps
