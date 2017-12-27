import numpy as np


class transferConfig:
    def __init__(self, gpsConfig):
        self.lamada1 = gpsConfig.sign_set_g_bo1;
        self.posi_GEO = np.zeros((3, int(gpsConfig.sign_set_data_length)))  # 地理坐标系下位置
        self.T = gpsConfig.sign_set_td_period;
        self.atti_M = gpsConfig.sign_set_init_atti_m;
        self.atti_rate_M = np.zeros((3, 1));
        self.posi_M = gpsConfig.sign_set_init_posi_m;
        self.veloB_M = gpsConfig.sign_set_init_velo_m;
        self.acceB_M = np.zeros((3, 1));
        self.atti_R = gpsConfig.sign_set_init_atti_r;
        self.atti_rate_R = np.zeros((3, 1));
        self.posi_R = gpsConfig.sign_set_init_posi_r;
        self.veloB_R = gpsConfig.sign_set_init_velo_r;
        self.acceB_R = np.zeros((3, 1));

        self.ph1_M = np.zeros((len(gpsConfig.sign_set_prnmat), int(gpsConfig.sign_set_data_length)));
        self.ph2_M = self.ph1_M
        self.ph5_M = self.ph1_M

        self.ph1_R = self.ph1_M
        self.ph2_R = self.ph1_M
        self.ph5_M = self.ph1_M

        # 卫星到载体双差方向余弦矢量
        self.Bk3 = np.zeros((len(gpsConfig.sign_set_prnmat) - 1, 3, int(gpsConfig.sign_set_data_length)));

        self.ddN1 = np.zeros((len(gpsConfig.sign_set_prnmat) - 1, int(gpsConfig.sign_set_data_length)));
        self.ddN1_dop = self.ddN1
        self.ddN2 = self.ddN1
        self.ddN2_dop = self.ddN1
        self.ddN5 = self.ddN1
        self.ddN5_dop = self.ddN1
        self.Xb_real = np.zeros((3, int(gpsConfig.sign_set_data_length)));
