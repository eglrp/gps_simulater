import numpy as np

from com.nuaa.gps.nrc.config.transConfig import transferConfig


class gpsConfig:
    def __init__(self):
        self.sign_set_c = 299792458
        self.sign_set_g = 9.7803698
        self.sign_set_ell_a_gps = 6378137
        self.sign_set_ell_a_bds = 6378136
        self.sign_set_ell_f_gps = 1 / 298.257223563
        self.sign_set_ell_f_bds = 1 / 298.257222101
        self.sign_set_ell_e_gps = np.sqrt(1 - np.square(1 - self.sign_set_ell_f_gps))
        self.sign_set_ell_e_bds = np.sqrt(1 - np.square(1 - self.sign_set_ell_f_bds))
        self.sign_set_gm_gps = 3.986005e14
        self.sign_set_gm_bds = 3.986004418e14
        self.sign_set_earth_wie_gps = 7.2921151467e-5
        self.sign_set_earth_wie_bds = 7.292115e-5
        self.sign_set_pi_orbit = 3.1415926535898
        self.sign_set_l1 = 1575.42e6
        self.sign_set_l2 = 1227.60e6
        self.sign_set_l5 = 1176.45e6
        self.sign_set_g_bo1 = 0.1903
        self.sign_set_g_bo2 = 0.2442
        self.sign_set_g_bo5 = 0.2548
        self.sign_set_b1 = 1561.098e6
        self.sign_set_b2 = 1207.140e6
        self.sign_set_b3 = 1268.520e6
        self.sign_set_b_bo1 = 0.1920
        self.sign_set_b_bo2 = 0.2483
        self.sign_set_b_bo5 = 0.2363
        self.sign_set_code = 29.3
        self.sign_set_phase = 0.08
        self.sign_set_std_code = 3
        self.sign_set_std_phase = 0.03
        self.sign_set_std_phase_if = 0.009
        self.sign_set_delta_clock = 4.47e-09
        self.sign_set_delta_r_clock = 31
        self.sign_set_eph_file = 'G:\matlab\GNSS_INS组合相对导航\RINEX\brdc1980.17n'
        self.sign_set_prnmat = np.array([2, 5, 7, 9, 10, 25, 30])
        self.sign_set_start_time = 0
        self.sign_set_end_time = 20
        self.sign_set_td_period = 0.2
        self.sign_set_data_length = 1 + (self.sign_set_end_time - self.sign_set_start_time) / self.sign_set_td_period;
        self.sign_set_carr_init_phase = 0
        self.sign_set_wie = 7.2921151467e-5
        self.sign_set_init_atti_m = [0, 0, 0]
        self.sign_set_init_atti_r = [0, 0, 0]
        self.sign_set_init_posi_m = [118, 32, 500]
        self.sign_set_init_posi_r = [118, 32, 800]
        self.sign_set_init_velo_m = [0, 0, 0]
        self.sign_set_init_velo_r = [0, 0, 0]
        self.sign_set_t_sate = 86400


global_gpsConfig = gpsConfig()

global_transfer_Config = transferConfig(global_gpsConfig)
