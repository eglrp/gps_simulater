import numpy as np

from com.nuaa.gps.nrc.common_func.gnss_carrier_dd import gnss_carrier_dd
from com.nuaa.gps.nrc.common_func.gnss_observations import gnss_observations
from com.nuaa.gps.nrc.common_func.lamcacamb import LamCacAmb
from com.nuaa.gps.nrc.common_func.lsbaseline_errormodel import LSbaseline_ErrorModel
from com.nuaa.gps.nrc.common_func.satellite_positions import satellite_positions
from com.nuaa.gps.nrc.common_func.user_positions import user_positions
from com.nuaa.gps.nrc.entity.eph import eph
from com.nuaa.gps.nrc.config.configuration import global_gpsConfig as sign_set, global_transfer_Config as trans_set
from com.nuaa.gps.nrc.utils.file_reader_utils import readFile

data_file_directory = '../../../../data/rinex/'


class main_service:
    eph_list_brdc0060 = []
    eph_list_brdc1980 = []
    eph_list_ipil1910 = []

    def __init__(self, data_file_directory):
        self.dic = {'brdc0060.10n': self.eph_list_brdc0060,
                    'brdc1980.17n': self.eph_list_brdc0060,
                    'lpil1910.09n': self.eph_list_brdc0060}
        self.parse_directory(data_file_directory)
        self.orbit = self.eph_list_brdc0060[0:32]

    def parse_directory(self, data_file_directory):
        [self.aa(data_file_directory, key) for key in self.dic.keys()]

    def aa(self, data_file_directory, name):
        file_path = data_file_directory + name;
        obj = self.dic[name]
        content = self.load_data(file_path)
        self.parse_file(content, obj)
        print('load data from {} size {}'.format(file_path, len(obj)))

    def parse_file(self, content, save_to_list):
        list = []
        for i in range(len(content) // 8):
            sub_content = " ".join(content[i * 8:(i + 1) * 8])
            eph_obj = eph(sub_content)
            list.append(eph_obj)
        save_to_list.extend(list)

    def load_data(self, file_path):
        content = readFile(file_path)
        return content

    def get_guan(self):
        for t in np.arange(start=sign_set.sign_set_start_time, stop=sign_set.sign_set_end_time,
                           step=sign_set.sign_set_td_period):
            xs_tx, vs_tx, time_m, time_r = satellite_positions(trans_set, t, self.orbit);
            xr_m, vr_m = user_positions(t, trans_set.T, trans_set.atti_M,
                                        trans_set.atti_rate_M, trans_set.veloB_M,
                                        trans_set.acceB_M, trans_set.posi_M);
            xr_r, vr_r = user_positions(t, trans_set.T, trans_set.atti_R,
                                        trans_set.atti_rate_R, trans_set.veloB_R,
                                        trans_set.acceB_R, trans_set.posi_R);
            i = round(1 + 5 * t);
            trans_set.Xb_real[:, i] = xr_m - xr_r;
            gnss_obs = gnss_observations(t, xs_tx, vs_tx, xr_m, vr_m, xr_r, vr_r)
            gnss_dd = gnss_carrier_dd(t, gnss_obs)

        for i in range(1, sign_set.sign_set_data_length):
            for k in range(1, sign_set.sign_set_prnmat):
                trans_set.Bk3[k, :, i] = [gnss_dd.dde1x[k, i],
                                          gnss_carrier_dd.dde1y[k, i],
                                          gnss_carrier_dd.dde1z[k, i]];

        Xb, ddN1 = LSbaseline_ErrorModel(gnss_dd, trans_set, gnss_obs);

        ddN0 = ddN1[:, 1] * np.ones((1, np.mat(ddN1, 2).shape));  # 基线初始模糊度矩阵
        ddN1_dop = ddN1 - ddN0;
        afixed, t1, Ps1 = LamCacAmb(gnss_dd.ddph1, ddN1, gnss_obs.Bk3, ddN1_dop, sign_set.sign_set_g_bo1, Xb);

        # 利用Lambda法求解的初始模糊度求取所有模糊度
        ddNF1 = afixed[:, 1] * np.ones((1, np.mat(ddN1, 2).shape)) + ddN1_dop;  # 所有历元模糊度  ???
        # 模糊度固定后再求基线向量
        Xb, Q_Xb = LSbaseline_Ambfixed(gnss_dd.ddph1, ddNF1, gnss_obs.Bk3, Xb);

        Baselinelength = np.mat(np.zeros((sign_set.sign_set_data_length, 1)))
        Baseline_real = Baselinelength
        Baseline_error = Baselinelength
        for i in range(sign_set.sign_set_data_length):
            Baselinelength[i, 1] = np.sqrt(np.square(Xb(1, i)) + np.square(Xb(2, i)) + np.square(Xb(3, i)));  # 基线矢量的长度
            Baseline_real[i, 1] = np.square(
                np.sqrt(np.square(trans_set.Xb_real(1, i)) + np.square(trans_set.Xb_real(2,
                                                                                         i))
                        + np.square(trans_set.Xb_real(
                    3, i))));
            Baseline_error[i, 1] = Baseline_real[i, 1] - Baselinelength[i, 1];

            if __name__ == '__main__':
                service = main_service(data_file_directory)
