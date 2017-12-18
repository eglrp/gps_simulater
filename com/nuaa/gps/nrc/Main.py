import numpy as np

from com.nuaa.gps.nrc.common_func.gnss_carrier_dd import gnss_carrier_dd
from com.nuaa.gps.nrc.common_func.gnss_observations import gnss_observations
from com.nuaa.gps.nrc.common_func.satellite_positions import satellite_positions
from com.nuaa.gps.nrc.common_func.user_positions import user_positions
from com.nuaa.gps.nrc.entity.eph import eph
from com.nuaa.gps.nrc.config.configuration import global_gpsConfig, global_transfer_Config
from com.nuaa.gps.nrc.entity.gnss_carrier_dd import GNSS_Carrier_DD
from com.nuaa.gps.nrc.entity.gnss_observations import GNSSObservations
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
        for t in np.arange(start=global_gpsConfig.sign_set_start_time, stop=global_gpsConfig.sign_set_end_time,
                           step=global_gpsConfig.sign_set_td_period):
            xs_tx, vs_tx, time_m, time_r = satellite_positions(global_transfer_Config, t, self.orbit);
            xr_m, vr_m = user_positions(t, global_transfer_Config.T, global_transfer_Config.atti_M,
                                             global_transfer_Config.atti_rate_M, global_transfer_Config.veloB_M,
                                             global_transfer_Config.acceB_M, global_transfer_Config.posi_M);
            xr_r, vr_r = user_positions(t, global_transfer_Config.T, global_transfer_Config.atti_R,
                                        global_transfer_Config.atti_rate_R, global_transfer_Config.veloB_R,
                                        global_transfer_Config.acceB_R, global_transfer_Config.posi_R);
            i = round(1 + 5 * t);
            global_transfer_Config.Xb_real[:, i] = xr_m - xr_r;
            var_gnss_observations = gnss_observations(t, xs_tx, vs_tx, xr_m, vr_m, xr_r, vr_r)
            var_gnss_carrier_dd = gnss_carrier_dd(t, var_gnss_observations)

        for i in range(1, global_gpsConfig.sign_set_data_length):
            for k in range(1, global_gpsConfig.sign_set_prnmat):
                global_transfer_Config.Bk3[k, :, i] = [global_transfer_Config.self.gnss_carrier_dd.dde1x[k, i],
                                                       gnss_carrier_dd.dde1y[k, i],
                                                       gnss_carrier_dd.dde1z[k, i]];



if __name__ == '__main__':
    service = main_service(data_file_directory)
