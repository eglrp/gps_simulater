import numpy as np
import math

from com.nuaa.gps.nrc.config.transConfig import transferConfig
from com.nuaa.gps.nrc.entity.eph import eph
from com.nuaa.gps.nrc.config.configuration import gpsConfig
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

        self.gps_config = gpsConfig();
        self.transfer_Config = transferConfig(self.gps_config)

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
        for t in np.arange(start=self.gps_config.sign_set_start_time, stop=self.gps_config.sign_set_end_time,
                           step=self.gps_config.sign_set_td_period):
            xs_tx, vs_tx, time_m, time_r = self.satellite_positions(t, self.orbit);
            xr_m, vr_m = self.user_positions(t, 0);
            xr_r, vr_r = self.user_positions(t, 1);
            i = round(1 + 5 * t);
            self.transfer_Config.Xb_real[:, i] = xr_m - xr_r;
            gnss_observations = self.gnss_observations(t, xs_tx, vs_tx, xr_m, vr_m, xr_r, vr_r)
            self.gnss_carrier_dd = self.gnss_carrier_dd(t, gnss_observations)

        for i in range(1, self.gps_config.sign_set_data_length):
            for k in range(1, self.gps_config.sign_set_prnmat):
                self.transfer_Config.Bk3[k, :, i] = [self.transfer_Config.self.gnss_carrier_dd.dde1x[k, i],
                                                     self.gnss_carrier_dd.dde1y[k, i],
                                                     self.gnss_carrier_dd.dde1z[k, i]];

    def user_positions(self, t, flag):
        xr_m = ''
        vr_m = ''
        return xr_m, vr_m

    def gnss_observations(self, t, xs_tx, vs_tx, xr_m, vr_m, xr_r, vr_r):
        return GNSSObservations

    def gnss_carrier_dd(self, t):
        return GNSS_Carrier_DD

    def satellite_positions(self, t, orbit):
        nsat = len(self.gps_config.sign_set_prnmat);
        time_m = np.zeros(nsat, 1);
        xs_tx = np.zeros(nsat, 3);
        vs_tx = np.zeros(nsat, 3);
        # position in GEOREF : long (deg)  lati (deg)  alti (m)
        posi_m = self.gps_config.sign_set_init_posi_m
        posi_r = self.gps_config.sign_set_init_posi_r

        # Master
        trace_m_posi__geo = posi_m;
        trace_m_posi__ecef = self.geo2ecef(posi_m);

        # Rover
        traceR_posi_GEO = posi_r;
        traceR_posi_ECEF = self.geo2ecef(posi_r);

        Xr1 = traceR_posi_ECEF;
        t_orbit = self.gps_config.sign_set_t_sate + t;
        time_r = ''

        for j in range(nsat):
            # 解算卫星位置速度
            Xs, Vs = self.sate_posivelo(t_orbit, orbit[j, :]);  # 信号接收时刻卫星位置

        return xs_tx, vs_tx, time_m, time_r

    def sate_posivelo(self, t, orbit):
        GM = 3.986005e14;  # 地球引力常数
        wie = 7.2921151467e-5;  # 地球自转角速度
        satp = ''
        satv = ''

        return satp, satv

    def geo2ecef(self, posi_geo):
        Re = 6378137;  # GPS(WGS - 84) Ellipsoid semi - major axis[m]
        f = 1 / 298.257223563;  # GPS(WGS - 84) Ellipsoid flattening
        Wie = 7.2921151467e-5;  # GPS Angular velocity of the Earth rotation[rad / s]
        g = 9.7803698;

        # user position(rad)
        long = posi_geo[1] * math.pi / 180.0;
        lati = posi_geo[2] * math.pi / 180.0;
        heig = posi_geo[3];

        # earth ellipticity
        Rn = Re * (1 + f * math.sin(lati) * math.sin(lati));

        X = (Rn + heig) * math.cos(lati) * math.cos(long);
        Y = (Rn + heig) * math.cos(lati) * math.sin(long);
        Z = (Rn * np.square(1 - f) + heig) * math.sin(lati);
        return np.array([X, Y, X])


if __name__ == '__main__':
    service = main_service(data_file_directory)
