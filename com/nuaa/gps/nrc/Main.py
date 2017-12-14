import numpy as np

from com.nuaa.gps.nrc.entity.eph import eph
from com.nuaa.gps.nrc.utils.file_reader_utils import readFile

data_file_directory = '../../../../data/rinex/'


class main_service:
    eph_list_brdc0060 = []
    eph_list_brdc1980 = []
    eph_list_ipil1910 = []

    def __init__(self, data_file_directory):
        self.parse_directory(data_file_directory)

    def parse_directory(self, data_file_directory):

        file_path = data_file_directory + 'brdc0060.10n';
        content = self.load_data(file_path)
        list = self.parse_file(content)
        self.eph_list_brdc0060 = list
        print('load data from {} size {}'.format(file_path, len(self.eph_list_brdc0060)))

        file_path = data_file_directory + 'brdc1980.17n';
        content = self.load_data(file_path)
        list = self.parse_file(content)
        self.eph_list_brdc1980 = list
        print('load data from {} size {}'.format(file_path, len(self.eph_list_brdc1980)))

        file_path = data_file_directory + 'lpil1910.09n';
        content = self.load_data(file_path)
        list = self.parse_file(content)
        self.eph_list_ipil1910 = list
        print('load data from {} size {}'.format(file_path, len(self.eph_list_ipil1910)))

    def parse_file(self, content):
        list = []
        for i in range(len(content) // 8):
            sub_content = " ".join(content[i * 8:(i + 1) * 8])
            eph_obj = eph(sub_content)
            list.append(eph_obj)
        return list

    def load_data(self, file_path):
        content = readFile(file_path)
        return content


if __name__ == '__main__':
    service = main_service(data_file_directory)
