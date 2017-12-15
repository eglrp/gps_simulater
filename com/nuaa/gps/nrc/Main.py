import numpy as np

from com.nuaa.gps.nrc.entity.eph import eph
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


if __name__ == '__main__':
    service = main_service(data_file_directory)
