import re

import numpy as np
from pandas import Series


class eph:
    def __init__(self, text):
        if text is not None and len(text) > 0:
            regex = re.compile('\s+')
            s = regex.split(str.strip(text))
            self.temp_svprn = s[0]
            self.year = s[1]
            self.month = s[2]
            self.day = s[3]
            self.hour = s[4]
            self.minute = s[5]
            self.second = s[6]
            self.position = s[7:-1]
