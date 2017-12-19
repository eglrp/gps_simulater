import math
import numpy as np


def geo_ccef(posi_geo):
    re = 6378137;  # GPS(WGS - 84) Ellipsoid semi - major axis[m]
    f = 1 / 298.257223563;  # GPS(WGS - 84) Ellipsoid flattening
    wie = 7.2921151467e-5;  # GPS Angular velocity of the Earth rotation[rad / s]
    g = 9.7803698;

    # user position(rad)
    long = posi_geo[0] * math.pi / 180.0;
    lati = posi_geo[1] * math.pi / 180.0;
    heig = posi_geo[2];

    # earth ellipticity
    Rn = re * (1 + f * math.sin(lati) * math.sin(lati));

    X = (Rn + heig) * math.cos(lati) * math.cos(long);
    Y = (Rn + heig) * math.cos(lati) * math.sin(long);
    Z = (Rn * np.square(1 - f) + heig) * math.sin(lati);
    return np.array([X, Y, Z])
