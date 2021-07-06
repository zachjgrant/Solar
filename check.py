#!/usr/bin/env python

import numpy as np
#import matlibplot.pyplot as plt
import pandas as pd

data = np.array(pd.read_excel('solar_data.xlsx'));
lat = data[:,3];
loc = data[:,1];

a = np.cos(90*np.pi/180);

b = np.array([[1,2,3],[4,5,6]]);
c = np.array([[2,3,4],[5,6,7]]);

d = np.multiply(b,c);

e = np.concatenate((b,c),axis=0);

f = e[0:2,:];

print(f);
