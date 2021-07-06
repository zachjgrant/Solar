#!/usr/bin/env python

#Serial version of solar energy calculations

import numpy as np
import math as m
import statistics as stat
#import matlibplot.pyplot as plt
import pandas as pd

data = pd.read_excel('solar_data.xlsx');

mon = np.linspace(1,12,12);
loc = np.linspace(1,239,239);

for i in mon:
	for j in loc:
		if mon == 1:
			day = np.linspace(1,31,31);
		elif mon == 2:
			day = np.linspace(32,59,28);
		elif mon == 3:
			day = np.linspace(60,90,31);
		elif mon == 4:
			day = np.linspace(91,120,30);
		elif mon == 5:
                        day = np.linspace(121,151,31);
		elif mon == 6:
			day = np.linspace(152,181,30);
		elif mon == 7:
			day = np.linspace(182,212,31);
		elif mon == 8:
			day = np.linspace(213,243,31);
		elif mon == 9:
			day = np.linspace(244,273,30);
		elif mon == 10:
			day = np.linspace(274,304,31);
		elif mon == 11:
			day = np.linspace(305,334,30);
		else:
			day = np.linspace(335,366,31);

		for n in day:
			B = (n-1) * 360/365;
			delta = (180/np.pi)*(0.006819 - 0.399912*np.cos(B) + 0.070257*np.sin(B) - 0.006758*np.cos(2*B) + 0.000907*np.sin(2*B) - 0.002697*np.cos(3*B) + 0.00148*np.sin(3*B));


printf(data);			
