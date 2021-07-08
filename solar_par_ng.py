#!/usr/bin/env python

#Parallel version of solar energy calculations

#Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import geopandas as gpd
import seaborn as sb
import time
from mpi4py import MPI

#Get MPI characteristics
comm = MPI.COMM_WORLD;
rank = comm.Get_rank();
size = comm.Get_size();

#Start timer
start_time = time.time();

#Reading the excel data using pandasi
data = np.array(pd.read_excel('solar_data.xlsx'));

#Initializing the step operations by month and location
#mon = np.linspace(1,12,12,dtype=int);
loc = np.linspace(0,236,237,dtype=int);
st = np.linspace(0,49,50,dtype=int);

#Parsing the latitude and state data
lat = data[:,3];
state_raw = data[:,1];

#Parsing the clearness index values
KT = data[:,4:];

#Useful constant
pi = np.pi;
dtr = pi/180;

#Matrix for land area of each state in sq. meters
area = np.array([1477953,131171,134771,294207,403466,268431,12542,5047,138887,148959,16635,144669,214945,143793,92789,211754,102269,111898,20201,25142,79883,146435,206232,178040,121531,376962,125920,178711,198974,23187,19047,314161,284332,122057,105829,177660,248608,115883,2678,77857,196350,106798,676587,212818,102279,23871,172119,140268,62259,251470])*1000000;

#Initializing the storage variables
H_o_bar_loc = np.zeros((237,1));

#Beginning calculations by month by location
i = rank;
for j in loc:
    phi = lat[j];
    if i == 1:
        day = np.linspace(1,31,31);
        l = 31;
    elif i == 2:
        day = np.linspace(32,59,28);
        l = 28;
    elif i == 3:
        day = np.linspace(60,90,31);
        l = 31;
    elif i == 4:
        day = np.linspace(91,120,30);
        l = 30;
    elif i == 5:
        day = np.linspace(121,151,31);
        l = 31;
    elif i == 6:
        day = np.linspace(152,181,30);
        l = 30;
    elif i == 7:
        day = np.linspace(182,212,31);
        l = 31;
    elif i == 8:
        day = np.linspace(213,243,31);
        l = 31;
    elif i == 9:
        day = np.linspace(244,273,30);
        l = 30;
    elif i == 10:
        day = np.linspace(274,304,31);
        l = 31;
    elif i == 11:
        day = np.linspace(305,334,30);
        l = 30;
    else:
        day = np.linspace(335,366,31);
        l = 31;

        #Calculating the radiation on a horizontal surface for each day of the month
        H_o = np.linspace(0,l-1,l);
        a = 0; #step variable for storage
        for n in day:
                delta = 23.45 * np.sin(2*pi*(284+n)/365); #declination
                if -np.tan(phi*dtr)*np.tan(delta*dtr) > 1 or -np.tan(phi*dtr)*np.tan(delta*dtr) < -1:
                        H_o[a] = 0;
                else:
                        omega_s = np.degrees(np.arccos(-np.tan(phi*dtr)*np.tan(delta*dtr))); #sunset hour angle
                        H_o[a] = (24*3600*1367/pi) * (1+0.033*np.cos(2*pi*n/365)) * ((np.cos(phi*dtr)*np.cos(delta*dtr)*np.sin(omega_s*dtr))+(omega_s*dtr*np.sin(phi*dtr)*np.sin(delta*dtr))); #radiation on a horizontal surface outside the Earth's atmosphere
                a = a + 1;

        #Take the daily average for the month and save
        H_o_bar_loc[j] = np.mean(H_o);

if rank == 0:
	H_o_bar_loc_comb = np.zeros((12,237));
else:
	H_o_bar_loc_comb = None;

H_o_bar_loc_comb = np.array(comm.gather(H_o_bar_loc,root=0));

if rank == 0:
	H_o_bar_loc_comb = np.transpose(H_o_bar_loc_comb);
	#Calculate surface radiation by location
	H_bar_loc = np.zeros((237,12));
	H_bar_loc = np.multiply(KT,H_o_bar_loc_comb);
	H_bar_state = np.zeros((50,12));
	#Combine the values by state
	mon = np.linspace(1,12,12,dtype=int);
	for k in mon:
		p = -1;
		for m in st:
			n = 0;
			p = p+1;
			while state_raw[p] == state_raw[p+1] and p < 235:
				p = p+1;
				n = n+1;
			H_bar_state[m,k-1] = np.mean(H_bar_loc[p-n:p+1,k-1])/1000000;

	#Sorting by state
	state = [];
	for q in state_raw:
		if q not in state:
            		state.append(q);
	state = np.transpose(state);

	#Combine for all months and create data frame
	H_bar_state_year = np.mean(H_bar_state, axis=1);
	d_year = {'State': state, 'Value': H_bar_state_year};
	df_year = pd.DataFrame(data=d_year);
	print(df_year);

	#End timer
	end_time = time.time();
	print("---%s seconds---" %(end_time-start_time));
