#!/usr/bin/env python

import numpy as np
#import matlibplot.pyplot as plt
import pandas as pd
from mpi4py import MPI

comm = MPI.COMM_WORLD;
rank = comm.Get_rank();
size = comm.Get_size();

if rank == 0:
	a = [1,2,3];
else:
	a = [4,5,6];

b = np.zeros((2,3));
b = comm.gather(a);

if rank == 0:
	b = np.transpose(b);
	b = np.sum(b,axis=1);
	print(b);
