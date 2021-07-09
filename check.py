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
b = comm.gather(a,root=0);

if rank == 0:
	c = [[1,4],[2,5],[3,6]];
	b = np.transpose(b);
	d = np.multiply(b,c);
	print(b);
	print(c);
	print(d);
