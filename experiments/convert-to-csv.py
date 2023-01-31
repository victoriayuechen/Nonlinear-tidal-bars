import h5py
import numpy as np
import pandas as pd
from vtk.util import numpy_support
from pyevtk.hl import gridToVTK, unstructuredGridToVTK
from pyevtk.vtk import VtkGroup


hf_grid = h5py.File('fine-mesh-h5/grid.h5', 'r')
x = np.array(hf_grid.get('x'))

y = np.array(hf_grid.get('y'))
x_u = np.array(hf_grid.get('xu'))
x_v = np.array(hf_grid.get('xv'))
y_u = np.array(hf_grid.get('yu'))
y_v = np.array(hf_grid.get('yv'))

nx = x.shape[0] - 2
ny = y.shape[0] + 1

dset = h5py.File('fine-mesh-h5/swec.h5', 'r')
str_t = sorted([float(x) for x in list(dset.keys())])
Nt = len(str_t)
saved_grid = np.zeros((Nt,  nx*(ny-1)))


for idx, n in enumerate(str_t):
    grid_array = np.array(dset[str(n)]['zeta'])[:-1,1:-1].T
    saved_grid[idx, :] = grid_array.flatten()

np.savetxt("swec-fine.csv", saved_grid[1:], delimiter=",")

