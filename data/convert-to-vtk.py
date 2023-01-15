import h5py
import numpy as np
from vtk.util import numpy_support
from pyevtk.hl import gridToVTK, unstructuredGridToVTK
from pyevtk.vtk import VtkGroup


hf_grid = h5py.File('data-h5/grid.h5', 'r')
x = np.array(hf_grid.get('x'))
y = np.array(hf_grid.get('y'))
x_u = np.array(hf_grid.get('xu'))
x_v = np.array(hf_grid.get('xv'))
y_u = np.array(hf_grid.get('yu'))
y_v = np.array(hf_grid.get('yv'))

nx = x.shape[0]
ny = y.shape[0] + 1

dset = h5py.File('data-h5/swec.h5', 'r')
str_t = sorted([float(x) for x in list(dset.keys())])
Nt = len(str_t)

x_mgrid, y_mgrid = np.mgrid[100:10000:100, 25:975:50]

x_grid = np.array([x_mgrid])
y_grid = np.array([y_mgrid])
z_grid = np.zeros_like(x_grid)

g = VtkGroup("./v-group")
for idx, n in enumerate(str_t):
    grid_array = np.array([np.array(dset[str(n)]['v'])[:-1,:].T])
    filename = "t{time:.2f}".format(time=n)
    gridToVTK("./v-group/"+filename, x_grid, y_grid, z_grid, pointData={"v": grid_array})
    g.addFile(filepath=filename+".vts", sim_time=n)    
    
g.save()
