import h5py
import numpy as np
from vtk.util import numpy_support
from pyevtk.hl import gridToVTK, unstructuredGridToVTK
from pyevtk.vtk import VtkGroup


hf_grid = h5py.File('new-data/grid.h5', 'r')
x = np.array(hf_grid.get('x'))
y = np.array(hf_grid.get('y'))
x_u = np.array(hf_grid.get('xu'))
x_v = np.array(hf_grid.get('xv'))
y_u = np.array(hf_grid.get('yu'))
y_v = np.array(hf_grid.get('yv'))

nx = x.shape[0]
ny = y.shape[0] + 1

dset = h5py.File('hepkema-100/swec.h5', 'r')
str_t = sorted([float(x) for x in list(dset.keys())])
Nt = len(str_t)

x_mgrid, y_mgrid = np.mgrid[0:10100:100, 0:1050:50]

x_grid = np.array([x_mgrid])
y_grid = np.array([y_mgrid])
z_grid = np.zeros_like(x_grid)

# Change to v, u, or zeta
g = VtkGroup("./u-group")
for idx, n in enumerate(str_t):
    grid_array = np.array([np.array(dset[str(n)]['u'])[:,:].T])
    filename = "t{time:.2f}".format(time=n)
    gridToVTK("./u-group/"+filename, x_grid, y_grid, z_grid, pointData={"u": grid_array})
    g.addFile(filepath=filename+".vts", sim_time=n)    
    
g.save()

