import h5py
import numpy as np
import pandas as pd
from vtk.util import numpy_support
from pyevtk.hl import gridToVTK, unstructuredGridToVTK
from pyevtk.vtk import VtkGroup
import matplotlib.pyplot as plt
from matplotlib import animation
EPS = 1e-10

## Read in Hepkema's grid
hf_grid = h5py.File('hepkema-100/grid.h5', 'r')
x = np.array(hf_grid.get('x'))

y = np.array(hf_grid.get('y'))
x_u = np.array(hf_grid.get('xu'))
x_v = np.array(hf_grid.get('xv'))
y_u = np.array(hf_grid.get('yu'))
y_v = np.array(hf_grid.get('yv'))

nx = x.shape[0] - 2
ny = y.shape[0] + 1

dset = h5py.File('hepkema-100/swec.h5', 'r')
str_t = sorted([float(x) for x in list(dset.keys())])
Nt = len(str_t)
zeta_hepkema = np.zeros((Nt, nx*(ny - 1)))


for idx, n in enumerate(str_t):
    # 20, 101 (y, x)
    # 101, 20
    grid_array = np.array(dset[str(n)]['zeta'])[:-1,1:-1]
    zeta_hepkema[idx, :] = grid_array.flatten()

# Obtain grids from Hepkema and Diffusion
zeta_diff = pd.read_csv('csv-files/diff-100-post.csv', delimiter=',').to_numpy()[1:]
zeta_hepkema = zeta_hepkema[1:-1]
zeta_hepkema[np.abs(zeta_hepkema) < EPS] = 0
diff_hepkema = zeta_hepkema - zeta_diff


# ## Save average difference over time  
diff_average_time = np.mean(np.abs(diff_hepkema), axis=1)
np.savetxt("mean-diff-diff.csv", diff_average_time, delimiter=",")

apvm =  pd.read_csv('csv-files/zeta-apvm-dt5.csv', delimiter=',').to_numpy()[1:]
apvm_diff = zeta_diff - apvm
test = np.mean(np.abs(apvm_diff), axis=1)
np.savetxt("mean-diff-apvm-test.csv", test, delimiter=",")


# Animation 
initial_frame = 100 *  diff_hepkema[0].reshape((ny - 1, nx))  / zeta_hepkema[0].reshape((ny - 1, nx))
initial_frame[np.isnan(initial_frame) < EPS] = 0
plt.ion()
fig = plt.figure()
plt.clf()
img = plt.imshow(initial_frame, extent=[0, x[-1], y[-1], 0], interpolation='none')
plt.gca().invert_yaxis()
plt.colorbar(img, orientation='horizontal',label='Percentage Difference')
tlt = plt.title(r"Difference in $\zeta(x,y,t)$ - Diffusion vs. Hepkema")


def animate(frame):
    image_array = None 
    if frame == 1:
        image_array = np.copy(initial_frame)
    else: 
        image_array = 100 * (diff_hepkema[frame].reshape((ny - 1, nx))) / zeta_hepkema[frame].reshape((ny - 1, nx))
    
    # Remove any NANs resulting from div by zero 
    image_array[np.isnan(image_array)] = 0
    img.set_array(image_array)
    tlt.set_text(r"Diffusion vs. Hepkema: $\zeta(x,y,t), t =$ " + str(np.round(frame*100,3)))

    return img

anim = animation.FuncAnimation(fig,func=animate, frames=range(1,diff_hepkema.shape[0]),interval=20,repeat=False)
plt.pause(Nt*5e-2)

# Save the animation and snapshot
writergif = animation.PillowWriter(fps=30) 
anim.save("zeta-diff-hepkema.gif", writer=writergif)
plt.savefig("zeta-diff-hepkema.svg")

plt.close()

