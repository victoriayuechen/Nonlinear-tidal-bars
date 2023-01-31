import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

apvm = pd.read_csv('mean-diff-apvm.csv', delimiter=',').to_numpy()
diff = pd.read_csv('mean-diff-diff.csv', delimiter=',').to_numpy()
against_eachother = pd.read_csv('mean-diff-apvm-test.csv' , delimiter=',').to_numpy()

plt.plot(np.arange(apvm.shape[0]), apvm, label="APVM order 0 vs. Hepkema")
plt.plot(np.arange(diff.shape[0]), diff, label="Diffusion vs. Hepkema")
plt.plot(np.arange(against_eachother.shape[0]), against_eachother, label="Diffusion vs. APVM")

plt.xlabel('Time $(\\times 10^2)$ (s)')
plt.ylabel('Absolute Mean Difference')
plt.legend()
plt.title('Mean Absolute Difference over Time')
plt.show()
plt.savefig("mean-difference.svg")
plt.close()