import matplotlib.pyplot as plt
import numpy as np
import matplotlib as m
from scipy import optimize
m.rcParams.update({'font.size':10})

f6, f8 =np.genfromtxt('palomar2.txt', dtype=float, comments='#', usecols=(2,3), unpack=True)
c = f6-f8

fig, axs = plt.subplots(1, 2, sharey=True)

axs[0].scatter(c, f6, color="red", s=0.2)
axs[1].scatter(c, f8, color="blue", s=0.2)

fig.tight_layout()
ax = plt.gca()
ax.invert_yaxis()
plt.show()