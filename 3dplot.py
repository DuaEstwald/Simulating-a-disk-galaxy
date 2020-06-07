from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

coord = np.load('cEXTRA100010.npy')

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')


for n in range(1000):
    ax.plot(coord[n,:,0],coord[n,:,1],coord[n,:,2],color='olivedrab',alpha=0.5,linewidth = 0.4)
#    cset = ax.contour(coord[n,:,0], coord[n,:,1], coord[n,:,2], zdir='z')
#cset = ax.contour(X, Y, Z, zdir='x', cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, zdir='y', cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_xlim(-40, 40)
ax.set_ylabel('Y')
ax.set_ylim(-40, 40)
ax.set_zlabel('Z')
ax.set_zlim(-20, 20)

plt.savefig('cEXTRA100010.png')
