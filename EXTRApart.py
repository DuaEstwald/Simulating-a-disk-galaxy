# This exercise pretends to recreate the 'Orbit of Planets' pdf


import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
import time
import sys

# =================================================================
# ============= INITIAL CONDITION =================================
# =================================================================


# ------------- UNITS ------------
G = 4.302e-6 # Gravitatonal constant in units kpc(km/s)2Msun-1
#G = 6.6708e-8 # Gravitatonal constant in CGS units [cm^3/g/s^2,]


msun = 1.98892e33 # Msol[grams]
kpc_cm = 3.085678e21 # kpc[cm]
kpc_km = 3.085678e16 #kpc[km]
yr = 31557600. # s
tint = 9.8e8 #yr
kms = 1.e5 #cm/s

# Its important to know that all of this units determines the internal unit of time:
# one internal time unit corresponding to 9.8e8 yr.

Npart = sys.argv[1]
neigh = sys.argv[3]
m,x0,y0,z0,vx0,vy0,vz0 = np.loadtxt('disk'+str(Npart)+'.txt',unpack=True)


# -------------------------------

Rs = 20.#*kpc_cm
rho0 = 5932371.#*msun/(kpc_cm**3)
eps = float(sys.argv[2])

N = 1000

t0 = 0.
tf = 2. # years

dt = (tf-t0)/N




x = np.zeros((N,len(m)))
x[0] = x0

y = np.zeros((N,len(m)))
y[0] = y0

z = np.zeros((N,len(m)))
z[0] = z0



#vphi = np.sqrt(G*m[0]/a0) # au/year, initial velocity
vx = np.zeros((N,len(m)))
vx[0] = vx0
vy = np.zeros((N,len(m)))
vy[0] = vy0
vz = np.zeros((N,len(m)))
vz[0] = vz0



# =======================================================
# ========== DERIVATION METHODS =========================
# =======================================================


def RK4(c,v,fc,fv,h):
    k1_c = fc(v)
    k1_v = fv(c)
    k2_c = fc(v + k1_v*h*0.5) # Tienes que ponerlo todo con las mismas dimensiones
    k2_v = fv(c + k1_c*h*0.5)
    k3_c = fc(v + k2_v*h*0.5)
    k3_v = fv(c + k2_c*h*0.5)
    k4_c = fc(v + k3_v*h)
    k4_v = fv(c + k3_c*h)
    return c + h/6. * (k1_c+2.*k2_c+2.*k3_c+k4_c), v + h/6. * (k1_v+2.*k2_v+2.*k3_v+k4_v)
# Euler
def fcoord(vel):
    return vel

def fvel(r):
    dr = r-coord[index[0][1:],i]
    return - (G*mcenter[j]/(np.sqrt(np.sum(r**2))**3.))*r-G*np.sum((m[index[0][1:]]/(np.sqrt(np.sum(dr**2,axis=1)+eps**2)**3.))*dr.T,axis=1)


coord = np.array((x,y,z)).T
vel = np.array((vx,vy,vz)).T

start = time.time()

for i in range(N-1): 
    Rmax = np.sqrt(np.sum(coord[:,i]**2,axis=1))
    mcenter = 4.*np.pi*rho0*(Rs**3.)*(np.log((Rs+Rmax)/Rs)-(Rmax/(Rs+Rmax)))
    tree = KDTree(coord[:,i],leaf_size=2)
    for j in range(m.shape[0]): # Esto solo funciona para coord.T, vel.T
        dist, index = tree.query(coord[j,i].reshape(1,-1),k=int(neigh))
        coord[j,i+1],vel[j,i+1] = RK4(coord[j,i],vel[j,i],fcoord,fvel,dt)        
#    coord[0] = np.zeros((N,3)) # Es importante poner esto a la hora de resolver el sistema porque sino estas restando en fvel con nan y no se resuelven las orbitas

t_sim = time.time() - start
with open('timesimulation.txt','a') as tfile:
    tfile.write('\t'+str(Npart)+'\t'+str(t_sim)+'\n')

# ================= PLOTS ==================
np.save('cEXTRA'+str(Npart)+str(neigh),coord)
np.save('vEXTRA'+str(Npart),vel)



plt.ion()
plt.close('all')

for n in range(m.shape[0]):
    plt.plot(coord[n,:,0],coord[n,:,1])
plt.grid()
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title(str(Npart)+' particles')
plt.savefig(str(Npart)+'_EXTRA.png')






