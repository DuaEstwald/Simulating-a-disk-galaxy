# This exercise pretends to recreate the 'Orbit of Planets' pdf


import numpy as np
import matplotlib.pyplot as plt

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

# -------------------------------

N = 1000

t0 = 0.
tf = 2.#*tint#*yr # internal units corresponding to 18.6e8yr

dt = (tf-t0)/N

r0 = 8.#*kpc_cm #kpc
phi0 = 0.
vphi = 200.#*kms #km/s

Rs = 20.#*kpc_cm
rho0 = 5932371.#*msun/(kpc_cm**3)
#mcenter = 5.e10#*msun #Msun # EL FALLO ESTA AQUI, TIENES QUE CALCULAR LA MASA DENTRO DEL RADIO DEL SOL, ESTO ES SOLO UNA ESTIMACION, CALCULALA EN CONDICIONES


x = np.zeros(N)
x[0] = r0*np.cos(phi0)

y = np.zeros(N)
y[0] = r0*np.sin(phi0)

z = np.zeros(N)
z[0] = 0.



#vphi = np.sqrt(G*m[0]/a0) # au/year, initial velocity
vx = np.zeros(N)
vx[0] = -vphi*np.sin(phi0)
vy = np.zeros(N)
vy[0] = vphi*np.cos(phi0)
vz = np.zeros(N)
vz[0] = 0.





# =======================================================
# ========== DERIVATION METHODS =========================
# =======================================================

def Euler(c,v,fc,fv,h):
    return c + fc(v)*h, v + fv(c)*h
    
def RK2(c,v,fc,fv,h):
    cmid = c + 0.5*fc(v)*h
    vmid = v + 0.5*fv(c)*h
    return c + fc(vmid)*h, v + fv(cmid)*h

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
    return - (G*mcenter/(np.sqrt(r[0]**2+r[1]**2+r[2]**2)**3.))*r


#def fvel(r):
#    ind = np.arange(len(m))
#    dr = r-coord[ind !=j,i]
#    return -G*np.sum((m[m!=m[j]]/(np.sqrt(dr[:,0]**2.+dr[:,1]**2.+dr[:,2]**2.)**3.))*dr.T,axis=1)


coorde = np.array((x,y,z))
vele = np.array((vx,vy,vz))

coordrk4 = np.array((x,y,z))
velrk4 = np.array((vx,vy,vz))
for i in range(N-1): 
#    for j in range(len(m)): # Esto solo funciona para coord.T, vel.T
    Rmaxe = np.sqrt(np.sum(coorde[:,i]**2))
    Rmaxrk4 = np.sqrt(np.sum(coordrk4[:,i]**2))
    mcenter = 4.*np.pi*rho0*(Rs**3.)*(np.log((Rs+Rmaxe)/Rs)-(Rmaxe/(Rs+Rmaxe)))
    coorde[:,i+1],vele[:,i+1] = Euler(coorde[:,i],vele[:,i],fcoord,fvel,dt)
    mcenter = 4.*np.pi*rho0*(Rs**3.)*(np.log((Rs+Rmaxrk4)/Rs)-(Rmaxrk4/(Rs+Rmaxrk4)))
    coordrk4[:,i+1],velrk4[:,i+1] = RK4(coordrk4[:,i],velrk4[:,i],fcoord,fvel,dt)        
#    coord[0] = np.zeros((N,3)) # Es importante poner esto a la hora de resolver el sistema porque sino estas restando en fvel con nan y no se resuelven las orbitas


# ================= PLOTS ==================

plt.ion()
plt.close('all')

#plt.plot(coord[0,:,0],coord[0,:,1],'x')
#for n in range(len(m-1)):
plt.plot(coorde[0],coorde[1],'g',label='Euler')
plt.plot(coordrk4[0],coordrk4[1],'b',label='RK4')
plt.legend()
plt.grid()
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title(str(N)+' steps')
plt.savefig(str(N)+'steps.png')






