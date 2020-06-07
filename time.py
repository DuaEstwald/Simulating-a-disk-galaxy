import numpy as np
import matplotlib.pyplot as plt

Npart, time = np.loadtxt('timesimulation.txt',unpack = True)

m1,n1 = np.polyfit(np.log10(Npart[:3]),np.log10(time[:3]),1)
m2,n2 = np.polyfit(np.log10(Npart[3:6]),np.log10(time[3:6]),1)
m3,n3 = np.polyfit(np.log10(Npart[6:9]),np.log10(time[6:9]),1)

m4,n4 = np.polyfit(np.log10(Npart[9:11]),np.log10(time[9:11]),1)


plt.loglog(Npart[:3],time[:3],'x',color='sienna')
plt.loglog(Npart[:3],10**(m1*np.log10(Npart[:3])+n1),label='DM only, m = %.2f'%m1,color='sienna')


plt.loglog(Npart[3:6],time[3:6],'x',color='firebrick')
plt.loglog(Npart[3:6],10**(m2*np.log10(Npart[3:6])+n2),label = 'All interactions, m = %.2f'%m2,color='firebrick')

plt.loglog(Npart[6:9],time[6:9],'x',color='olivedrab')
plt.loglog(Npart[6:9],10**(m3*np.log10(Npart[6:9])+n3),label='10 nearest neighbours, m = %.2f'%m3,color='olivedrab')

plt.loglog(Npart[9:11],time[9:11],'x',color='darkslategray')
plt.loglog(Npart[9:11],10**(m4*np.log10(Npart[9:11])+n4),color='darkslategray',label='50 nearest neighbours, m = %.2f'%m4)


plt.legend()
plt.grid()


plt.xlabel('Number of particles')
plt.ylabel('Code speed')

plt.savefig('timerunning.png')
