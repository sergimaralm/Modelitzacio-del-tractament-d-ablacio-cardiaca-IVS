import matplotlib.pyplot as plt
import numpy as np

L = 2
N = 101
Pext=(0.472 * 40**2)/(2*L**2)
deltax=1/(N-1)
ta=0.025
x = np.linspace(0,2,N)

#Constants

deltat=0.5*deltax**2
To=309.65*((0.56)/(Pext * L**2))
rang_temps=round(ta/deltat)
gamma = deltat/deltax**2
T0 = To*np.ones(N)
T1 = np.zeros(N)

#Jacobi per la matriu Mx = b, on x serà T_(n+1)
M = np.identity(N) + (gamma)*(np.diag(2*np.ones(N), 0) + np.diag(-np.ones(N-1), -1) + np.diag(-np.ones(N-1), 1))

def jacobi(M, T0):
    b = deltat*np.ones(N) + T0
    for i in range(round(50)):        
        T1 = (1/(1+2*gamma))*((-np.diag(np.diag(M,-1),-1) - np.diag(np.diag(M,1), 1)) @ T0 + b)
        T1[0] = To
        T1[-1] = To
        T0 = T1.copy()
    return T1
#Metode implicit
for n in range(rang_temps):
    T1 = jacobi(M, T0)
    T0 = T1.copy()
    # if any(((T1[0:36]/((0.56)/(Pext * L**2)))-273.15) > 50) == True:
    #     temps = n*deltat
    #     print(temps)
    #     break

plt.plot(x,(T1/((0.56)/(Pext * L**2)))-273.15)
plt.grid(True)
plt.xlabel('Posició (cm)')
plt.ylabel('Graus Centígrads')
plt.hlines(50,0,2,color='r',linestyle='--')
plt.hlines(80,0,2,color='r',linestyle='--')
plt.vlines(0.75,36.5,55,color='r',linestyle='--')
plt.vlines(1.25,36.5,55,color='r',linestyle='--')

plt.show()


