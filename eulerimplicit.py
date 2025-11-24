import matplotlib.pyplot as plt
import numpy as np

L = 1
N = 101
Pext=(0.472 * 40**2)/(2*L**2)
deltax=1/(N-1)
ta=0.025
x = np.linspace(0,2,N)

#Constants

deltat=deltax**2
To=309.65*((0.56)/(Pext * L**2))
rang_temps=round(ta/deltat)
gamma = deltat/deltax**2
T0 = To*np.ones(N)
T1 = np.zeros(N)
T2 = np.zeros(N)

#Jacobi per la matriu Mx = b, on x serà T_(n+1)

def jacobi(T1, T0):
    M = np.identity(N) + (gamma/2)*(np.diag(2*np.ones(N), 0) + np.diag(-np.ones(N-1), -1) + np.diag(-np.ones(N-1), 1))
    b = (gamma*(deltax**2))*np.ones(N) + (np.identity(N)) @ T1  - (2/gamma) * (M - np.identity(N)) @ T0
    l = 0
    for i in range(round(25)):        
        T2 = (1/(1+gamma))*((np.diag(-np.diag(M,-1),-1) + np.diag(-np.diag(M,1), 1)) @ T0 + b)
        T2[0] = To
        T2[-1] = To
        l += 1
        T0 = T2.copy()
    return T2, T0, l



# Per al primer pas temporal fem el mètode explícit, doncs necessitem dos punts per començar, escollim deltat = 0.25*deltax**2
for i in range(1,N-1):
        T1[i]= T0[i] + (0.49) * (T0[i+1] - 2*T0[i] + T0[i-1] + (deltax)**2)
        T1[0]=To
        T1[-1]=To
#Metode implicit
l = []
for n in range(rang_temps-1):
    T2, T0, l_temp = jacobi(T1, T0)
    l.append(l_temp)
    T2[0] = To
    T2[-1] = To
    T1[0] = To
    T1[-1] = To
    T1 = T2.copy()

print(l, rang_temps)

plt.plot(x,(T2/((0.56)/(Pext * L**2)))-273.15)
plt.grid(True)
plt.xlabel('Posició (cm)')
plt.ylabel('Graus Centígrads')
plt.hlines(50,0,2,color='r',linestyle='--')
plt.hlines(80,0,2,color='r',linestyle='--')
plt.vlines(0.75,36.5,55,color='r',linestyle='--')
plt.vlines(1.25,36.5,55,color='r',linestyle='--')

plt.show()


