import numpy as np
import matplotlib.pyplot as plt

L=1
N=101
Pext=(0.472 * 40**2)/(2*L**2)
deltax=1/(N-1)
deltat=0.49*((deltax)**2)
To=309.65*((0.56)/(Pext * L**2))
ta=0.025
rang_temps=round(ta/deltat)

T =np.zeros(N)
T_nuevo=np.zeros(N)
x = np.linspace(0,2,N)
T.fill(To)

for n in range(rang_temps):
    for i in range(1,N-1):
        T_nuevo[i]= T[i] + (deltat)/(deltax)**2 * (T[i+1] - 2*T[i] + T[i-1] + (deltax)**2)
        T_nuevo[0]=To
        T_nuevo[-1]=To
    
    T = T_nuevo.copy()
    if any(((T[0:36]/((0.56)/(Pext * L**2)))-273.15) > 50) == True:
        break

print(n)

plt.plot(x,(T/((0.56)/(Pext * L**2)))-273.15)
plt.grid(True)
plt.xlabel('Posició (cm)')
plt.ylabel('Graus Centígrads')
plt.hlines(50,0,2,color='r',linestyle='--')
plt.hlines(80,0,2,color='r',linestyle='--')
plt.vlines(0.75,36.5,80,color='r',linestyle='--')
plt.vlines(1.25,36.5,80,color='r',linestyle='--')
plt.show()