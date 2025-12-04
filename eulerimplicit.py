import matplotlib.pyplot as plt
import numpy as np
from funcions import analitica
from funcions import jacobi

#Paràmetres
L = 2
N = 101
Pext=(0.472 * 40**2)/(2*L**2)
To=309.65*((0.56)/(Pext * L**2))
ta=0.025

#Discretitzacions
x = np.linspace(0,2,N)
deltax=1/(N-1)
deltat_list = [(deltax)**2, 0.5*((deltax)**2)]  

T_list = []
err_list = []
for deltat in deltat_list:

    gamma = deltat/deltax**2
    T = To*np.ones(N) # Condició inicial

    M = np.identity(N) + (gamma)*(np.diag(2*np.ones(N), 0) + np.diag(-np.ones(N-1), -1) + np.diag(-np.ones(N-1), 1)) #matriu M
    #Bucle
    t = 0.0
    while t < ta:
        b = deltat*np.ones(N) + T #Vector b
        T = jacobi(M, b, T)
        
        t += deltat
    T_list.append((T/((0.56)/(Pext * L**2)))-273.15)
    err_list.append(abs(analitica(ta) - ((T/((0.56)/(Pext * L**2)))-273.15)))
    

# Gràfica per cada cas
plt.plot(x, T_list[0], label = f"$\\Delta T = (\\Delta X)²$")
plt.plot(x, T_list[1], label = f"$\\Delta T = 0.5 (\\Delta X)²$")
plt.hlines(50,0,2,color='r',linestyle='--')
plt.hlines(80,0,2,color='r',linestyle='--')
plt.vlines(0.75,36.5,80,color='r',linestyle='--')
plt.vlines(1.25,36.5,80,color='r',linestyle='--')
plt.xlabel('Posició (cm)')
plt.ylabel('Temperatura  (ºC)')
plt.gca().tick_params(direction="in")
plt.legend(loc=1, frameon=False, borderaxespad= 1)
plt.savefig('figures/implicit.png', bbox_inches='tight')
plt.show()

# Gràfica d'errors absoluts
plt.plot(x, err_list[0], label = f"$\\Delta T = (\\Delta X)²$")
plt.plot(x, err_list[1], label = f"$\\Delta T = 0.5 (\\Delta X)²$")
plt.xlabel('Posició (cm)')
plt.ylabel('Error absolut (ºC)')
plt.gca().tick_params(direction="in")
plt.legend(loc=1, frameon=False, borderaxespad= 0)
plt.savefig('figures/err_implicit.png', bbox_inches='tight')
plt.show()

# Gràfica d'errros relatius percentuals
plt.plot(x, (err_list[0] / analitica(ta))*100, label = f"$\\Delta T = (\\Delta X)²$")
plt.plot(x, (err_list[1] / analitica(ta))*100, label = f"$\\Delta T = 0.5 (\\Delta X)²$")
plt.xlabel('Posició (cm)')
plt.ylabel('Error relatiu (%)')
plt.gca().tick_params(direction="in")
plt.legend(loc=1, frameon=False, borderaxespad= 0)
plt.savefig('figures/relerr_implicit.png', bbox_inches='tight')
plt.show()