import matplotlib.pyplot as plt
import numpy as np
import math
from codigoanalitica import analitica

#Jacobi per la matriu Ax = b, on x serà T_(n+1)
def jacobi(A, b, x0, max_iter=50):

    # Descomposició A = D + E + F
    D = np.diag(np.diag(A))
    E = np.tril(A, k=-1)     
    F = np.triu(A, k=1)      
   
    # Matriu D^{-1} i P
    D_inv = np.diag(1 / np.diag(D))
    P = -D_inv @ (E + F)

    x = x0.copy()
    for _ in range(max_iter):
        x_new = P @ x + D_inv @ b
        x_new[0] = x0[0]
        x_new[-1] = x0[-1]
        x = x_new
    return x

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
        # if any(((T1[0:36]/((0.56)/(Pext * L**2)))-273.15) > 50) == True:
        #     print(t)
        #     break
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
plt.ylabel('Error relatiu')
plt.gca().tick_params(direction="in")
plt.legend(loc=1, frameon=False, borderaxespad= 0)
plt.savefig('figures/relerr_implicit.png', bbox_inches='tight')
plt.show()


