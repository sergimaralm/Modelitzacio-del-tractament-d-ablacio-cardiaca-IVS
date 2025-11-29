import numpy as np
import matplotlib.pyplot as plt
from codigoanalitica import analitica

# Metóde Jacobi resoldre Ax = b
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

# Parametres
L = 2
N = 101
ta = 0.025    
Pext = (0.472*40**2)/(2*L**2)
To = 309.65*((0.56)/(Pext*L**2))  
x = np.linspace(0, 2, N)

deltax = 1 / (N - 1)
deltat_list = [(deltax)**2, 0.5*((deltax)**2)]    

T_list = []
err_list = []
for deltat in deltat_list:
    gamma = deltat / deltax**2
    # Condició inicial
    T = To*np.ones(N)

    # Matrius A i B
    A = np.zeros((N, N))
    B = np.zeros((N, N))

    for i in range(N):
        # Diagonal principal
        A[i, i] = 1 + gamma
        B[i, i] = 1 - gamma

        # Sub y super diagonals
        if i > 0:
            A[i, i-1] = -gamma / 2
            B[i, i-1] = gamma / 2
        if i < N-1:
            A[i, i+1] = -gamma / 2
            B[i, i+1] = gamma / 2

    # Bucle
    t = 0.0
    while t < ta:
        b = (B @ T) + deltat 
        T = jacobi(A, b, T)  
        # if any(((phi[0:36]/((0.56)/(Pext * L**2)))-273.15) > 50) == True:
        #     print(t)
        #     break
        t += deltat

    T_list.append((T/((0.56)/(Pext * L**2)))-273.15)
    err_list.append(abs(analitica(ta) - ((T/((0.56)/(Pext * L**2)))-273.15)))

# Gràfica per cada cas
plt.plot(x, T_list[0], label = f"$\\Delta T = \\Delta X)²$")
plt.plot(x, T_list[1], label = f"$\\Delta T = 0.5 \\Delta X)²$")
plt.hlines(50,0,2,color='r',linestyle='--')
plt.hlines(80,0,2,color='r',linestyle='--')
plt.vlines(0.75,36.5,80,color='r',linestyle='--')
plt.vlines(1.25,36.5,80,color='r',linestyle='--')
plt.xlabel('Posició (cm)')
plt.ylabel('Temperatura  (ºC)')
plt.legend()
plt.grid()
plt.show()

# Gràfica d'errors 
plt.plot(x, err_list[0], label = f"$\\Delta T = \\Delta X)²$")
plt.plot(x, err_list[1], label = f"$\\Delta T = 0.5 \\Delta X)²$")
plt.xlabel('Posició (cm)')
plt.ylabel('Error absolut (ºC)')
plt.legend()
plt.grid()
plt.show()