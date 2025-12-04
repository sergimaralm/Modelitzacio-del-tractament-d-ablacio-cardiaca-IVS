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

# Discretitzacions
deltax = 1 / (N - 1)
deltat = 0.5*((deltax)**2)   

T_list = []

gamma = deltat / deltax**2
# Condició inicial
T = To*np.ones(N)
T_list.append((T/((0.56)/(Pext * L**2)))-273.15)

# Matrius A i B
A = np.identity(N) + (gamma)*(np.diag(np.ones(N), 0) + np.diag(-0.5*np.ones(N-1), -1) + np.diag(-0.5*np.ones(N-1), 1))
B = np.identity(N) + (gamma)*(np.diag(-np.ones(N), 0) + np.diag(0.5*np.ones(N-1), -1) + np.diag(0.5*np.ones(N-1), 1))

# Bucle
t = 0.0
while t < ta:
    b = (B @ T) + deltat 
    T = jacobi(A, b, T)  
    T_list.append(((T/((0.56)/(Pext * L**2)))-273.15))
    # if any(((T[0:36]/((0.56)/(Pext * L**2)))-273.15) > 50) == True:
    #     print(t)
    #     break
    t += deltat

#for i in range(1,len(T_list)):
