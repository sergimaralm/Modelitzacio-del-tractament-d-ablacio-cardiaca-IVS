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
Ta = 0.025    
Pext=(0.472 * 40**2)/(2*L**2)
To=309.65*((0.56)/(Pext * L**2))  

dx = 1 / (N - 1)
dt = 0.5* dx**2         

r = dt / dx**2
x = np.linspace(0, 2, N)

# Condició inicial
phi = To*np.ones(N)

# Matrius A i B
A = np.zeros((N, N))
B = np.zeros((N, N))

for i in range(N):
    # Diagonal principal
    A[i, i] = 1 + r
    B[i, i] = 1 - r

    # Sub y super diagonals
    if i > 0:
        A[i, i-1] = -r / 2
        B[i, i-1] = r / 2
    if i < N-1:
        A[i, i+1] = -r / 2
        B[i, i+1] = r / 2

# Bucle
t = 0.0
while t < Ta:
    b = (B @ phi) + dt 
    phi = jacobi(A, b, phi)  
    # if any(((phi[0:36]/((0.56)/(Pext * L**2)))-273.15) > 50) == True:
    #     print(t)
    #     break
    t += dt

plt.plot(x,(phi/((0.56)/(Pext * L**2)))-273.15)
plt.grid(True)
plt.xlabel('Posició (cm)')
plt.ylabel('Graus Centígrads')
plt.hlines(50,0,2,color='r',linestyle='--')
plt.hlines(80,0,2,color='r',linestyle='--')
plt.vlines(0.75,36.5,55,color='r',linestyle='--')
plt.vlines(1.25,36.5,55,color='r',linestyle='--')

plt.show()