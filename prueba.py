import numpy as np

def jacobi(A, b, x0, tol=1e-10, max_iter=5000):

    # Descomposición A = D + E + F
    D = np.diag(np.diag(A))
    E = np.tril(A, k=-1)     # triangular inferior sin diagonal
    F = np.triu(A, k=1)      # triangular superior sin diagonal

    D_inv = np.diag(1 / np.diag(D))

    # Matriz P = D^{-1}(E+F)
    P = -D_inv @ (E + F)
    # Término constante c = D^{-1}b
    c = D_inv @ b

    x = x0.copy()
    for _ in range(max_iter):
        x_new = P @ x + c
        if np.linalg.norm(x_new - x, np.inf) < tol:
            return x_new
        x = x_new

    print("Jacobi no convergió")
    return x

# PARÁMETROS DEL PROBLEMA
L = 1.0             
N = 101          
dx = 1 / (N - 1)
dt = dx**2         

Ta = 0.025      

Pext=(0.472 * 40**2)/(2*L**2)
To=309.65*((0.56)/(Pext * L**2))
r = dt / dx**2
x = np.linspace(0, 2, N)

# CONDICIÓN INICIAL 
phi = To*np.ones(N)

# MATRICES A y B (tridiagonales)
A = np.zeros((N, N))
B = np.zeros((N, N))

for i in range(N):
    # Diagonal principal
    A[i, i] = 1 + r
    B[i, i] = 1 - r

    # Sub y super diagonales
    if i > 0:
        A[i, i-1] = -r / 2
        B[i, i-1] = r / 2
    if i < N-1:
        A[i, i+1] = -r / 2
        B[i, i+1] = r / 2

# BUCLE EN EL TIEMPO 
t = 0.0
while t < Ta:
    b = (B @ phi) + dt  # lado derecho
    phi = jacobi(A, b, phi)  # resolvemos Ax=b por Jacobi
    phi[0] = To
    phi[-1] = To
    t += dt

print("Solución final en T =", Ta)
print(phi)
