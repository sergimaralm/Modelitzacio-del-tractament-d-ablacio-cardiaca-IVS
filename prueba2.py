import numpy as np

def jacobi(A, b, x0, tol=1e-10, max_iter=5000):
    D = np.diag(np.diag(A))
    E = np.tril(A, k=-1)
    F = np.triu(A, k=1)
    D_inv = np.diag(1 / np.diag(D))

    P = -D_inv @ (E + F)
    c = D_inv @ b

    x = x0.copy()
    for _ in range(max_iter):
        x_new = P @ x + c
        if np.linalg.norm(x_new - x, np.inf) < tol:
            return x_new
        x = x_new

    print("Jacobi no convergió")
    return x


# PARÁMETROS
L = 1.0
N = 101
dx = 1/(N-1)
dt = dx**2
Ta = 0.025

Pext = (0.472 * 40**2)/(2*L**2)
To = 309.65*((0.56)/(Pext * L**2))

r = dt/dx**2
x = np.linspace(0, L, N)

# CONDICION INICIAL
phi = To*np.ones(N)

# MATRICES REDUCIDAS (solo nodos interiores)
M = N-2
A = np.zeros((M, M))
B = np.zeros((M, M))

for i in range(M):
    A[i, i] = 1 + r
    B[i, i] = 1 - r

    if i > 0:
        A[i, i-1] = -r/2
        B[i, i-1] = r/2
    if i < M-1:
        A[i, i+1] = -r/2
        B[i, i+1] = r/2


# BUCLE EN EL TIEMPO
t = 0.0
while t < Ta:

    # Vector interior
    phi_int = phi[1:-1]

    # Lado derecho reducido (interiores)
    b = B @ phi_int + dt * np.ones(M)

    # CORRECCIÓN: añadir la contribución completa de las fronteras
    # (r/2 desde B*T^n no está incluida porque B es MxM; hay que sumar
    # tanto la parte de B*T^n como la parte que viene de mover A*T^{n+1})
    b[0]  += r * To
    b[-1] += r * To

    # Resolver por Jacobi
    phi_int = jacobi(A, b, phi_int)

    # Ensamblar solución completa
    phi[0] = To
    phi[-1] = To
    phi[1:-1] = phi_int

    t += dt


print("Solución final en T =", Ta)
print(phi)
