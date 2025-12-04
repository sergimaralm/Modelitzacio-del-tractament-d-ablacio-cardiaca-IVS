import math
import numpy as np
import matplotlib.pyplot as plt

#Retorna temperatures en el temps ta
def analitica(ta):
    def Temperatura (x,t,n):
        To = (0.472*((40/(math.sqrt(2)))**2))/(0.56)
        T = 4*To*math.sin(n*math.pi*x)*((1-math.exp(-n**2 * math.pi**2 * t))/(n**3 * math.pi**3))
        return T
    
    beta = 36.5
    n=1
    r=0
    N=9999    
    
    x_list = np.linspace(0,1,101)
    T_list = []

    for x in x_list:
        while n < N:
            r = r + Temperatura(x,ta,n)
            n = n + 2
        T_list.append(r+beta)
        n=1
        r=0
    return T_list

# Metóde Jacobi resoldre Ax = b, donat x0
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