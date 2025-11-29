import math
import numpy as np
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    x_list = np.linspace(0,2,101)
    T_list = analitica(0.025)

    plt.plot(x_list,T_list)
    plt.grid()
    plt.xlabel('Posició (cm)')
    plt.ylabel('Graus Centígrads')
    plt.hlines(50,0,2,color='r',linestyle='--')
    plt.hlines(80,0,2,color='r',linestyle='--')
    plt.vlines(0.75,36.5,80,color='r',linestyle='--')
    plt.vlines(1.25,36.5,80,color='r',linestyle='--')
    print("La Temperatura és")
    print(T_list)
    plt.show()