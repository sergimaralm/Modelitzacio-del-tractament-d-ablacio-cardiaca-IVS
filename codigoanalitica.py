import math
beta = 36.5
To = (0.472*((40/(math.sqrt(2)))**2))/(0.56)
def Temperatura (x,t,n):
    T = 4*To*math.sin(n*math.pi*x)*((1-math.exp(-n**2 * math.pi**2 * t))/(n**3 * math.pi**3))
    return T
n=1
N=9999
r=0
while n < N:
    r = r + Temperatura(0.025,0.5,n)
    n = n + 2

print("La Temperatura és")
print(r+beta)
