import numpy as np
import matplotlib.pyplot as plt

from codigoanalitica import analitica
#condicions inicials
L=2
N=101 #mallat
Pext=(0.472 * 40**2)/(2*L**2) #normalització
deltax=1/(N-1)
deltat_list=[0.51*((deltax)**2), 0.49*((deltax)**2), 0.25*((deltax)**2)] #Llista de mallats
To=309.65*((0.56)/(Pext * L**2))
ta=0.025

T_analitic = np.array(analitica(ta))

for deltat in deltat_list: #iteració per als diferents mallats
    rang_temps=round(ta/deltat)
    T =np.zeros(N) #Creació de llista de 0
    T_nuevo=np.zeros(N)
    x = np.linspace(0,2,N) 
    T.fill(To) #Canvi de la llista de 0 per la nostra condició inicial

    for n in range(rang_temps):#iteració del métode euler explícit
        for i in range(1,N-1):
            T_nuevo[i]= T[i] + (deltat)/(deltax)**2 * (T[i+1] - 2*T[i] + T[i-1] + (deltax)**2)
            T_nuevo[0]=To #apliquem condicions de contorn
            T_nuevo[-1]=To
        
        T = T_nuevo.copy() 

    T_graus =  (T/((0.56)/(Pext * L**2)))-273.15 #desnormalitzem la temperatura que tenim
    Error = np.abs(T_graus-T_analitic) #Calcul d'error
    #Creació de les gràfiques de temperatura
    plt.figure()
    plt.plot(x,(T_graus))
    plt.grid(True)
    plt.xlabel('Posició (cm)')
    plt.ylabel('Graus Centígrads')
    plt.hlines(50,0,2,color='r',linestyle='--')
    plt.hlines(80,0,2,color='r',linestyle='--')
    plt.vlines(0.75,36.5,80,color='r',linestyle='--')
    plt.vlines(1.25,36.5,80,color='r',linestyle='--')
    plt.show()
    #Creació de les gràfiques d'errors
    plt.figure()
    plt.plot(x, Error, label=f'Error $\Delta$t={deltat:.6f}')
    plt.title(f'Error Numéric Absolut ($\Delta$t={deltat:.6f})')
    plt.xlabel('Posició (cm)')
    plt.ylabel('Error Absolut (°C)')
    plt.grid(True)
    plt.legend()
    plt.show()
        # if any(((T[0:36]/((0.56)/(Pext * L**2)))-273.15) > 50) == True:
        #         print(f"El temps mínim en segons és de {n*deltat}")
        #         break
    



    # print(f"El temps mínim en segons és de {n*deltat}")