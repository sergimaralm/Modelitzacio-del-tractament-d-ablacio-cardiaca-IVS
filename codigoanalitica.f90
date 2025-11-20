program practica1
    implicit none
    REAL, PARAMETER :: pi = 3.1415926535897932384626433832795028842
    REAL, PARAMETER :: beta = 36.5
    REAL, PARAMETER :: tiempo=0.025, x=0.5, To=(0.472*((40/SQRT(2.0))**2))/0.56
    INTEGER, PARAMETER :: numero=9999
    INTEGER :: n
    REAL :: r, T
    n=1
    r=0.0
    DO WHILE (n<numero)
        T = 4*To*SIN(n*pi*x)*((1-EXP(-(n**2)*(pi**2)*(tiempo)))/((n**3)*(pi**3)))
        r = r + T
        n = n + 2
    END DO
        
    WRITE(*,*) 'La Temperatura es:',r+beta
end program practica1