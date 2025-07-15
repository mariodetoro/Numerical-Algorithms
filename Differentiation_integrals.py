# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 18:55:05 2023

@author: mdeto
"""

import numpy as np
import scipy as sp
from math import pi, sin, cos
import matplotlib.pyplot as plt

def DerivadaProgresiva (f, x):
    h = 1e-3
    return (f(x + h) - f(x))/ h

def DerivadaRegresiva (f, x):
    h = 1e-3
    return(f(x) - f(x - h))/ h

def DerivadaCentrada (f, x):
    h = 1e-3
    return (f(x + h) - f(x - h)) / (2*h)

def Derivada2Centrada (f, x):
    h = 1e-3
    return (1/(h**2))*(f(x - h) -2*f(x) +f(x + h))

def DerivadaCentradaProgresiva3 (f, x):
    h = 1e-3
    return (1/(2*h))*(-3*f(x) +4*f(x + h) -f(x + 2*h))
    
def DerivacionNumerica(f, x1, x2, N):
    global puntos 
    puntos= np.linspace(x1, x2, N)  #Creamos una partición equiespaciada de puntos entre los extremos dados
    D = np.zeros ((N,N))
    dx = (x2 - x1)/ (N - 1)         #Calculamos la distancia entre los puntos
    for i in range(N):              #Hallamos la imagen de los puntos de la partición
        puntos[i] = f(puntos[i])
    D[0][0] = -1
    D[0][1] = 1
    for i in range(1,N -1):         #Creamos la matriz de derivación numérica
        for j in range(N -1):
            D[i][i - 1] = -1/2
            D[i][i] = 0
            D[i][i + 1] = 1 / 2
    D[N - 1][N - 2] = -1
    D[N - 1][N - 1] = 1
    D = 1/dx * D                  
    return np.matmul(D, puntos)
    
def f(x):
    return x**2

def IntegralRiemann(f, a, b, n):
    dx = (b - a) / (n - 1)                              #Calculamos la distancia entre los puntos
    integral = 0
    x_values = np.linspace(a + dx/2, b - dx/2, n - 1)   #Creamos una partición equiespaciada de puntos entre los extremos dados.  Para tomar el valor medio de cada trapecio, sumamos y restamos la mitad de la distancia al inicio y al final
    # xS_values = np.linspace(a, b - dx, n - 1)
    integral = dx*f(x_values).sum()                     #Multiplicamos la distancia entre valores por la suma de todos los puntos calculados    
    return integral

def TrapecioIntegral(f, a, b, n):
    dx = (b - a) / (n - 1)                                                             #Calculamos la distancia entre los puntos
    x_values = np.linspace(a, b, n)                                                    #Creamos una partición equiespaciada de puntos entre los extremos dados.  
    integral = dx*(f(x_values[1:-1]).sum() + f(x_values[0])/2 +f(x_values[-1])/2)      #Multiplicamos la distancia entre puntos por la suma de las imágenes de los mismos, excepto los extremos, que van divididos entre dos.
    return integral

def SimpsonIntegral(f, a, b, n):
    dx = (b - a) / (n - 1)                                                                                                    #Calculamos la distancia entre los puntos
    x_values = np.linspace(a, b, n)                                                                                           #Creamos una partición equiespaciada de puntos entre los extremos dado
    integral = (dx/3)*(f(x_values[0]) + 2*(f(x_values[2:-1:2])).sum() + 4*(f(x_values[1:-1:2])).sum() + f(x_values[-1]))      #Multiplicamos la distancia entre puntos dividida entre tres a las imágenes de los puntos, con la excepción de los puntos pares, que van multiplicados por 2; y los impares, por 4.
    return integral

def F(x):
    for i in range(x.shape[0]):
        x[i] = sin(x[i])
    return x


x = np.linspace(-1, 1, 21)
y = x** 2
        
print("Derivada Progresiva =",DerivadaProgresiva(sin, pi/3))
print("Derivada Regresiva =",DerivadaRegresiva(sin, pi/3))
print("Derivada Centrada =",DerivadaCentrada(sin, pi/3))
print("Derivada2Centrada =",Derivada2Centrada(sin,pi/2))
print("DerivadaCentradaProgresiva3 =",DerivadaCentradaProgresiva3(sin,pi/2))
print(DerivacionNumerica(sin, 0, pi, 11))
print("Riemann = ",IntegralRiemann(f, -1, 1, 21))
print("Trapecio =",TrapecioIntegral(f, -1, 1, 21))
print("Simpson =",SimpsonIntegral(f, -1, 1, 21))
print("TrapecioScipy = ",sp.integrate.trapezoid(y,x))
print("SimpsonScipy =",sp.integrate.trapezoid(y,x))

# ax = plt.plot(np.linspace(0,pi,30),f(np.linspace(0,pi,30)),color="red")
# bx = plt.plot(np.linspace(0,pi,30),DerivacionNumerica(sin,0,pi,30),color="orange")
# plt.xlabel('x', fontsize=12)
# plt.ylabel('y', fontsize=12)
# plt.title('Derivación Numérica de puntos de una función', fontsize=16)

