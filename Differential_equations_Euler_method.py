# -*- coding: utf-8 -*-
"""
Created on Thu May 11 18:32:16 2023

@author: Mario de Toro y Jorge de Torres
"""
import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos

# --------- CONDICIONES DEL PROBLEMA ---------
m = 1
k = 1
h_1 = 0.02
tf = 5
y0 = [0,0]
y1 = [0,-1]

# --------- FUNCIÓN PARA REALIZAR EL MÉTODO DE EULER ---------

def Euler(f,y0,tf,h,c):
    n = int(tf / h)                                     #Calculamos el número de intervalos de la partición
    y = np.zeros((2,n))                                 #Creamos la matriz en la que vamos a almacenar los valores de la posición y velocidad
    t = np.linspace(0,tf,n)                             #Creamos la partición del intervalo de tiempo
    y[:,0] = y0                                         #Establecemos las condiciones iniciales
    for i in range(n-1):
        y[:,i + 1] = y[:,i] + h*(f((y[:,i]),(t[i]),c))    #Aplicamos el método de Euler recurriendo a la función externa
    return t, y[0],y[1]


# 1) MATRIZ DE LA FUNCIÓN A EVALUAR

def f(u,t,c):
    F = 1
    A = np.array([[-(2*c)/m,-k/m],[1,0]])
    B = np.array([F/m,0])
    S = np.matmul(A,u) + B
    return S

# --------- SOLUCIÓN ANALÍTICA A F = SIN(T)) ---------
    
def fs(x):
    for i in range(np.shape(x)[0]):
        x[i] = -cos(x[i])
    return x

# 2) CONSTRUCCIÓN DE LA GRÁFICA PARA EL EJERCICIO 2

plt.title("∆t = 2*10-2      F = 1", fontsize=16)
plt.plot(Euler(f,y0,tf,h_1,0.25)[0],Euler(f,y0,tf,h_1,0.25)[2], label = "Position",color = "orange")
plt.plot(Euler(f,y0,tf,h_1,0.5)[0],Euler(f,y0,tf,h_1,0.5)[2],color = "orange")
plt.plot(Euler(f,y0,tf,h_1,0.75)[0],Euler(f,y0,tf,h_1,0.75)[2],color = "orange")
plt.plot(Euler(f,y0,tf,h_1,1.0)[0],Euler(f,y0,tf,h_1,1.0)[2],color = "orange")
plt.plot(Euler(f,y0,tf,h_1,1.25)[0],Euler(f,y0,tf,h_1,1.25)[2],color = "orange")
plt.plot(Euler(f,y0,tf,h_1,1.5)[0],Euler(f,y0,tf,h_1,1.5)[2],color = "orange")
plt.plot(Euler(f,y0,tf,h_1,0.25)[0],Euler(f,y0,tf,h_1,0.25)[1], label = "Speed",color = "blue")
plt.plot(Euler(f,y0,tf,h_1,0.5)[0],Euler(f,y0,tf,h_1,0.5)[1],color = "blue")
plt.plot(Euler(f,y0,tf,h_1,0.75)[0],Euler(f,y0,tf,h_1,0.75)[1],color = "blue")
plt.plot(Euler(f,y0,tf,h_1,1.0)[0],Euler(f,y0,tf,h_1,1.0)[1],color = "blue")
plt.plot(Euler(f,y0,tf,h_1,1.25)[0],Euler(f,y0,tf,h_1,1.25)[1],color = "blue")
plt.plot(Euler(f,y0,tf,h_1,1.5)[0],Euler(f,y0,tf,h_1,1.5)[1],color = "blue")

# 3) CONSTRUCCIÓN DE LA GRÁFICA PARA EL EJERCICIO 3

# plt.title("C = 0.5    ∆t = 10-4   F = sen(t)", fontsize=16)
# plt.plot(Euler(f,y1,tf,h_1,0.5)[0],Euler(f,y0,tf,h_1,0.5)[2], label = "Position")
# plt.plot(Euler(f,y1,tf,h_1,0.5)[0],Euler(f,y0,tf,h_1,0.5)[1], label = "Speed")

#---------  GRÁFICA SOLUCIÓN ANALÍTICA A F = SIN(T)) ---------

# cx = plt.plot(Euler(f,y0,tf,h_1)[0],fs(Euler(f,y0,tf,h_1)[0]),color="red")
plt.xlabel('Time (t)', fontsize=12)
plt.ylabel('Distance (x)', fontsize=12)
plt.legend(loc="upper left")

