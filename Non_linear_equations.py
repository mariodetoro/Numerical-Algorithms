# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 18:38:48 2023

@author: Mario
"""
from math import cos, exp, sin
from scipy.optimize import fsolve

def Bolzano (f,a,b):
    xi = a                                                                     #Definimos nuestros valores extremos
    x = b
    while abs((x - xi)/x) > 1e-14  or (abs(f(x))) > (1e-14):                   #Establecemos los criterios de convergencia
        xi = x
        x = (a + b)/2                                                          #Calculamos el punto medio entre los valores
        if f(x)*f(a) < 0:                                                      #Si el producto es negativo/positivo, sustiuimos nuestros extremos y seguimos operando
            b = x
        elif f(x)*f(b) < 0:
            a = x
        elif f(x) == 0:                                                        #Si la imagen del punto medio es 0, hemos encontrado nuestra solución.
            return x
    return x

def NewtonRaphson (f,x):
    xi = 1
    while abs((x - xi)/x) > 1e-14  or (abs(f(x))) > (1e-14):                   #Establecemos los criterios de convergencia
        xi = x
        x = xi - (f(xi)/fd(xi))                                                #Hallamos el punto de corte con el eje x de la recta tangente a la imagen del punto dado
    return x

def Secante(f,x1, x2):
    while abs((x2 - x1)/x2) > 1e-14  or (abs(f(x2))) > (1e-14):                #Establecemos los criterios de convergencia
        x_old = x2
        x2 = x2 - f(x2)*((x1 - x2)/(f(x1) - f(x2)))                            #Hallamos el punto de corte con el eje x de la recta que pasa por las imágenes de los puntos dados
        x1 = x_old
    return x2

def f(x):
    return (exp(3*x) - cos(x))
def fd(x):
    return (3*exp(3*x) + sin(x))

print("Bolzano =",Bolzano(f,-1,-3))
print("NewtonRaphson =",NewtonRaphson(f,-1))
print("Secante =",Secante(f,-1,-3))
print("Raíz exacta en el punto =", fsolve(f,-1))


    