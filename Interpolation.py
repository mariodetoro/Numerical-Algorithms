# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 18:32:59 2023

@author: mdeto
"""

import numpy as np
import sympy as sym
from scipy.interpolate import lagrange
from Tema_I_Sistemas_de_Ecuaciones import Gauss

def Lagrange(x,x_i,y_i):
    n = np.shape(x_i)[0]
    px = 0
    for i in range(n):
        L = 1
        for j in range(n):
            if i != j:
                L *= (x - x_i[j]) / (x_i[i] - x_i[j])
        px += L*y_i[i]
    if x == sym.Symbol('x'):
        return px.expand()
    else:
        return px
    
def MinimosCuadrados(x, y):
    n = np.shape(x)[0]
    xy = 0
    x2 = 0
    for i in range(n):
        xy += x[i] * y[i]
        x2 += x[i] ** 2
    x_m = sum(x) / n
    y_m = sum(y) / n
    sxy = xy - n*x_m*y_m
    sxx = x2 - n * (x_m**2)
    m = sxy / sxx
    c = y_m - m * x_m
    return m, c

def MinimosVandermonde(x, y):
    N = np.shape(x)[0]
    K = np.ones((N, N))
    for i in range(1,N):
        for j in range(N):
            K[j][i] = x[j] ** i
    # return np.linalg.lstsq(K, y, rcond= None)[0]
    return Gauss(K,y)
    

def RegresiónPolinómica(x, y, n):
    return np.polyfit (x, y, n)

x_i = np.array ([1, 2.5, 4])
y_i = np.array ([2, 7, 3])
poly = lagrange(x_i,y_i)
X = ([0,1,2,3])
Y = ([-1,0.2,0.9,2.1])



x = sym.Symbol('x')
print("------LAGRANGE------")
print("Punto/Polinomio =",Lagrange(x,x_i, y_i))
print("Polinomio =",lagrange(x_i, y_i))
# print(poly)
print("------MINIMOS CUADRADOS------")
print("m =",MinimosCuadrados(X, Y)[0])
print("c =",MinimosCuadrados(X, Y)[1])
print("------VANDERMONDE------")
print(MinimosVandermonde(X, Y))
print("------REGRESION POLINOMICA------")
print(RegresiónPolinómica(X, Y, 2))