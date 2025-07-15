# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 19:27:24 2023

@author: Mario
"""
import numpy as np

def Unitario(u):
    s = 0
    for item in u:
        s = s + (item)**2
    norma = s**0.5
    return u / norma

def LU (A):
    A = np.array(A, float)
    N = 0
    for i in range(len(A)):
        N += 1  
    L = np.identity(N)
    for k in range(0, N - 1):
        for i in range (k + 1, N):
            l = (A[i][k]) / (A[k][k])
            L[i][k] = l
            for j in range (k, N):
                A[i][j] = A[i][j] - l*A[k][j]
    U = A
    return L , U

def Inversa(A):  
    L = LU(A)[0]  
    U = LU(A)[1]
    N = 0
    for i in range(len(U)):
        N += 1  
    I = np.identity(N)
    Y = np.zeros((N, N))
    S = np.zeros((N, 1))
    for k in range(0, N):
        S[0] = I[0][k]
        B = I[:,k]
        for i in range(0, N - 1):
            F = 0
            for j in range(0, N - 1):
                F += L[i + 1][j] * S[j]
            S[i + 1] = B[i + 1] - F
        for t in range(0, N):
            Y[t][k] = S[t]
        S = np.zeros((N, 1))
    A1 = np.zeros((N,N))
    x = np.zeros((N, 1))
    T = np.zeros((N,1))
    N = N - 1
    for k in range(0, N + 1):
        T = Y[:,k]   
        x[N] = (1 / (U[N][N])) * T[N]
        for i in range(N - 1, -1, -1):
            F = 0
            for j in range(i, N):
                F += U[i][j + 1]*x[j + 1]
            x[i] = (1 / (U[i][i])) * (T[i] - F)
        for t in range(0, N + 1):
            A1[t][k] = x[t]
    return A1

def Potencia(A):
    N = 0
    for i in range(len(A)):
        N += 1  
    u = np.ones ((N, 1))
    Km = 100
    e = 1e-6
    λ_m = 10
    for k in range(1, Km):
        u = np.matmul(A, u)
        u = Unitario(u)
        λ = np.matmul(np.transpose(u), np.dot(A, u))
        if abs(λ - λ_m) < e:
            λ_m = λ
    return float(λ), u

def PotenciaInversa(A):
    A1 = Inversa(A)
    N = 0
    for i in range(len(A)):
        N += 1  
    u = np.ones ((N, 1))
    Km = 100
    e = 1e-6
    λ_m = 10
    for k in range(1, Km):
        u = np.matmul(A1, u)
        u = Unitario(u)
        λ = np.matmul(np.transpose(u), np.dot(A, u))
        if abs(λ - λ_m) < e:
            λ_m = λ
    return float(λ), u

def PotenciaDesplazamiento(A):
    N = 0
    for i in range(len(A)):
        N += 1 
    dI = -10 * np.identity(N)
    B = A - dI
    B1 = Inversa(B)
    N = 0
    for i in range(len(A)):
        N += 1  
    u = np.ones ((N, 1))
    Km = 100
    e = 1e-6
    λ_m = 10
    for k in range(1, Km):
        u = np.matmul(B1, u)
        u = Unitario(u)
        λ = np.matmul(np.transpose(u), np.dot(A, u))
        if abs(λ - λ_m) < e:
            λ_m = λ
    return float(λ), u
        
A = np. array ([[1, 2, -1],
                  [1,0, 1],
                  [4, -4,5]])
B = np. array ([[1, 1, 2],
                  [1,-1, 1],
                  [2, 1,-1]])
C = np. array ([[3, 2],
                  [7,-2]])
D = np. array ([[3, -1, 1],
                  [-2,4, 1],
                  [2, -1,2]])
E = np. array ([[2, -12],
                  [1,-5]])
F = np. array ([[1, 2, 3],
                  [4,5, 6],
                  [7, 8,9]])
G = np.array ([[3, 1, 0], [1, 2, -1], [0, -1, 1]])
H = np.array ([[1,1,1,1],
             [2,-1,3,-4],
             [3,2,-1,5],
             [1,-3,2,-4]])

print(Potencia(B))
print (PotenciaInversa(B))
print(PotenciaDesplazamiento(B))


            