# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 22:28:55 2023

@author: mdeto
"""
import numpy as np


def Gauss(A, B):
    A = np.array(A, float)
    B = np.array(B, float)
    N = 0
    for i in range(len(A)):
        N += 1     
    for k in range(0, N - 1):
        for i in range (k + 1, N):
            for t in range(k + 1, N):
                if (A[t][k] > A[k][k]):
                    S = np.copy(A[t,:])
                    A[t,:] = A[k,:]
                    A[k,:] = S
                    s = np.copy(B[t])
                    B[t] = B[k]
                    B[k] = s
            if (abs(A[k][k]) < (10**(-14))):
                print ("No solution")
            l = (A[i][k]) / (A[k][k])
            for j in range (k, N):
                A[i][j] = A[i][j] - l*A[k][j]
            B[i] = B[i] - l*B[k]      
    x = np.ones((N, 1), float)
    N = N - 1
    x[N] = (1 / (A[N][N])) * B[N]


    for i in range(N - 1, -1, -1):
        F = 0
        for j in range(i, N):
            F += A[i][j + 1]*x[j + 1]
        x[i] = ( (B[i] - F) / (A[i][i]) ) 
    return x

def LU (A):
    A = np.array(A, float)
    N = 0
    for i in range(len(A)):
        N += 1  
    L = np.identity(N)
    for k in range(0, N - 1):
        for i in range (k + 1, N):
            # for t in range(k + 1, N):
            #     if (A[t][k] > A[k][k]):
            #         S = np.copy(A[t,:])
            #         A[t,:] = A[k,:]
            #         A[k,:] = S
            #         s = np.copy(B[t])
            #         B[t] = B[k]
            #         B[k] = s
            # if (abs(A[k][k]) < (10**(-14))):
            #     print ("No solution")
            l = (A[i][k]) / (A[k][k])
            L[i][k] = l
            for j in range (k, N):
                A[i][j] = A[i][j] - l*A[k][j]
    U = A
    return L , U


def Det(A):
    U = LU(A)[1]
    N = 0
    D = 1
    for i in range(len(U)):
        N += 1  
    for k in range(0, N):
        D = D * U[k][k]
    return D

def Inversa(A):  
    A = np.array(A, float)
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

def Sistema(A, B):
    N = 0
    for i in range(len(A)):
        N += 1  
    L = LU(A)[0]
    U = LU(A)[1]
    Y = np.zeros((N, 1))
    Y[0] = B[0]
    for i in range(0, N - 1):
        F = 0
        for j in range(0, N - 1):
            F += L[i + 1][j] * Y[j]
        Y[i + 1] = B[i + 1] - F
    x = np.zeros((N, 1))
    N = N - 1
    x[N] = (1 / (U[N][N])) * Y[N]
    for i in range(N - 1, -1, -1):
        F = 0
        for j in range(i, N):
            F += U[i][j + 1]*x[j + 1]
        x[i] = (1 / (U[i][i])) * (Y[i] - F)
    return x

def Norma(A):
    AN = np.copy(A)
    fila = []
    columna = []
    N = 0
    t = 0
    s = 0
    for i in range(len(AN)):
        N += 1  
    for i in range(0, N):
        for j in range(N):
            t+= AN[i][j]
            s+= AN[j][i]
        fila.append(t)
        columna.append(s)
        t = 0
        s = 0
    return max(fila), max(columna)

def KA(A):
    if (Norma(A)[1] * Norma(Inversa(A))[1]) <= 50:
        return True
    else:
        return False

def Vandermonde(A, B):
    N = 0
    for i in range(len(A)):
        N += 1 
    K = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            K[i][j] = A[i][0] ** j
    if (KA(K)):
        X = Sistema(K, B)
        return X
    else:
        return "Not valid, K(A) =",(Norma(A)[1] * Norma(Inversa(A))[1])
    

A = np.array ([[2,3,4,5],
             [6,15,19,23],
             [8,42,60,70],
             [12,60,1,17]])
B = np.array ([[5],[30],[98],[144]])
C = np.array ([[0,2,1],
             [2,2,4],
             [-1,0.5,-1]])
D = np.array ([[1],[-2],[0]])
E = np.array ([[10,-1,2,0],
             [-1,11,-1,3],
             [2,-1,10,-1],
             [0,3,-1,8]])
F = np.array ([[6],[25],[-11],[15]])
G = np.array ([[2,1,0],
             [1,2,1],
             [0,1,2]])
H = np.array ([[1],[2],[3]])
I = np.array(([[1], [2], [3]]))
J = np.array(([[3], [5], [7]]))
K = np.array ([[3, 1, 0], 
               [1, 2, -1], 
               [0, -1, 1]])
L = np.array ([[1,1,1,1],
             [2,-1,3,-4],
             [3,2,-1,5],
             [1,-3,2,-4]])
M = np.array ([[10],[9],[13],[-3]])
n = np. array ([[1, 2, -1],
                  [1,0, 1],
                  [4, -4,5]])

# print(Gauss(A,B))
# print(LU(A))
# print(Det(A))
# print(Inversa(A))
# print(Sistema(A, B))
# print(Norma(A))
# print(Vandermonde(I,J))
