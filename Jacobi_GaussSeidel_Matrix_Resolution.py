# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 18:33:49 2023

@author: Mario
"""

import numpy as np
from math import sqrt

def LDU(A):
    N = np.shape(A)[0]
    L = np.zeros((N,N))
    D = np.zeros((N,N))
    U = np.zeros((N,N))
    Ac = np.copy(A)
    for i in range(N):
        D[i][i] = Ac[i][i]
        for j in range(i,N):
           L[j][i] = A[j][i]
    L = L -D
    U = Ac - D - L
    return L, D, U
        
def Jacobi(A, B):
    A = np.array(A, float)
    B = np.array(B, float)
    N = int(sqrt(np.size(A)))
    X = np.ones((N,1), float)
    X1 = np.ones((N,1), float)
    for k in range(0, 1000):
        for i in range(0, N):
            t = 0
            for j in range(0, N):
                if i != j:
                    t += A[i][j] * X[j]
            X1[i] = (B[i] - t) / A[i][i]
        X = X1
    return X1

def JACOBI(A,B):
    N = np.shape(A)[0]
    X = np.ones((N,1), float)
    Km = 100
    e = 1e-12
    L,D,U = LDU(A)
    T = np.matmul(-np.linalg.inv(D),(L + U))
    C = np.matmul(np.linalg.inv(D),B)
    n = 0
    for k in range(Km):
        X_OLD = X
        X = np.matmul(T,X_OLD) + C
        if (np.linalg.norm((X - X_OLD)) / np.linalg.norm(X)) < e:
            return X,n
        n += 1
    return "No converge"

def GAUSS_SEIDEL(A,B):
    N = np.shape(A)[0]
    X = np.ones((N,1), float)
    Km = 100
    e = 1e-12
    L,D,U = LDU(A)
    T = np.matmul(-np.linalg.inv(L + D),U)
    C = np.matmul(np.linalg.inv(L + D),B)
    n = 0
    for k in range(Km):
        X_OLD = X
        X = np.matmul(T,X_OLD) + C
        if (np.linalg.norm((X - X_OLD)) / np.linalg.norm(X)) < e:
            return X, n
        n += 1
    return "No converge"

def Gauss_Seidel(A, B):
    A = np.array(A, float)
    B = np.array(B, float)
    N = int(sqrt(np.size(A)))
    X = np.ones((N,1), float)
    for k in range(0, 1000):
        for i in range(0, N):
            t = 0
            for j in range(0, N):
                if i != j:
                    t += A[i][j] * X[j]
            X[i] = (B[i] - t) / A[i][i]
    return X
            
        
        
A = np.array ([[4, -1, 1],
              [2, 5, 2],
              [1, 2 ,4]])
B = np.array ([[8],[3],[11]])

# print(LDU(A))
print(JACOBI(A,B))
# print(Jacobi(A, B))
print(GAUSS_SEIDEL(A,B))
# print(Gauss_Seidel(A, B))
        