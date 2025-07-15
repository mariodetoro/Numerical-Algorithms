# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 18:46:15 2023

@author: Mario
"""
import numpy as np


def QR(A):
    n = np.shape(A)[0]
    Q = np.zeros ((n, n))
    R = np.zeros ((n, n))
    for i in range (0, n):
        v = A[:, i]
        s = 0
        for j in range (0, i):
            R[j, i] = np.dot(v, Q[:, j])
            s += (np.dot(v, Q[:, j]) * Q[:,j])
        u = v - s
        u = u / np.linalg.norm(u)
        Q[:, i] = u
        R[i, i] =  np.dot(v, Q[:, i])
    return Q, R
    # return np.linalg.qr(A)
 
def PotenciaDesplazamiento(A,d):
    N = 0
    for i in range(len(A)):
        N += 1 
    dI = d * np.identity(N)
    B = A - dI
    B1 = np.linalg.inv(B)
    N = 0
    for i in range(len(A)):
        N += 1  
    u = np.ones ((N, 1))
    Km = 100
    e = 1e-6
    λ_m = 10
    for k in range(1, Km):
        u = np.matmul(B1, u)
        u = u / np.linalg.norm(u)
        λ = np.matmul(np.transpose(u), np.dot(A, u))
        if abs(λ - λ_m) < e:
            λ_m = λ
    return float(λ), u

def IteracionSimultanea(A):  
    u = np.zeros(A.shape[1])
    Q_0 = np.identity(A.shape[0])
    e = 1e-6
    λ_m = 10
    for k in range(1, 1000):
        Z = np.dot(A, Q_0)
        Q, R = QR(Z)
        Q_0 = Q
        for j in range(Q_0.shape[0]):
            u[j] = np.dot(np.transpose(Q_0[:,j]), np.dot(A, Q_0[:,j]))
            if abs(u[j]- λ_m) < e:
                break
    return u, Q_0

def AlgoritmoQR(A):
    N = np.shape(A)[1]
    autovectores = np.zeros((N,N))
    autovalores = np.zeros((N,1))
    Z = np.copy(A)
    for i in range(1000):
        Q = QR(Z)[0]
        R = QR(Z)[1]
        Z = np.matmul(R,Q)
    for i in range(N):
        for j in range(N):
            autovalores[j] = Z[j][j]
            autovectores[j][i] = PotenciaDesplazamiento(A, autovalores[i])[1][j]
    return autovalores, autovectores
    # return np.linalg.eig(A)
        
        
    

A = np.array([[1.,-1., -2.],[4., 0., 3.],[-1., 1., 1.]])
B = np. array ([[1, 2, -1],
                  [1,0, 1],
                  [4, -4,5]])
C = np. array ([[1, 1, 2],
                  [1,-1, 1],
                  [2, 1,-1]])

# print(QR(B))
print(AlgoritmoQR(C))
print(IteracionSimultanea(C))
# print(PotenciaDesplazamiento(A,2))
    