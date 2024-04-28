#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np

def elim_gaussiana(A):
    cant_op = 0
    m = A.shape[0]
    n = A.shape[1]
    Ac = A.copy()
    for i in range(0,n):
        pivote = A[i][i]
        if pivote == 0:
            for j in range(i + 1, n):
                c = A[j][i]
                if c != 0:
                    print("no es posible realizar descomposicion LU sin intercambiar filas")
                else:
                    continue
        for j in range(i + 1, n):
            c = A[j][i]
            if c != 0:
                x = -c / pivote
                A[j] += x * A[i]
                Ac[j]+= x * A[i]
                Ac[j][i]= -x
                cant_op += 1
    
    L = np.tril(Ac, -1) + np.eye(m)
    U = np.triu(Ac)
    
    return L, U, cant_op

def matrices_bn(n:int): 
    n=n
    B = np.eye(n) - np.tril(np.ones((n,n)),-1) 
    B[:n,n-1] = 1

    return B



def probar_descomposicion (n):
    i=2
    res:list =[]
    while i <=n:
        A= matrices_bn(i)
        L, U, cant_op = elim_gaussiana(A)
        if np.allclose(np.linalg.norm(A - L@U, 1), 0):
            res.append("SI")
        else:
            res.append("NO")
      
        i+=1
    
    return res
resultados= probar_descomposicion(100)

def main():
    n = 7
    B = np.eye(n) - np.tril(np.ones((n,n)),-1) 
    B[:n,n-1] = 1
    print('Matriz B \n', B)
    
    L,U,cant_oper = elim_gaussiana(B)
    
    print('Matriz L \n', L)
    print('Matriz U \n', U)
    print('Cantidad de operaciones: ', cant_oper)
    print('B=LU? ' , 'Si!' if np.allclose(np.linalg.norm(B - L@U, 1), 0) else 'No!')
    print('Norma infinito de U: ', np.max(np.sum(np.abs(U), axis=1)) )

if __name__ == "__main__":
    main()
    
    