#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np

def elim_gaussiana(A):
    cant_op = 0
    m=A.shape[0] #filas
    n=A.shape[1] #columnas
    Ac= A.copy() 
    
    if m!=n:
        print('Matriz no cuadrada')
        return 
    
    ## desde aqui -- CODIGO A COMPLETAR
    
    ## Aca calculo los valores de los multiplicadores y lo actualizo en A
    for j in range (n): #columnas
        pivote = A[j,j]
        for i in range(1, m): #filas
            if j < i:
                k = calculo_k(A[i], pivote, j) #calculo k
                if k != 0:
                    A[i] = A[i] - k*A[j]
                    Ac[i][j] = k #agrego k 
                    cant_op += 1
                else:
                    continue
                
    ## Aca actualizo los valores arriba de la diagonal
    for j in range (n): #columnas
        for i in range(1, m): #filas
            if j >= i:
                Ac[i][j] = A[i][j]
    
    ## hasta aqui

      
    L = np.tril(Ac,-1) + np.eye(A.shape[0]) 
    U = np.triu(Ac)
    return L, U, cant_op
    

def calculo_k(fila_actual, divisor, iterador):
    multiplicador = 0
    if divisor != 0:
        multiplicador = fila_actual[iterador] / divisor
    return multiplicador
   

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
    
    