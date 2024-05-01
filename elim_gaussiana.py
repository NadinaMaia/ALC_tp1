#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np

def elim_gaussiana(A):
    cant_op = 0
    m=A.shape[0]
    n=A.shape[1]
    Ac = A.copy()
    L= np.eye(m)
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    ## desde aqui -- CODIGO A COMPLETAR
    if m!=n: 
        print('Matriz no cuadrada')
        return
    else: 
        for i in range (0, (n)): 
            m=i 
            l=1 
            pivote= Ac[m][i] 
            Mi= np.zeros((n,n))
            while m+l < n and pivote!=0:
                c= Ac[m+l][i] 
                if c <0 and pivote > 0:
                    x= c/pivote
                elif c>0 and pivote>0: 
                    x= -c/pivote
                elif c<0 and pivote<0:
                    x= -c/pivote
                else  :
                    x= c/pivote 
                if x!=0 :
                    cant_op+=1
                    Mi[m+l][i]= x
                Ac[m+l]= Ac[m+l]+ x*Ac[m]
                l+=1
            L= L-Mi
    
    ## hasta aqui
            
    L = L
    U = np.triu(Ac)
    
    return L, U, cant_op

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
    
    