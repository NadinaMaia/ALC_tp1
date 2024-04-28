#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

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
        for i in range(0,n):
            Mi= np.zeros((n,n))
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
                    cant_op+=1
                    Mi[j][i]= x
                    Ac[j]= Ac[j]+ x*Ac[i]
            L= L - Mi
    ## hasta aqui      
    L = L
    U = Ac
    return L, U, cant_op

A= np.array([[1,2,2],[2,1,1],[2,2,1]])
B = np.array([[1, 2, 3, 4],
    [5, 6, 7, 8],
    [9, 10, 11, 12],
    [13, 14, 15, 16]])


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

def dataGrafico( n: int):
    i=2
    tamaño:list=[]
    operaciones:list=[]
    while i <=n:
        A= matrices_bn(i)
        L, U, cant_op = elim_gaussiana(A)
        tamaño.append(i)
        operaciones.append(cant_op)
        i+=1
    data = {'tamaño': tamaño , 'operaciones': operaciones}
    return data
        
data = dataGrafico(100)
df = pd.DataFrame(data)
plt.bar(data= df, x='tamaño', height='operaciones')
# Genera el grafico de cantidad de operaciones de la elimincacion gausiana en funcion del tamaño de una mtriz
fig, ax = plt.subplots()    
plt.rcParams['font.family'] = 'sans-serif'                
ax.bar(data=data, 
       x='tamaño', 
       height='operaciones', 
       color='hotpink'
      )
ax.set_title('cantidad de operaciones sugun el tamaño de la matriz ')                    # Agrega el titulo en el margen superior
ax.set_xlabel('tamaño', fontsize='medium')                       # Agrega una etiqueta sobre el eje x
ax.set_ylabel('Cantidad de operaciones', fontsize='medium')    # Agrega una etiqueta sobre el eje y
ax.set_ylim(0, 5000)                                         # Asigna rango al eje y

 
plt.scatter(data = df, x='tamaño', y='operaciones')

fig, ax = plt.subplots()
 
plt.rcParams['font.family'] = 'sans-serif'           
ax.scatter(data = data, 
             x='tamaño', 
             y='operaciones',
             s=5,                         # Tamano de los puntos
            color= 'indigo')

ax.set_title('cantidad de operaciones sugun el tamaño de la matriz ')                    # Agrega el titulo en el margen superior
ax.set_xlabel('tamaño', fontsize='medium')                       # Agrega una etiqueta sobre el eje x
ax.set_ylabel('Cantidad de operaciones', fontsize='medium')    # Agrega una etiqueta sobre el eje y
ax.set_ylim(0, 5000)  
ax.set_xlim(0,100)



 
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
    
    