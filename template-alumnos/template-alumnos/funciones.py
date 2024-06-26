# -*- coding: utf-8 -*-
"""
Created on Tue May  7 20:17:06 2024

autores: Natali Biasoni, Nadina Soler

descripcion: En este archivo se encuentran las funciones para 
                el TP1 de ALC dividido en funciones.
"""
import numpy as np
import networkx as nx
import scipy
import matplotlib.pyplot as plt
import time
import random

# =============================================================================
# FUNCION PARA LEER LOS ARCHIVOS
# =============================================================================

def leer_archivo(input_file_path):

    f = open(input_file_path, 'r')
    n = int(f.readline())
    m = int(f.readline())
    W = np.zeros(shape=(n,n))
    for _ in range(m):
    	line = f.readline()
    	i = int(line.split()[0]) - 1
    	j = int(line.split()[1]) - 1
    	W[j,i] = 1.0
    f.close()
    
    return W

# =============================================================================
# IMPORTAMOS LOS ARCHIVOS
# =============================================================================
carpeta = "./tests/"
I= leer_archivo(carpeta + "instagram_famosos_grafo.txt") 
M= leer_archivo(carpeta + "mathworld_grafo.txt") 
DS= leer_archivo(carpeta + "test_dosestrellas.txt") 
A= leer_archivo(carpeta+ "test_aleatorio.txt")
TS= leer_archivo(carpeta + "test_30_segundos.txt")

# =============================================================================
# FUNCIONES PARA DIBUJAR GRAFO
# =============================================================================

def dibujarGrafo(W, print_ejes=True):
    
    options = {
    'node_color': 'yellow',
    'node_size': 200,
    'width': 3,
    'arrowstyle': '-|>',
    'arrowsize': 10,
    'with_labels' : True}
    
    N = W.shape[0]
    G = nx.DiGraph(W.T)
    
    #renombro nodos de 1 a N
    G = nx.relabel_nodes(G, {i:i+1 for i in range(N)})
    if print_ejes:
        print('Ejes: ', [e for e in G.edges])
    
    nx.draw(G, pos=nx.spring_layout(G), **options)

# =============================================================================
# FUNCIONES AUXILIARES AL CALCULO DEL RANKING
# =============================================================================

def calculo_D(A):
    m=A.shape[0] #filas
    n=A.shape[1] #columnas
    D = np.zeros((m,n))
    for i in range(n): #columnas
        grado = 0
        for j in range(m): #filas
            grado += A[j][i]
        if grado != 0:
            D[i][i] = 1/grado
    return D

def calculo_A(W,p):
    D = calculo_D(W)
    m = W.shape[0] #filas
    I = np.diag(np.ones(m),0)
    A = I - p*W@D
    return A

def elim_gaussiana(A):
    cant_op = 0
    m=A.shape[0] #filas
    n=A.shape[1] #columnas
    Ac= A.copy() 
    Ad =  A.copy() 
    
    if m!=n:
        print('Matriz no cuadrada')
        return 

    ## Aca calculo los valores de los multiplicadores y lo actualizo en A
    for j in range (n): #columnas
        pivote = Ad[j,j]
        for i in range(1, m): #filas
            if j < i:
                k = calculo_k(Ad[i], pivote, j) #calculo k
                if k != 0:
                    Ad[i] = Ad[i] - k*Ad[j]
                    Ac[i][j] = k #agrego k 
                    cant_op += 1
                else:
                    continue
                
    ## Aca actualizo los valores arriba de la diagonal
    for j in range (n): #columnas
        for i in range(1, m): #filas
            if j >= i:
                Ac[i][j] = Ad[i][j]
     
    L = np.tril(Ac,-1) + np.eye(A.shape[0]) 
    U = np.triu(Ac)
    return L, U, cant_op
    

def calculo_k(fila_actual, divisor, iterador):
    multiplicador = 0
    if divisor != 0:
        multiplicador = fila_actual[iterador] / divisor
    return multiplicador
   
def ranking(score):
   sorted_indices = sorted(range(len(score)), key=lambda i: score[i], reverse=True)
   
   # Crear una lista para almacenar el ranking
   rnk = [0] * len(score)
   
   rank = 1
   for i in range(len(score)):
        if i > 0 and (np.isclose(score[sorted_indices[i]],score[sorted_indices[i-1]]))== False:
            rank += + 1
        rnk[sorted_indices[i]] = rank
    
   return rnk
  

            
# =============================================================================
# FUNCIONES PRINCIPALES PARA EL CALCULO DE RANKING
# =============================================================================

def calcularRanking(M, p):
    npages = M.shape[0]
    rnk = np.arange(0, npages) # ind{k] = i, la pagina k tienen el iesimo orden en la lista.
    scr = np.zeros(npages) # scr[k] = alpha, la pagina k tiene un score de alpha 
    # COIDGO
    A = calculo_A(M,p)
    n = A.shape[1]
    e = np.ones(n)
    L, U, cant_op = elim_gaussiana(A)
    y = scipy.linalg.solve_triangular(L,e, lower=True)
    x = scipy.linalg.solve_triangular(U,y)
    norma = scipy.linalg.norm(x,1)
    scr = x/norma
    rnk = ranking(scr)
    #
    return rnk, scr

def obtenerMaximoRankingScore(M, p):
    output = -np.inf
    # calculo el ranking y los scores
    rnk, scr = calcularRanking(M, p)
    output = np.max(scr)
    
    return output

# =============================================================================
# FUNCIONES PRINCIPALES PARA ANALISIS CUALITATIVO
# =============================================================================


def rankings_segunP(M):
    P = []  # Guardo los distintos valores de p
    mejores_paginas = []  # Guardo las páginas mejor rankeadas
    p = 0.95
    
    while p > 0:
        P.append(p)
        mejor_pagina = []  # Guardo las páginas mejor rankeadas con un p por si hay empate
        ranking, scores = calcularRanking(M, p)
        
        for j in range(len(scores)):
            if ranking[j]==1:
                mejor_pagina.append(j)
        p -= 0.05
        mejor_pagina = [x + 1 for x in mejor_pagina]  # Convertir a índice 1-based
        mejores_paginas.append(mejor_pagina)
    
    return P, mejores_paginas


def Graf_scores(W, nombre:str, p):
    sp= str(p)
    rnk, scores= calcularRanking(W, p)
    pagina= list(range(1,len(scores)+1))
    plt.figure(figsize=(10, 6))
    plt.scatter(pagina, scores, s=100,color = 'hotpink')
    plt.xlabel('paginas')
    plt.ylabel('puntaje')
    plt.xticks(pagina) 
    plt.title('puntaje de cada pagina del test ' + nombre + ' utilizando p =' + sp )
    plt.grid(True)
    plt.show()  
    

def Graf_MejoresPaginas_segunP (M, test): 
    d = M.shape[0]
    paginas = list(range(d))
    # Datos
    Altura = []
    P, mejores_paginas = rankings_segunP(M)
    for pagina in paginas:
        contador = 0
        for m in mejores_paginas:
            if pagina in m:
                contador += 1
        Altura.append(contador)
    print(Altura)
    # Crear el gráfico de barras horizontales
    plt.barh(paginas, Altura, color='darkseagreen', edgecolor='black')
    # Agregar etiquetas y título
    plt.xlabel('Cantidad de veces que obtuvieron el mejor Score')
    plt.ylabel('Páginas')
    plt.title('Cantidad de veces que una página obtuvo el mejor score variando el p en el test '+ test)
    # Establecer límites de los ejes
    plt.xlim(0, max(Altura) + 1)  
    plt.ylim(-0.5, d - 0.5) 
    # Mostrar el gráfico
    plt.show()
    
def graf_rankingP2(M, test): #matriz de conexiones, nombre del test 
    # Graficar quién fue la página mejor rankeada para cada p
    P, mejores_paginas = rankings_segunP(M)
    ylim= M.shape[0]
    # Inicializar listas para almacenar datos porsi hay empates 
    y1 = []
    x = []
    # Recorrer los valores de p y las páginas mejor rankeadas
    for i in range(len(mejores_paginas)):
        paginas = mejores_paginas[i]
        for pagina in paginas:
            y1.append(pagina)
            x.append(P[i]) 
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y1, s=100, color = "darkseagreen")
    plt.xlabel('P')
    plt.ylabel('Página Mejor Rankeada')
    plt.xticks(P)
    plt.yticks(range(0, ylim, 1)) 
    plt.title('Página Mejor Rankeada según el P utilizado  el test ' + test)
    plt.grid(True)
    plt.show()  
          
    #grafico que me muestre el porcentaje de veces que una pagina fue la mejor rankeada variando el P
def ranking_P3(M,test):
    P, mejores_paginas = rankings_segunP(M)
    paginas = list(range(0, M.shape[0])) 
    porcentajes=[]
    p=[]
    for pagina in paginas:
        contador = 0
        for m in mejores_paginas:
            if pagina in m: 
                contador += 1
        porcentaje = contador / len(mejores_paginas) * 100
        if porcentaje != 0:  # Excluir porcentajes de cero
            porcentajes.append(porcentaje)
            p.append(pagina)
    # Crear el gráfico de torta
    plt.pie(porcentajes, labels=p, autopct='%1.1f%%', wedgeprops={'linewidth': 2})
    plt.title('Porcentaje de veces que una página fue mejor rankeada con el test ' + test) 
    plt.show()


def variacionScore_p(M):
    rnk, scr= calcularRanking(M, 0.5)
    primeras3=[]
    j=0
    for i in range (0, len(rnk)):
        if rnk[i]== 1:
            primeras3.append(i)
    while len(primeras3)<3:
        if rnk[j]== 2 or rnk[j]==3:
            primeras3.append(j)
        j+=1
    a= primeras3[0]
    b= primeras3[1]
    c= primeras3[2]    
 
    P = []  # Guardo los distintos valores de p
    scr_a=[]
    scr_b=[]
    scr_c=[]
    p = 0.95
    
    while p > 0:
        P.append(p)
        ranking, scores = calcularRanking(M, p)
        scr_a.append(scores[a])
        scr_b.append(scores[b])
        scr_c.append(scores[c])
        p -= 0.05
    a+=1
    b+=1
    c+=1
    return P, a, scr_a, b, scr_b, c, scr_c

def graf_variacionScore(M, nombre):
    X,a, scr_a, b, scr_b, c, scr_c= variacionScore_p(M)
    # Gráfico de líneas
    fig, ax = plt.subplots(figsize=(12,8))
    fig.size= (8,8)
    ax.plot(X, scr_a, marker = "o", color= 'hotpink', label = "pagina"+ str(a))
    ax.plot(X, scr_b, marker = "o",color = 'darkseagreen', label = "pagina"+ str(b))
    ax.plot(X, scr_c, marker = "o", color= 'indigo', label = "pagina"+ str(c))
    ax.legend()
    plt.xlabel('P')
    plt.ylabel('puntaje')
    plt.xticks(X) 
    plt.title('puntaje las  3 paginas mejores rankeadas en el test ' + nombre + ' variando el P' )
    plt.grid(True)
    plt.show()  
# =============================================================================
# FUNCIONES PRINCIPALES PARA ANALISIS CUANTITATIVO
# =============================================================================
def tiempo_de_ejecucion(f, W, p):
    inicio = time.time()
    f(W,p)
    fin = time.time()
    tiempo_transcurrido = fin - inicio
    return tiempo_transcurrido

def tiempo_de_ejecucion_tamaño (n,p):
    i=2
    tiempos=[]
    tamaños=[]
    while i <=n:
        W= np.random.choice([0, 1], size=(i,i))
        np.fill_diagonal(W, 0)
        tiempo= tiempo_de_ejecucion(calcularRanking, W, p)
        tiempos.append(tiempo)
        tamaños.append(i)
        i+=1
    return tamaños, tiempos

def tiempo_ejecucion_tamaño_2 (n,p):
    i=2
    tiempos=[]
    tamaños=[]
    while i <=n:
        W= np.zeros((i, i))
        W[1][0]=1
        W[0][1]=1
        tiempo= tiempo_de_ejecucion(calcularRanking, W, p)
        tiempos.append(tiempo)
        tamaños.append(i)
        i+=1
    return tamaños, tiempos


def tiempo_ejecucion_densidad(n,p):
    W=  np.zeros((n, n))
    tiempo= tiempo_de_ejecucion(calcularRanking, W, p)
    tiempos= []
    nodos=[]
    conexiones = 0
    tiempos.append(tiempo)
    nodos.append(conexiones)
    for i in range (0,n):
        for j in range (0,n):
            if i!=j:
                W[i][j]=1
                tiempo= tiempo_de_ejecucion(calcularRanking, W, p)
                tiempos.append(tiempo)
                conexiones+=1
                nodos.append(conexiones)
    return nodos, tiempos


def graf_tiempo_tamaño():
    tamaño1, tiempo1= tiempo_de_ejecucion_tamaño(100, 0.5)
    tamaño2, tiempo2= tiempo_de_ejecucion_tamaño(100, 0.25)
# Crear el gráfico de dispersión con múltiples conjuntos de datos
    plt.scatter(tamaño1,tiempo1, color='darkseagreen', label='p=0.5')
    plt.scatter(tamaño2, tiempo2, color='darkseagreen', label='p=0.25')
    plt.plot(tamaño1,tiempo1, color='darkgreen', linestyle='-')
    plt.plot(tamaño2, tiempo2, color='darkseagreen', linestyle='-')
# Añadir etiquetas y leyenda
    plt.xlabel('dimensiones del grafo ')
    plt.ylabel('tiempo de ejecucion tardado [s]')
    plt.title('Tiempo de ejecucion del calculo del rankingpage segun el tamaño del grafo')
    plt.legend()
   # Mostrar el gráfico
    plt.grid(True)
    plt.show() 
    
def graf_tiempo_densidad(n):
    nodos,tiempo= tiempo_ejecucion_densidad(n,0.5)
# Crear el gráfico de dispersión con múltiples conjuntos de datos
    size= str(n)
    plt.scatter(nodos,tiempo, color='darkseagreen', label=size+"*"+size+" ,p=0.5")
# Añadir etiquetas y leyenda
    plt.xlabel('conexiones dentro del grafo ')
    plt.ylabel('tiempo de ejecucion tardado [s]')
    plt.title('Tiempo de ejecucion del calculo del rankingpage segun las conexiones entre paginas')
    plt.legend()
   # Mostrar el gráfico
    plt.grid(True)
    plt.show()    



# =============================================================================
# FUNCIONES PRINCIPALES PARA ANALISIS TEST DOS ESTRELLAS
# =============================================================================      

def test_dos_estrellas (p):
    W= leer_archivo(carpeta + "test_dosestrellas.txt")
    contador=0
    paginas= [1,2,3,4,8,9,10,11]
    resto= [5,6,7]
    for j in paginas:
        W[0][j]=1
        contador+=1
        ranking, scores= calcularRanking(W,p)
        maximoScore= np.max(scores)
        if scores[0]== maximoScore:
            print("se agregaron "+ str(contador)+" conexiones")
            print(W)
            return contador, W
        else:
            continue

def test_ds_2 (p):
    W= leer_archivo(carpeta + "test_dosestrellas.txt")
    contador=1
    l = W.shape[0]
    W[0][6]=1
    for j in range (0,l):
        if j!=5 and j!=6 and j!=8:
            W[0][j]=1
            contador+=1
            ranking, scores= calcularRanking(W,p)
            maximoScore= np.max(scores)
            if scores[0]== maximoScore:
                print("se agregaron "+ str(contador)+" conexiones")
                print(W)
                return contador, W
            else:
                continue

def comparacion_DS(p):
    W= leer_archivo(carpeta + "test_dosestrellas.txt")
    contador, W2= test_dos_estrellas(p)
    rnk, scoresW = calcularRanking (W, p)
    rank, scoresW2= calcularRanking (W2, p)
    paginas = [1,2,3,4,5,6,7,8,9,10,11,12]
    plt.figure(figsize=(10, 6))
    plt.scatter(paginas, scoresW, s=100, color= 'hotpink')
    plt.scatter(paginas, scoresW2, s=100,  color='darkseagreen')
    plt.xlabel('Paginas del test dos estrellas')
    plt.ylabel('scores de las paguinas del test dos estrellas')
    plt.xticks(paginas)
    plt.title('Comparacion de los puntajes obtenidos con el test Dos estrellas original y agregando conexiones a la pagina 1')
    plt.grid(True)
    plt.show()

def comparacion_DS2(p):
    W= leer_archivo(carpeta + "test_dosestrellas.txt")
    contador, W2= test_ds_2 (p)
    rnk, scoresW = calcularRanking (W, p)
    rank, scoresW2= calcularRanking (W2, p)
    paginas = [1,2,3,4,5,6,7,8,9,10,11,12]
    plt.figure(figsize=(10, 6))
    plt.scatter(paginas, scoresW, s=100, color= 'hotpink')
    plt.scatter(paginas, scoresW2, s=100,  color='darkseagreen')
    plt.xlabel('Paginas del test dos estrellas')
    plt.ylabel('scores de las paguinas del test dos estrellas')
    plt.xticks(paginas)
    plt.title('Comparacion de los puntajes obtenidos con el test Dos estrellas original y agregando conexiones a la pagina 1')
    plt.grid(True)
    plt.show()
    
    
# =============================================================================
# FUNCIONES PRINCIPALES PARA ANALISIS TEST PROPIOS
# =============================================================================      
#ninguno conectado
ninguno_conectado= np.zeros((20, 20))
rankingNingunoConectados, scoresNingunoConectados= calcularRanking (ninguno_conectado, 0.5)
 

#todos conectados 1
#la pagina i solo es linkeada y linkea a la pagina j

#todos conectados 

def todos_conectados ():
    W = np.zeros((20, 20))
    for j in range (0, 20):
         for i in range(0,20):
             if i !=j:
                 W[i][j]=1
    return W
       
#todos conectados 2
def matriz_t2():
    W = np.zeros((20, 20))
    for j in range (0, 20):
        if j!=5 and j!=16 and j!=10: 
            W[5][j]= 1
        if j!= 7 and j!=16 and j!=10:
            W[7][j]=1
    W[15][5]=1
    W[15][7]=1
    W[15][16]=1
    W[15][10]=1
    return W

def experimento_cond(W):
## Calculo los valores de la condicion
    cond = []
    ps=[]
    p= 0.95
    while p >0:
        matriz = calculo_A(W, p)
        res = np.linalg.cond(matriz)
        cond.append(res)
        ps.append(p)
        p-=0.05
## Armo el grafico
    plt.plot(ps, cond, color='darkseagreen', marker = 'o') ##(valores de x, valores de y)
    
    # Etiquetar los ejes
    plt.xlabel('Valor de p')
    plt.ylabel('Cond(A)')
    plt.title('Numero de condicion en funcion del valor p')
    plt.grid(True)
    
    # Mostrar el gráfico
    plt.show()


