import numpy as np
import networkx as nx
import scipy
import matplotlib.pyplot as plt
import time
# =============================================================================
# FUNCION PARA LEER EL ARCHIVO
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
    
    ## desde aqui -- CODIGO A COMPLETAR
    
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
    
    ## hasta aqui

      
    L = np.tril(Ac,-1) + np.eye(A.shape[0]) 
    U = np.triu(Ac)
    return L, U, cant_op
    

def calculo_k(fila_actual, divisor, iterador):
    multiplicador = 0
    if divisor != 0:
        multiplicador = fila_actual[iterador] / divisor
    return multiplicador
   
def ranking(score):
    sorted_score = sorted(scr)
    rnk = [sorted_score.index(x) + 1 for x in scr]
    res = []
    for elemento in rnk:
        if elemento in res:
            res.append(elemento + 1)
            
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
    rnk = [sorted_score.index(x) + 1 for x in scr]
    #
    return rnk, scr

def obtenerMaximoRankingScore(M, p):
    output = -np.inf
    # calculo el ranking y los scores
    rnk, scr = calcularRanking(M, p)
    output = np.max(scr)
    
    return output

#ARCHIVOS DE ENTRADA
archivo_test = './tests/test_dosestrellas.txt'
    
#CARGA DE ARCHIVO EN GRAFO
W = leer_archivo(archivo_test)

dibujarGrafo(W, print_ejes=False)
# defino la probabilidad de salto de continuar los links de la pagina actual
p = 0.5
# Realizo el test unitario para el calculo del mayor score, que pruebe que el codigo funciona correctamente.
print('*'*50)
print('Test unitario 1')
try:
    compara = obtenerMaximoRankingScore(W, p)
    assert(np.isclose(compara, 0.1811, atol= 0.0001))
except:
    print('OUCH!! - No paso el test unitario')
else:
    print('BIEN! - Paso correctamente el test unitario')
print('*'*50)

rnk, scr = calcularRanking(W, p)
print(rnk)
print(scr)

# =============================================================================
# FUNCIONES PRINCIPALES PARA ANALISIS CUALITATIVO
# =============================================================================

def rankings_segunP(M, n):
    P = [] #guardo los distintos P
    mejores_paginas = [] #guardo las paginas mejores ranqueadas
    for i in range(1, n+1):
        p = 1 / i
        P.append(p)
        mejor_pagina = [] #guardo las paginas mejores rankeadas con un P por si hay empate
        ranking, scores = calcularRanking(M, p)
        maxScore = np.max(scores)
        
        for j in range(len(scores)):
            score = scores[j]
            if score == maxScore:
                m = j 
                mejor_pagina.append(m)
                
        mejores_paginas.append(mejor_pagina)

    return P, mejores_paginas

    
def Graf_MejoresPaginas_segunP (M, n): 
    d = M.shape[0]
    paginas = list(range(d))
    # Datos
    Altura = []
    P, mejores_paginas = rankings_segunP(M, n)
    for pagina in paginas:
        contador = 0
        for m in mejores_paginas:
            if pagina in m:
                contador += 1
        Altura.append(contador)
    print(Altura)
    paleta = "Pastel1"
    # Crear el gráfico de barras horizontales
    plt.barh(paginas, Altura, color=plt.cm.get_cmap(paleta)(range(len(paginas))), edgecolor='black')
    # Agregar etiquetas y título
    plt.xlabel('Cantidad de veces que obtuvieron el mejor Score')
    plt.ylabel('Páginas')
    plt.title('Cantidad de veces que una página obtuvo el mejor score variando el p')
    # Establecer límites de los ejes
    plt.xlim(0, max(Altura) + 1)  
    plt.ylim(-0.5, d - 0.5) 
    # Mostrar el gráfico
    plt.show()
  

    
def graf_rankingP2(M, n):  
    # Graficar quién fue la página mejor rankeada para cada p
    P, mejores_paginas = rankings_segunP(M, n)
    
    # Inicializar listas para almacenar datos porsi hay empates 
    y1 = []
    x = []
    # Recorrer los valores de p y las páginas mejor rankeadas
    for m, paginas in enumerate(mejores_paginas):
        for pagina in paginas:
            y1.append(pagina)
            x.append(P[m])
    
    plt.scatter(x, y1, s=100, c=x, cmap='spring')
    
    plt.xlabel('P')
    plt.ylabel('Página Mejor Rankeada')
    plt.title('Página Mejor Rankeada según el P utilizado')
    plt.colorbar(label='X')
    plt.show()  
              
    #grafico que me muestre el porcentaje de veces que una pagina fue la mejor rankeada variando el P
def ranking_P3(M, n):   
    P, mejores_paginas = rankings_segunP(M, n)
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
    paleta = "Pastel1"
    # Crear el gráfico de torta
    plt.pie(porcentajes, labels=p, colors=plt.cm.get_cmap(paleta)(range(len(porcentajes))), autopct='%1.1f%%', wedgeprops={'linewidth': 2})
    plt.title('Porcentaje de veces que una página fue mejor rankeada') 
    plt.show()

def medir_tiempo(f, M, n):
    inicio = time.time()
    resultado = f (M,n)
    fin = time.time()
    tiempo_transcurrido = fin - inicio
    print(f"Tiempo de ejecución de {f.__name__}: {tiempo_transcurrido} segundos")
    return resultado

# =============================================================================
# FUNCIONES PRINCIPALES PARA ANALISIS CUANTITATIVO
# =============================================================================
def tiempo_de_ejecucion(f, W, p):
    inicio = time.time()
    resultado = f(W,p)
    fin = time.time()
    tiempo_transcurrido = fin - inicio
    return tiempo_transcurrido

def tiempo_de_ejecucion_tamaño (n,p):
    i=2
    tiempos=[]
    tamaños=[]
    while i <=n:
        W= np.random.choice([0, 1], size=(i,i))
        tiempo= tiempo_de_ejecucion(calcularRanking, W, p)
        tiempos.append(tiempo)
        tamaños.append(i)
        i+=1
    return tiempos, tamaños

def tiempo_ejecucion_densidad (n,p):
    W=  np.zeros((n, n))
    tiempo= tiempo_de_ejecucion(calcularRanking, W, p)
    tiempos= []
    tiempos.append(tiempo)
    for i in range (0,n):
        for j in range (0,n):
            W[i][j]=1
            tiempo= tiempo_de_ejecucion(calcularRanking, W, p)
            tiempos.append(tiempo)
    return tiempos