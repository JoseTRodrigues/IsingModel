# -*- coding: utf-8 -*-
"""
MUDANÇA DE FASE DE MATERIAIS FERROMAGNÉTICOS

@author: José Rodrigues, nº 2019246536

Fontes: https://codereview.stackexchange.com/questions/199754/codewars-n-dimensional-von-neumann-neighborhood-in-a-matrix
        https://github.com/lukepolson/youtube_channel/blob/main/Python%20Metaphysics%20Series/vid14.ipynb
"""
#%% IMPORTS
from random import random,choice,randint
from copy import deepcopy
from numpy import array,exp
from matplotlib.pyplot import plot,imshow,show

#%%######################## ENERGIA DO SISTEMA ################################

#desconsiderando campo externo
#Cálculo de H/J, J é desconhecido

def H(matrix):
    viz=[]
    'INETRAÇÕES ENTRE DOMÍNIOS DE SPIN'
    for i in range(n):
        viz.append([])
        for j in range(n):
            if i==n-1 and j==n-1:
                viz[i].append([])
            elif j==n-1:
                viz[i].append([matrix[i+1][j]])
            elif i==n-1:
                viz[i].append([matrix[i][j+1]])  
            elif i!=n-1 and j!=n-1:
                viz[i].append([matrix[i][j+1],matrix[i+1][j]])
                
    E=[]           
    'HAMILTONIANO da INTERAÇÃO INTERNA'
    for i in range(n):
        for j in range(n):
            E.append(-sum(matrix[i][j]*array(viz[i][j])))
    return sum(E)

def von_neumann(matrix, coordinates, distance=1): #coordenadas (m,n)
    'VIZINHANÇA DE VON NEUMANN'
    dimensions = len(coordinates)
    neigh = []
    app = neigh.append

    def recc_von_neumann(arr, curr_dim=0, remaining_distance=distance, isCenter=True):
        #the breaking statement of the recursion
        if curr_dim == dimensions:
            if not isCenter:
                app(arr)
            return

        dimensions_coordinate = coordinates[curr_dim]
        if not (0 <= dimensions_coordinate < len(arr)):
            return 

        dimesion_span = range(dimensions_coordinate - remaining_distance, 
                              dimensions_coordinate + remaining_distance + 1)
        for c in dimesion_span:
            if 0 <= c < len(arr):
                recc_von_neumann(arr[c], 
                                  curr_dim + 1, 
                                  remaining_distance - abs(dimensions_coordinate - c), 
                                  isCenter and dimensions_coordinate == c)
        return
    
    recc_von_neumann(matrix)
    return neigh


#%%######################## MATRIZES DE SPINS  ################################

def mspin(n):
    'DISTRIBUIÇÃO INICIAL DE SPIN'
    matrix=[]
    for i in range(n):
        matrix.append([])
        for j in range(n):
            matrix[i].append(choice(spin))
    return matrix     

def mprint(matrix):
    'print de matriz'
    for i in range(len(matrix)):
        print(array(matrix)[i])
    return

def up(n):
    up=[]
    for i in range(n):
        up.append([])
        for j in range(n):
            up[i].append(1) #ou -1 (down)
    return up


#%%######################## MODELO DE ISING    ################################

n=20 #dimensão da matriz
spin=[-1,1] #valores possíveis de spin

mspin1=mspin(n) #matriz inicial
imshow(mspin1)
show()
# B=0.7
# steps=10000

def ising(mspin1,B,steps):
    
    E1=H(mspin1) #energia da matriz inicial
    
    test=deepcopy(mspin1) #matriz de teste
    
    E=[]
    dE=[]
    
    for s in range(steps):
        
        'ALTERAÇÃO ALEATÓRIA DE 1 SPIN DE COORD (x,y)'
        x=randint(0,len(mspin1)-1)
        y=randint(0,len(mspin1)-1)
        
        viz=von_neumann(test, (x,y)) #VIZINHANÇA DO PONTO (x,y)
        
        test[x][y]=-1*mspin1[x][y] #matriz com 1 spin alterado
        
        'ENERGIA DO DOMÍNIO DO SPIN DE COORD (x,y)'
        Hf=sum(test[x][y]*array(viz))
        Hi=sum(mspin1[x][y]*array(viz))
        dH=Hf-Hi
        
        if dH > 0 and random() < exp(-dH*B):    
            mspin1[x][y]=-1*mspin1[x][y]
            # E.append(H(mspin1))
            E1+=dH
            E.append(E1)
            dE.append(dH)
            
        elif dH <= 0:
            mspin1[x][y]=-1*mspin1[x][y]
            # E.append(H(mspin1))
            E1+=dH
            E.append(E1)
            dE.append(dH)
    
    return E, dE

#%%

E,dE = ising(mspin1, 0.7, 10000)
plot(E)
show()

# plot(dH)
# show()
