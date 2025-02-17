# -*- coding: utf-8 -*-
"""
MUDANÇA DE FASE DE MATERIAIS FERROMAGNÉTICOS

@author: José Rodrigues, nº 2019246536

"""
#%%
from random import random,choice,randint
from copy import deepcopy
from numpy import array,exp
from matplotlib.pyplot import plot,scatter,hist,show




n=5 #nº de linhas/colunas da matriz de spins
spin=[-1,1] #valores possíveis de spin

#%%
"""
########################## ENERGIA DO SISTEMA #################################
"""
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
            E.append(matrix[i][j]*array(viz[i][j]))
    return sum(E)

# up=[]
# for i in range(n):
#     up.append([])
#     for j in range(n):
#         up[i].append(1) #ou -1 (down)
        
# Hmin=H(up)

#%%
"""
########################## MATRIZES DE SPINS ##################################
"""

def mspin():
    'DISTRIBUIÇÃO INICIAL DE SPIN'
    matrix=[]
    for i in range(n):
        matrix.append([])
        for j in range(n):
            matrix[i].append(choice(spin))
    return matrix     

def flip(matrix):
    'ALTERAÇÃO ALEATÓRIA DE 1 SPIN'
    matrix=deepcopy(matrix)
    i=randint(0,n-1)
    j=randint(0,n-1)
    matrix[i][j]=-1*matrix[i][j]

    return matrix     

def mprint(matrix):
    'print de matriz'
    for i in range(n):
        print(array(matrix)[i])
    return

#%%
"""
######################## MODELO DE ISING ######################################
"""

steps=5000

B=0.7 #Beta=1/KT

E=[]

mspin1=mspin()
test=flip(mspin1)

# print('\nDistribuição inicial de spins:\n')
# mprint(mspin1)

for i in range(steps):
    Hi=H(mspin1)
    Hf=H(test)
    dH=Hf-Hi
    
    if Hf > Hi and random() < exp(-(dH)*B):    
        E.append(Hf)
        mspin1=test
        test=flip(test)
    
    elif Hf <= Hi:
        E.append(Hf)
        mspin1=test
        test=flip(test)
        
    print(i)    

plot(E)
show()