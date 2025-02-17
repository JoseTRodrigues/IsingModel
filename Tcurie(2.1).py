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
from numpy import array,exp,linspace,savetxt
from scipy.constants import k
from scipy.ndimage import convolve, generate_binary_structure
from matplotlib.pyplot import plot,title,axis,xlabel,ylabel,yscale,\
    imshow,colorbar,show



#%%################### ENERGIA DO SISTEMA e VIZINHANÇA ########################

#desconsiderando campo externo
#Cálculo de H/J, J é desconhecido

def H(lattice):
    # applies the nearest neighbours summation
    kern = generate_binary_structure(2, 1) 
    kern[1][1] = False
    arr = -lattice * convolve(lattice, kern, mode='constant', cval=0)
    return arr.sum()

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
    spin=[-1,1] #valores possíveis de spin
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

def mup(n):
    up=[]
    for i in range(n):
        up.append([])
        for j in range(n):
            up[i].append(1) #ou -1 (down)
    return up


#%%######################## MODELO DE ISING    ################################

def ising(mspin1,B,steps):
    
    test=deepcopy(mspin1) #matriz de teste
    
    Hmin=H(array(mup(n)))
    
    E=[H(array(mspin1))] #energia da matriz inicial

    for step in range(steps):
        
        'ALTERAÇÃO ALEATÓRIA DE 1 SPIN DE COORD (x,y)'
        x=randint(0,len(mspin1)-1)
        y=randint(0,len(mspin1)-1)
        
        viz=von_neumann(test, (x,y)) #VIZINHANÇA DO SPIN de coord (x,y)
        
        test[x][y]=-1*mspin1[x][y] #matriz com 1 spin alterado
        
        'ENERGIA DO DOMÍNIO DO SPIN DE COORD (x,y)'
        Hf=sum(test[x][y]*array(viz))
        Hi=sum(mspin1[x][y]*array(viz))
        dH=Hf-Hi
        
        if dH > 0 and random() < exp(-dH*B):    
            mspin1[x][y]=-1*mspin1[x][y]
            
        elif dH <= 0:
            mspin1[x][y]=-1*mspin1[x][y]
        
        # elif H(array(mspin1))==Hmin:
        #     break
            
        E.append(H(array(mspin1)))
    
    # imshow(test) # Matriz final
    # imshow(mspin1)
    # title(f'\u03B2 = {B}')
    # colorbar(imshow(mspin1),shrink=0.3)
    # axis('off')
    # show()
    
    return E,m



#%%####################### TEMPERATURA DE CURIE ###############################

if __name__  == "__main__":
    
    B=linspace(0.1,5,30) #Beta <=> Temperatura
    steps=1000000
    
    n=20 #dimensão da matriz
    mspin1=mspin(n) #matriz inicial
    imshow(mspin1)
    title('Distribuição Inicial')
    colorbar(imshow(mspin1),shrink=0.3)
    axis('off')
    show()
    
    E_med=[]
    for i in B:
        E = ising(mspin1 , i , steps)
        E_med.append(array(E[-500000:]).mean())
        # print(i*30)
    
    T=1/B
    plot(T,E_med)
    # yscale('log')
    title('Evolução da Energia com a Temperatura')
    xlabel(r'$(\frac{k}{J})$ T')
    ylabel('E')
    show()

    T_curie=T[E_med.index(max(E_med))]/k
    print('\nTemperatura de Curie = ',T_curie,'[K/J]')

    """
    EXPORT DADOS FINAIS
    """
    savetxt('Energia.csv', E_med, delimiter=',', fmt='%s')
    savetxt("Temperatura.csv", 1/B ,delimiter =", ", fmt ='% s')
