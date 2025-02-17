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
from matplotlib.pyplot import plot,xlabel,ylabel,subplots,title,axis,imshow,colorbar,show



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

#%%#######################  MATRIZES DE SPINS     #############################

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


#%%#######################   MODELO DE ISING      #############################

def ising(mspin1,B,steps):
    
    test=deepcopy(mspin1)   #matriz de teste
    mspin2=deepcopy(mspin1) #matriz final
    
    # Hmin=H(array(mup(n)))
    
    E=[H(array(mspin1))] #energia da matriz inicial
    
    m=[] #magnetização
    
    for step in range(steps):
        
        'ALTERAÇÃO ALEATÓRIA DE 1 SPIN DE COORD (x,y)'
        x=randint(0,len(mspin2)-1)
        y=randint(0,len(mspin2)-1)

        viz=von_neumann(test, (x,y)) #VIZINHANÇA DO SPIN de coord (x,y)
        
        test[x][y]=-1*mspin2[x][y] #matriz com 1 spin alterado
        
        'ENERGIA DO DOMÍNIO DO SPIN DE COORD (x,y)'
        Hf=sum(test[x][y]*array(viz))
        Hi=sum(mspin2[x][y]*array(viz))
        dH=Hf-Hi
        
        if dH <= 0:
            mspin2[x][y]=-1*mspin2[x][y]
        
        elif dH > 0 and random() < exp(-dH*B):    
            mspin2[x][y]=-1*mspin2[x][y]
            
        # elif H(array(mspin2))==Hmin:
        #     break
            
        E.append(H(array(mspin2)))
        # m.append(array(mspin2).sum()/n**2)

    # imshow(mspin2)
    # title('\u03B2 = ' '{:.2f}'.format(float(f'{B}')))
    # colorbar(imshow(mspin2),shrink=0.3)
    # axis('off')
    # show()

    return E
    # return E,m
    # return


#%%####################### TEMPERATURA DE CURIE   #############################

if __name__  == "__main__":
    
    steps=1000000
    pontos=30
    n=50 #dimensão da matriz
    
    T=linspace(0.02,10,pontos) # Temperatura
    B=1/T
    
    mspin1=mspin(n) #matriz inicial
    imshow(mspin1)
    title('Distribuição Inicial de Spins')
    colorbar(imshow(mspin1),shrink=0.3)
    axis('off')
    show()
    
    E_med=[]
    m_med=[]
    for i in B:
        E,m = ising(mspin1 , i , steps)
        E_med.append(array(E[-int(steps/2):]).mean())
        m_med.append(array(m[-int(steps/2):]).mean())
        print(1/i*pontos)
        
#%%#######################    PLOTS E DADOS       #############################

    fig, (ax1,ax2) = subplots(2,sharex=True)
    
    fig.suptitle('Energia e Magnetização')
    ax1.plot(T,abs(array(E_med)),'o--')
    ax1.set_ylabel(r'$\bar{E}$')
    ax1.set_ylim(-10000)
    
    # ax2.scatter(T,m_med,c='r')
    ax2.plot(T,abs(array(m_med)),'o--',c='r')
    ax2.set_xlabel(r'$(\frac{k}{J})$ T')
    ax2.set_ylabel(r'$\bar{m}$')
    ax2.set_xlim(-0.2,11)
    ax2.set_ylim(-1,1)
    show()

    """
    T de Curie
    """
    T_curie=(T[7]+T[6])/2 #obtido graficamente
    # T_curie=T[E_med.index(max(E_med))]
    # print('\nTemperatura de Curie = ',T_curie,'[K/J]')


    # """
    # EXPORT DADOS FINAIS
    # """
    # savetxt('Energia.csv', E_med, delimiter=',', fmt='%s')
    # savetxt('Magnetização.csv', m_med, delimiter=',', fmt='%s')
    # savetxt("Temperatura.csv", T ,delimiter =", ", fmt ='% s')


    # plot(E)
    # title('Energia para \u03B2 = 0.75 ')
    # xlabel('Steps')
    # ylabel('E')