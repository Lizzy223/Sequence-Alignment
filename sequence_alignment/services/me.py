
# bpy24.py


from __future__ import division, print_function
import numpy as np
import pandas as pd
from numpy import full

from IPython.core.display import display, HTML
display(HTML('<h1>Hello, world!</h1>'))

def pretty_table_from_array(data_array, row_labels,col_labels):
    """Show an HTML table from a 2d numpy array"""
    df = pd.DataFrame(data_array,index=row_labels,columns=col_labels)
    return df
row_labels = []

pt ={'match': 1, 'mismatch': -1, 'gap': -1}

def mch(alpha, beta):
    if alpha == beta:
        return pt['match']
    elif alpha == '-' or beta == '-':
        return pt['gap']
    else:
        return pt['mismatch']

def water(s1, s2):
    m, n = len(s1), len(s2)
    H = np.zeros((m+1, n+1))    
    T = np.zeros((m+1, n+1))
    max_score = 0
    scoring_array = full([n,m],0)
    # Score, Pointer Matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            sc_diag = H[i-1][j-1] + mch(s1[i-1], s2[j-1])
            sc_up = H[i][j-1] + pt['gap']
            sc_left = H[i-1][j] + pt['gap']
            H[i][j] = max(0,sc_left, sc_up, sc_diag)
            if H[i][j] == 0: T[i][j] = 0
            if H[i][j] == sc_left: T[i][j] = 1
            if H[i][j] == sc_up: T[i][j] = 2
            if H[i][j] == sc_diag: T[i][j] = 3
            if H[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = H[i][j];
    
    print('H=\n',H,'\n')
    # display(water(scoring_array,s1,s2))
    # print('T=\n',T,'\n')
    align1, align2 = '', ''
    i,j = max_i,max_j
    
   
    #Traceback
    while T[i][j] != 0:
        if T[i][j] == 3:
            a1 = s1[i-1]
            a2 = s2[j-1]
            i -= 1
            j -= 1
        elif T[i][j] == 2:
            a1 = '-'
            a2 = s2[j-1]
            j -= 1
        elif T[i][j] == 1:
            a1 = s1[i-1]
            a2 = '-'
            i -= 1
        print('%s ---> a1 = %s\t a2 = %s\n' % ('Add',a1,a2))
        align1 += a1
        align2 += a2

    align1 = align1[::-1]
    align2 = align2[::-1]
    sym = ''
    iden = 0
    for i in range(len(align1)):
        a1 = align1[i]
        a2 = align2[i]
        if a1 == a2:                
            sym += a1
            iden += 1
        elif a1 != a2 and a1 != '-' and a2 != '-': 
            sym += ' '
        elif a1 == '-' or a2 == '-':          
            sym += ' ' 3sl
    
    identity = iden / len(align1) * 100
    print('Identity = %f percent' % identity)
    print('Score =', max_score)
    print(align1)
    # print(sym)
    # print(align2)
# display('H=\n','\n','T=\n','\n')
# def table_array(np_array, H, T):
#         df = pd.DataFrame(np_array, index=H, columns=T)
#         return df

# display(table_array(H,T))




if __name__ == '__main__':
    water('CTCCCGCT','TATCC')