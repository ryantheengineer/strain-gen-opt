# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:21:06 2024

@author: Ryan.Larson
"""

import numpy as np
import plotly.express as px

def principal_strains(strainxx, strainyy, strainshear):
    avg_strain = (strainxx + strainyy)/2
    sqrt_term = np.sqrt(((strainxx - strainyy)/2)**2 + strainshear**2)
    
    principal_strain1 = avg_strain + sqrt_term
    principal_strain2 = avg_strain - sqrt_term
    return principal_strain1, principal_strain2

strainxx = np.linspace(0,10000,11,endpoint=True)
strainyy = np.linspace(0,10000,11,endpoint=True)
strainshear = np.linspace(0,10000,11,endpoint=True)

principal_strain1 = np.zeros((len(strainxx), len(strainyy), len(strainshear)))
principal_strain2 = np.zeros((len(strainxx), len(strainyy), len(strainshear)))

for i in range(len(strainxx)):
    for j in range(len(strainyy)):
        for k in range(len(strainshear)):
            principal_strain1[i,j,k], principal_strain2[i,j,k] = principal_strains(strainxx[i], strainyy[j], strainshear[k])
            
hue_values1 = principal_strain1.flatten()
hue_values2 = principal_strain2.flatten()

data1 = {'strainxx': np.repeat(strainxx, len(strainyy)*len(strainshear)),
        'strainyy': np.tile(np.repeat(strainyy, len(strainshear)), len(strainxx)),
        'strainshear': np.tile(strainshear, len(strainxx)*len(strainyy)),
        'hue': hue_values1}
data2 = {'strainxx': np.repeat(strainxx, len(strainyy)*len(strainshear)),
        'strainyy': np.tile(np.repeat(strainyy, len(strainshear)), len(strainxx)),
        'strainshear': np.tile(strainshear, len(strainxx)*len(strainyy)),
        'hue': hue_values2}

fig1 = px.scatter_3d(data1, x='strainxx', y='strainyy', z='strainshear', color='hue', opacity=0.8)
# fig2 = px.scatter_3d(data2, x='strainxx', y='strainyy', z='strainshear', color='hue', opacity=0.8)

fig1.show(renderer='browser')
            

# ncombos = 1000
# strain1_greater = []
# strain2_greater = []
# for i in range(ncombos):
#     xxindex = np.random.choice(strainxx.shape[0],1)
#     yyindex = np.random.choice(strainyy.shape[0],1)
#     shearindex = np.random.choice(strainshear.shape[0],1)
    
#     xx = strainxx[xxindex]
#     yy = strainyy[yyindex]
#     shear = strainshear[shearindex]
    
#     strain1, strain2 = principal_strains(xx, yy, shear)
    
#     if all(val > strain1 for val in [xx, yy, shear]):
#         strain1_greater.append(True)
#     else:
#         strain1_greater.append(False)
        
#     if all(val > strain2 for val in [xx, yy, shear]):
#         strain2_greater.append(True)
#     else:
#         strain2_greater.append(False)

# pct_strain1 = strain1_greater.count(True)/len(strain1_greater)
# pct_strain2 = strain2_greater.count(True)/len(strain2_greater)

# print(f'{pct_strain1} of Principal Strain 1 values were greater than all inputs')
# print(f'{pct_strain2} of Principal Strain 2 values were greater than all inputs')
