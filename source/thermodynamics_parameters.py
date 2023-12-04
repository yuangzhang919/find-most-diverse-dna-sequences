#!/usr/bin/env python3

#### input the parameter lookup table from the unified, SL, 1998 paper

import pandas as pd
import numpy as np

R = 1.987 / 1000 # ideal gas constant, kcal/(K mol)
beta = 1. # factor in the partition function for binding initialization.
omega0 = 1.26e-4 # Cooperativity, sigma,in the power law of exact loop entropy. Value from Blossey and Carlon 2003.
alpha = 2.15 # power law exponent of exact loop entropy. Value from Blossey and Carlon 2003.


seqs = ["AA", "TT", 
        "AT", 
        "TA", 
        "CA", "TG",
        "GT", "AC",
        "CT", "AG",
        "GA", "TC",
        "CG", 
        "GC",
        "GG", "CC",
        "term_GC",
        "term_AT"]

delta_G_37s = [-1.00, -1.00,
               -0.88, 
               -0.58,
               -1.45, -1.45,
               -1.44, -1.44,
               -1.28, -1.28,
               -1.30, -1.30,
               -2.17, 
               -2.24,
               -1.84, -1.84,
               0.98,
               1.03] 

delta_H = [-7.9, -7.9,
           -7.2,
           -7.2,
           -8.5, -8.5,
           -8.4, -8.4, 
           -7.8, -7.8,
           -8.2, -8.2,
           -10.6, 
           -9.8,
           -8.0, -8.0,
           0.1,
           2.3]

delta_S = [-22.2, -22.2,
           -20.4,
           -21.3,
           -22.7, -22.7,
           -22.4, -22.4,
           -21.0, -21.0,
           -22.2, -22.2,
           -27.2, 
           -24.4,
           -19.9, -19.9, 
           -2.8,
           4.1]

dimer_parameter_lookup_list = []

for i,seq in enumerate(seqs):
    dimer_parameter_lookup_dict = {}
    dimer_parameter_lookup_dict["seq"] = seq
    dimer_parameter_lookup_dict["delta_G_37"] = delta_G_37s[i]
    dimer_parameter_lookup_dict["delta_H"] = delta_H[i]
    dimer_parameter_lookup_dict["delta_S"]  = delta_S[i]/1000
    
    dimer_parameter_lookup_list.append(dimer_parameter_lookup_dict)
    
df_parameters = pd.DataFrame(dimer_parameter_lookup_list)