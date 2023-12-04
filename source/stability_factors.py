#!/usr/bin/env python3

from thermodynamics_parameters import df_parameters as df_parameters
from thermodynamics_parameters import R, beta, omega0, alpha

import numpy as np



def stability_factor_lookup(T, 
                            c_mono = 0.02, 
                            model = "SL_oligo", 
                            N = 100):
    """
    Function to buildup the stability factor lookup table for each oligo dimer, termial base, and 010 isolated base.
    
    Input:
    ========================================================================
    T: float - temperature in Kelvin, the temperature points to calculate the stability factors.
    model: Thermodynamic parameter sets to describe stacking energy. "SL_oligo" for the DNA oligo model in SantaLucia PNAS 1998, or "SL_polymer" for the DNA polymer model in SantaLucia PNAS 1998. Default is "SL_oligo"
    c_mono: Salt concentration, salt correction is given by SantaLucia PNAS 1998 
    
    Output:
    ========================================================================
    s010_lookup, send_lookup, s11_lookup: dictionaries - Lookup tables.
    """
    
    s010_lookup = {}
    send_lookup = {}
    s11_lookup = {}
    
    for nt in "AGCT":
        term = "term_GC" if nt in "GC" else "term_AT"
        delta_H = df_parameters[df_parameters["seq"] == term]["delta_H"].iloc[0]
        delta_S = df_parameters[df_parameters["seq"] == term]["delta_S"].iloc[0]
        
        if model == "SL_polymer":
            delta_G = delta_H - T*delta_S - 0.175*np.log(c_mono) - 0.20
        
        elif model == "SL_oligo":
            delta_S = delta_S + 0.368*(N-1)/N*np.log(c_mono)/1000
            delta_G = delta_H - T*delta_S 
            
        s010_lookup[nt] = np.exp(-(2*(delta_G)/R/T))
        send_lookup[nt] = np.exp(-((delta_G)/R/T))
        
    for dimer in df_parameters["seq"].to_numpy()[:-2]:
        delta_H = df_parameters[df_parameters["seq"] == dimer]["delta_H"].iloc[0]
        delta_S = df_parameters[df_parameters["seq"] == dimer]["delta_S"].iloc[0]
        
        if model == "SL_polymer":
            delta_G = delta_H - T*delta_S - 0.175*np.log(c_mono) - 0.20
        
        elif model == "SL_oligo":
            delta_S = delta_S + 0.368*(N-1)/N*np.log(c_mono)/1000
            delta_G = delta_H - T*delta_S 
            
        s11_lookup[dimer] = np.exp(-((delta_G)/R/T))
    
    return s010_lookup, send_lookup, s11_lookup

