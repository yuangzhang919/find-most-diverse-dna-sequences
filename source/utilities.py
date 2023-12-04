#!/usr/bin/env python3
import numpy as np
from recursion import recursion
from thermodynamics_parameters import R, beta, omega0, alpha


def helicity(sequence, 
             Ts,
             model = "SL_oligo",
             c_mono = 0.02):
    """
    Function to calculate helicity values at a given range of Temperatures.
    
    Input:
    ========================================================================
    sequence: str - the sequence of DNA to calculate the melting curve.
    Ts: list or 1D np.array - temperature in Kelvin, the temperature points to calculate helicity values
    model: Thermodynamic parameter sets to describe stacking energy. "SL_oligo" for the DNA oligo model in SantaLucia PNAS 1998, or "SL_polymer" for the DNA polymer model in SantaLucia PNAS 1998. Default is "SL_oligo"
    c_mono: Salt concentration, salt correction is given by SantaLucia PNAS 1998 
    
    Output:
    ========================================================================
    helicities: 1D np.array - helicity values corresponding to the input temperature values
    """
    p_profiles_stacks_in_T =np.zeros([len(Ts),len(sequence)])
    
    helicities = []
    for i,T in enumerate(Ts):
        p_profiles_stacks_in_T[i] = p_closed_profile(sequence, 
                                                     T = T,
                                                     model = model,
                                                     c_mono = c_mono)

        if T//10 == 0:
            print(f"{T} done")
    
    helicities = np.average(p_profiles_stacks_in_T, axis = 1)
    
    
    
    return helicities
    
def p_closed_profile(sequence, 
                     T,
                     model = "SL_oligo",
                     c_mono = 0.02):
    
    """
    Function to calculate the probabiliy profile of close bases.
    
    Input:
    ========================================================================
    sequence: str - the sequence of DNA to calculate the probability profile.
    T: float - temperature in Kelvin, the temperature points to calculate the probability profile.
    model: Thermodynamic parameter sets to describe stacking energy. "SL_oligo" for the DNA oligo model in SantaLucia PNAS 1998, or "SL_polymer" for the DNA polymer model in SantaLucia PNAS 1998. Default is "SL_oligo"
    c_mono: Salt concentration, salt correction is given by SantaLucia PNAS 1998 
    
    Output:
    ========================================================================
    p_closed: 1D np.array - closed base probability values corresponding to each base of the input DNA sequence.
    """
    
    
    V10_LR,  U01_LR, U11_LR, Q_total_LR, s010_dict, send_dict = recursion(sequence, 
                                                                          T = T,
                                                                          model = model,
                                                                          c_mono = c_mono)
    V10_RL,  U01_RL, U11_RL, Q_total_RL,_ , _ = recursion(sequence.reverse_complement(), 
                                                          T,
                                                          model = model,
                                                          c_mono = c_mono)

    ## calculate p_closed for locations 2 to N - 1 (Tostesen, 2003 indexing)

    p_term_1 = np.array(list(U01_LR.values()))[:-1] * np.array(list(s010_dict.values()))[1:-1] * np.array(list(U01_RL.values()))[:-1][::-1]

    p_term_2 = np.array(list(U01_LR.values()))[:-1] * np.array(list(send_dict.values()))[1:-1] * np.array(list(U11_RL.values()))[:-1][::-1]

    p_term_3 = np.array(list(U11_LR.values()))[:-1] * np.array(list(send_dict.values()))[1:-1] * np.array(list(U01_RL.values()))[:-1][::-1]

    p_term_4 = np.array(list(U11_LR.values()))[:-1]  * np.array(list(U11_RL.values()))[:-1][::-1]

    p_closed_2_Nm1 = (p_term_1 + p_term_2 + p_term_3 + p_term_4 )/beta/Q_total_LR

    p_closed_1 = V10_RL[len(sequence) + 1]/Q_total_LR
    p_closed_N = V10_LR[len(sequence) + 1]/Q_total_LR

    p_closed = np.hstack([p_closed_1, p_closed_2_Nm1 , p_closed_N])

   

    return p_closed