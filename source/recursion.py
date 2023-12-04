#/usr/bin/env python3 

from stability_factors import stability_factor_lookup
from loop_entropy import omega, W
from thermodynamics_parameters import R, beta, omega0, alpha


def recursion(sequence, 
              T, 
              model = "SL_oligo",
              c_mono = 0.02):
    
    """
    Function to calculate the partition functions of all subchains of the forward and reverse passes. This function is used in the p_closed_profile function in the utilities profile to further calculate the per-base probability for being in the closed state. This function is an implementation of the Tostesen, 2003 recursion algorithm for partition function calculation of the DNA melting problem. 
    
    Input:
    ========================================================================
    sequence: str - the sequence of DNA to calculate the partition functions
    T: float - temperature in Kelvin, the temperature to calculate partition functions
    model: Thermodynamic parameter sets to describe stacking energy. "SL_oligo" for the DNA oligo model in SantaLucia PNAS 1998, or "SL_polymer" for the DNA polymer model in SantaLucia PNAS 1998. Default is "SL_oligo"
    c_mono: Salt concentration, salt correction is given by SantaLucia PNAS 1998 
    
    Output:
    ========================================================================
    p_closed: 1D np.array - closed base probability values corresponding to each base of the input DNA sequence.
    """
    
    ## Calculate stability factors based on the choice of thermodynamic stacking parameter set and salt correction
    s010_lookup, send_lookup, s11_lookup = stability_factor_lookup(T, 
                                                                   c_mono = c_mono, 
                                                                   model = model, 
                                                                   N = len(sequence) - 1)
    

    ## Initialize recursion
    V10 = {}
    V10[1] = 1 ## V10(i + 1) = V10(1), i = 0, looking at the sub sequence [0]
    
    s010_dict = {}
    s010_dict[1] = s010_lookup[sequence[1 - 1]]
    V10[2] = beta*s010_dict[1]  ## V10(i + 1) = V10(2), loc is consistent with tostensen indexing, i = 1, looking at the sub sequence [01]
 
    U01 = {}
    U01[2] = beta  ## U01(2)
    U11 = {}
    send_dict = {}
    send_dict[1] = send_lookup[sequence[1 - 1]]
    U11[2] = beta*send_dict[1]*s11_lookup[sequence[2-2:2]] ## U11(2), loc is consistent with tostesen indexing
    s010_dict[2] = s010_lookup[sequence[2 - 1]]
    send_dict[2] = send_lookup[sequence[2 - 1]]
    V10[3] = s010_dict[2]*U01[2] + send_dict[2]* U11[2]

    Q_total = V10[1] + V10[2] + V10[3]

    ## Recursion
    for i in range(3, len(sequence)+1):
        s010_dict[i] = s010_lookup[sequence[i - 1]]
        send_dict[i] = send_lookup[sequence[i - 1]]
        
        U01[i] = beta*V10[1] + W(i, V10)
        U11[i] = s11_lookup[sequence[i-2:i]]*(U01[i - 1]*send_dict[i-1] + U11[i - 1])
        V10[i+1] = s010_dict[i]*U01[i]+ send_dict[i]* U11[i]

        Q_total += V10[i+1]
        
    return V10, U01, U11, Q_total, s010_dict, send_dict


