#!/usr/bin/env python3

from utilities import p_closed_profile, helicity
from Bio.Seq import Seq
from matplotlib import pyplot as plt
import numpy as np
import os
import pandas as pd
import time

import figure_setting

def create_results_folder(sequence):
    """
    Function to create the results folder once the user starts the program by inputting a DNA sequence.
    """
    # Define the folder name
    current_time = time.localtime()
    time_string = time.strftime("%Y%m%d%H%M%S", current_time)
    folder_name = "results_" + time_string

    # Get the current working directory
    current_directory = os.getcwd()

    # Create the full path for the new folder
    new_folder_path = os.path.join(current_directory, folder_name)

    # Check if the folder already exists before creating it
    if not os.path.exists(new_folder_path):
        # Create the new folder
        os.makedirs(new_folder_path)
        print(f"Folder '{folder_name}' created successfully.")
    else:
        print(f"The result folder 'results' already exists. The files in this folder will be overwritten")
        
    save_to_fasta(sequence, time_string)

    return folder_name

def save_to_fasta(data_string, time_string):
    
    # Define the file name
    file_name = f"results_{time_string}/{time_string}_DNA_melt.fasta"

    # Write the data string to the FASTA file
    with open(file_name, "w") as file:
        file.write(f">DNA_melt_{time_string}\n")
        file.write(str(data_string))


if __name__ == "__main__":
    
    ## Ask user to input a sequence. 
    sequence = input("Please input a oligo sequence").upper()
    sequence = Seq(sequence)

    
    ## Create results folder
    folder_name = create_results_folder(sequence)

    ## Plot probability profile under a given temperature
    # Ask user to input a temperature.
    T_plot_profile = float(input("Please input a temperature (Celcius) to plot the probability profile"))
    # Calculate probability profile
    probability_profile = p_closed_profile(sequence, 
                                           T_plot_profile+273.,
                                           model = "SL_oligo",
                                           c_mono = 0.02)
    # Plot and save probability profile
    fig, ax = plt.subplots(figsize = (8,5))
    ax.plot(probability_profile,'o-')
    ax.set_xlabel("NT_position")
    ax.set_ylabel("Closed state probability")
    ax.set_ylim(0,1)
    ax.set_title(f"Probability Profile at T = {T_plot_profile} C")
    plt.tight_layout()
    plt.savefig(folder_name + f"/T{T_plot_profile}C_probability_profile.png")
    pd.DataFrame({"NT position":np.arange(1, len(sequence)+1), 
                  "Closed state probability": probability_profile}).to_csv(folder_name + f"/T{T_plot_profile}C_probability_profile.csv")


    ## Plot melting curve
    # Ask user to input a temperature range
    T_low = float(input("Now please input the temperature range to plot the melting curve. \
    What is the lowest temperature you want to plot?"))
    T_high = float(input("What is the highest temperature you want to plot?"))
    Ts = np.arange(T_low, T_high+ 0.5, 0.5) + 273 ## change to Kelvin
    # Calculate melting curve as helicity values in the user-defined temperature range
    helicities = helicity(sequence, 
                          Ts,
                          model = "SL_oligo",
                          c_mono = 0.02)
    # Plot and save melting curve
    fig, ax = plt.subplots(figsize = (8,5))
    ax.plot(Ts - 273, helicities, 'o-') ## visualization, use Celcius
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("Helicity (%)")
    ax.set_ylim(0,1)
    ax.set_title(f"Melting Curve")
    plt.tight_layout()
    plt.savefig(folder_name + "/melting_curve.png")


    # Calculate minus derivative of melting curve
    dFdT = np.gradient(helicities, Ts)
    # Plot and save minus derivative of melting curve
    fig, ax = plt.subplots(figsize = (8,5))
    ax.plot(Ts - 273, -dFdT, 'o-') ## visualization, use Celcius for Ts
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("-dF/dT")
    plt.tight_layout()
    plt.savefig(folder_name + "/negative_derivative.png")

    pd.DataFrame({"Temperature (C)":Ts - 273., 
                  "Helicity": helicities,
                  "Minus Derivative (1/C)": -dFdT}).to_csv(folder_name + f"/melt_curve_derivative.csv")

    print("Calculation done")