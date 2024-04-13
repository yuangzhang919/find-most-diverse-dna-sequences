from utilities import helicity
from Bio.Seq import Seq
from matplotlib import pyplot as plt
from multiprocessing import Pool
import numpy as np
import os
import time
import random

# Function to generate a random DNA sequence with a given length range
def generate_a_random_DNA_sequence(min_length, max_length):
    bases = ['A', 'T', 'C', 'G']
    length = random.randint(min_length, max_length)
    sequence = ''.join(random.choice(bases) for _ in range(length))
    sequence = Seq(sequence)
    return sequence

def mean_squared_error(curve1, curve2):
    return np.mean((np.array(curve1) - np.array(curve2)) ** 2)

def find_most_diverse_sequences_simulated_annealing(curves, num_most_diverse_sequences, num_sequences, max_iterations, initial_temperature, cooling_rate):
    current_solution = random.sample(range(num_sequences), num_most_diverse_sequences)
    current_min_diff = min(mean_squared_error(curves[i], curves[j]) for i in current_solution for j in current_solution if i != j)

    best_solution = current_solution
    best_min_diff = current_min_diff

    temperature = initial_temperature

    min_diversity_history = [best_min_diff]

    for _ in range(max_iterations):
        new_solution = current_solution[:]
        while True:
            idx = random.randint(0, num_most_diverse_sequences - 1)
            new_idx = random.randint(0, num_sequences - 1)
            if new_idx not in new_solution:
                new_solution[idx] = new_idx
                break
        
        new_min_diff = min(mean_squared_error(curves[i], curves[j]) for i in new_solution for j in new_solution if i != j)

        if new_min_diff > current_min_diff or random.random() < np.exp((new_min_diff - current_min_diff) / temperature):
            current_solution = new_solution
            current_min_diff = new_min_diff

            if current_min_diff > best_min_diff:
                best_solution = current_solution
                best_min_diff = current_min_diff

        temperature *= cooling_rate
        min_diversity_history.append(best_min_diff)

    result_sequences = [sequences[idx] for idx in best_solution]
    return result_sequences, min_diversity_history

if __name__ == "__main__":
    start_time = time.time()  # Record the start time

    num_most_diverse_sequences = 10
    num_sequences = 100
    min_length = 25
    max_length = 25
    T_low = 20
    T_high = 100
    Ts = np.arange(T_low, T_high+ 0.5, 1) + 273
    max_iterations=10000
    initial_temperature=100
    cooling_rate=0.98

    # Generate random DNA sequences
    sequences = [generate_a_random_DNA_sequence(min_length, max_length) for _ in range(num_sequences)]
    generate_time = time.time() - start_time
    print(f"Time for generating random DNA sequences: {generate_time:.2f} seconds")

    # Parallelize the calculation of melting curves using multiprocessing
    with Pool(8) as pool:
        curves = pool.starmap(helicity, [(seq, Ts, "SL_oligo", 0.02) for seq in sequences])
    calculate_time = time.time() - start_time - generate_time
    print(f"Time for calculating melting curves: {calculate_time:.2f} seconds")
    ############################################################################################################

    most_diverse_sequences, min_diversity_history = find_most_diverse_sequences_simulated_annealing(curves, num_most_diverse_sequences, num_sequences, max_iterations, initial_temperature, cooling_rate)

    end_time = time.time()  # Record the end time
    running_time = end_time - start_time
    print(f"Simulated annealing algorithm running time: {(running_time - calculate_time - generate_time):.2f} seconds")
    print(f"Total running time: {(running_time):.2f} seconds")

    # Find the two sequences with the minimum difference
    min_diff = float('inf')
    min_diff_seq_indices = None
    for i in range(len(most_diverse_sequences)):
        for j in range(i+1, len(most_diverse_sequences)):
            diff = mean_squared_error(helicity(most_diverse_sequences[i], Ts, model="SL_oligo", c_mono=0.02),
                                      helicity(most_diverse_sequences[j], Ts, model="SL_oligo", c_mono=0.02))
            if diff < min_diff:
                min_diff = diff
                min_diff_seq_indices = (i, j)

    print(f"The two sequences with the minimum difference are:")
    print(f"Sequence {min_diff_seq_indices[0]+1}: {most_diverse_sequences[min_diff_seq_indices[0]]}")
    print(f"Sequence {min_diff_seq_indices[1]+1}: {most_diverse_sequences[min_diff_seq_indices[1]]}")
    print(f"The minimum difference is: {min_diff:.4f}")

    # Create a new folder and save results
    current_time = time.localtime()
    time_string = time.strftime("%Y%m%d%H%M%S", current_time)
    folder_name = "results_" + time_string
    os.makedirs(folder_name, exist_ok=True)

    # Save sequences to a fasta file
    fasta_file = os.path.join(folder_name, "most_diverse_sequences.fasta")
    with open(fasta_file, "w") as f:
        for i, seq in enumerate(most_diverse_sequences, start=1):
            f.write(f">Sequence_{i}\n")
            f.write(f"{seq}\n")

    # Plot and save melting curves for most_diverse_sequences
    fig, ax = plt.subplots(figsize=(8, 5))

    # Define a list of distinct colors
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

    for i, seq in enumerate(most_diverse_sequences, start=1):
        helicities = helicity(seq, Ts, model="SL_oligo", c_mono=0.02)
        ax.plot(Ts - 273, helicities, 'o-', color=colors[(i-1) % len(colors)], label=f"Sequence {i}")

    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("Helicity (%)")
    ax.set_ylim(0, 1)
    ax.set_title(f"Melting Curves for Top '{num_most_diverse_sequences}' Most Diverse Sequences")
    ax.legend(loc='best')  # Add a legend to the plot
    plt.tight_layout()
    plt.savefig(os.path.join(folder_name, "melting_curves.png"))  # Save the figure in the new folder

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(len(min_diversity_history)), min_diversity_history, '-', color='blue')
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Min Diversity")
    ax.set_title("Min Diversity over Iterations")
    plt.tight_layout()
    plt.savefig(os.path.join(folder_name, "min_diversity_over_iterations.png"))