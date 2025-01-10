import pandas as pd
import matplotlib.pyplot as plt

# Function to read and process data from a CSV file
def read_convergence_data(file_name):
    data = pd.read_csv(file_name)
    data['walltime'] = data['walltime'].apply(lambda x: float(x.replace('s', '')))
    return data

# Files for convergence results
files = ["convergence_results_1000.csv", "convergence_results_400.csv"]

# Loop over the files and plot the data
for file in files:
    data = read_convergence_data(file)

    # Separate data for convergence (CV) and divergence (DV)
    converged = data[data['status'] == 'CV']
    diverged = data[data['status'] == 'DV']

    # Plot iterations vs omega
    plt.figure(figsize=(10, 6))
    plt.plot(converged['omega'], converged['iterations'], marker='o', linestyle='-', color='blue', label='Converged (CV)')
    plt.scatter(diverged['omega'], [0] * len(diverged), color='red', label='Diverged (DV)', marker='x')
    plt.xlabel(r'$\omega$')
    plt.ylabel('Number of Iterations')
    plt.title(f'Convergence Analysis (Iterations) - {file}')
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)  # Add horizontal line for divergence
    plt.grid(True)
    plt.legend()
    plt.savefig(f"Convergence_Iterations_{file.split('.')[0]}.png")
    plt.clf()
