import pandas as pd
import matplotlib.pyplot as plt

# Function to read and process data from a CSV file
def read_convergence_data(file_name):
    data = pd.read_csv(file_name)
    data['walltime'] = data['walltime'].apply(lambda x: float(x.replace('s', '')))
    return data

### SIZE 1000
# Files for SOR and Red-Black SOR results
files = {
    "SEQ-SOR 400": "convergence_seq_results_400.csv",
    "RB-SOR 400": "convergence_results_400.csv"
}

# Plot data
plt.figure(figsize=(12, 8))

for label, file in files.items():
    data = read_convergence_data(file)

    # Separate data for convergence (CV) and divergence (DV)
    converged = data[data['status'] == 'CV']
    diverged = data[data['status'] == 'DV']

    # Plot converged points (iterations vs omega)
    plt.plot(
        converged['omega'],
        converged['iterations'],
        marker='o',
        linestyle='-',
        label=f'{label} (Converged)'
    )

    # Plot diverged points (mark as 0 iterations)
    plt.scatter(
        diverged['omega'],
        [0] * len(diverged),
        color='red',
        label=f'{label} (Diverged)',
        marker='x'
    )

# Customize plot
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel('Number of Iterations', fontsize=14)
# plt.title('Comparison of SEQ-SOR and RB-SOR Solvers', fontsize=16)
plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)  # Horizontal line for divergence
plt.grid(True)
plt.legend()
plt.savefig("Comparison_SEQ-SOR_RB-SOR_400.png")
plt.clf()

### SIZE 1000
# Files for SOR and Red-Black SOR results
files = {
    "SEQ-SOR 400": "convergence_seq_results_1000.csv",
    "RB-SOR 400": "convergence_results_1000.csv"
}

# Plot data
plt.figure(figsize=(12, 8))

for label, file in files.items():
    data = read_convergence_data(file)

    # Separate data for convergence (CV) and divergence (DV)
    converged = data[data['status'] == 'CV']
    diverged = data[data['status'] == 'DV']

    # Plot converged points (iterations vs omega)
    plt.plot(
        converged['omega'],
        converged['iterations'],
        marker='o',
        linestyle='-',
        label=f'{label} (Converged)'
    )

    # Plot diverged points (mark as 0 iterations)
    plt.scatter(
        diverged['omega'],
        [0] * len(diverged),
        color='red',
        label=f'{label} (Diverged)',
        marker='x'
    )

# Customize plot
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel('Number of Iterations', fontsize=14)
# plt.title('Comparison of SEQ-SOR and RB-SOR Solvers', fontsize=16)
plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)  # Horizontal line for divergence
plt.grid(True)
plt.legend()
plt.savefig("Comparison_SEQ-SOR_RB-SOR_1000.png")
plt.clf()
