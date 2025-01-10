import pandas as pd
import matplotlib.pyplot as plt

# Function to read and extract data from a CSV file
def read_data(file_name):
    data = pd.read_csv(file_name)
    return data['tasks'], data['walltime'].apply(lambda x: float(x.replace('s', '')))

# Function to compute speed-up
def compute_speedup(tasks, walltime):
    T1 = walltime.iloc[0]  # Time with one task
    speedup = T1 / walltime
    return tasks, speedup

# Scalability comparison for a single CPU frequency
file_omega1_2_cpu2000000 = 'scalability_results_omega1.2_cpu2000000_time6_size1000.csv'
tasks, walltime = read_data(file_omega1_2_cpu2000000)

# Plot scalability (time vs. tasks)
plt.plot(tasks, walltime, marker='o', label="CPU 2000000", linestyle='-', markersize=6)
plt.xlabel('Number of tasks')
plt.ylabel('Time (s)')
plt.title('Scalability for CPU 2000000')
plt.grid(True)  # Add grid
plt.legend()
plt.savefig("Scalability_CPU2000000.png")
plt.clf()

# Plot speed-up for the single frequency
tasks, speedup = compute_speedup(tasks, walltime)
plt.plot(tasks, speedup, marker='o', label="CPU 2000000", linestyle='-', markersize=6)
plt.xlabel('Number of tasks')
plt.ylabel('Speed-up')
plt.title('Speed-up for CPU 2000000')
plt.grid(True)  # Add grid
plt.legend()
plt.savefig("Speedup_CPU2000000.png")
plt.clf()

# Comparison of CPU frequency impact
files = [
    'scalability_results_omega1.2_cpu2000000_time6_size1000.csv',
    'scalability_results_omega1.2_cpu2400000_time6_size1000.csv',
    'scalability_results_omega1.2_cpu2800000_time6_size1000.csv',
    'scalability_results_omega1.2_cpu3200000_time6_size1000.csv'
]

plt.figure(figsize=(10, 6))

# Plot time comparison for multiple frequencies
for file in files:
    tasks, walltime = read_data(file)
    plt.plot(tasks, walltime, marker='o', label=file.split('_')[3], linestyle='-', markersize=6)

plt.xlabel('Number of tasks')
plt.ylabel('Time (s)')
plt.title('Scalability for CPU (20-32)')
plt.grid(True)  # Add grid
plt.legend(title="CPU Frequency")
plt.savefig("Scalability_CPUS.png")
plt.clf()

# Plot speed-up comparison for multiple frequencies
plt.figure(figsize=(10, 6))
for file in files:
    tasks, walltime = read_data(file)
    tasks, speedup = compute_speedup(tasks, walltime)
    plt.plot(tasks, speedup, marker='o', label=file.split('_')[3], linestyle='-', markersize=6)

plt.xlabel('Number of tasks')
plt.ylabel('Speed-up')
plt.title('Speed-up for CPU (20-32)')
plt.grid(True)  # Add grid
plt.legend(title="CPU Frequency")
plt.savefig("Speedup_CPUS.png")
plt.clf()
