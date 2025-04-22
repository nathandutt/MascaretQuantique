import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read the configuration file and extract parameters
def read_config(filename):
    config = {}
    with open(filename, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespace and skip empty lines or comments
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            # Try to split the line into key and value
            parts = line.split()
            if len(parts) == 2:  # Make sure we have exactly two parts (key, value)
                key, value = parts
                try:
                    config[key] = float(value)
                except ValueError:
                    print(f"Warning: Skipping invalid line (non-numeric value) in config: {line}")
            else:
                print(f"Warning: Skipping invalid line (expected key-value pair) in config: {line}")
    return config


# Read parameters from the config file
config = read_config('config.txt')

# Parameters from config
x_range = config.get('x_range', 20)
dx = config.get('dx', 0.1)
x = np.arange(-x_range, x_range, dx)
x = x[:int(2 * x_range / dx)]

# Read data from CSV (skip 99 out of 100)
filename = "evolution.csv"
stride = 25

# Prepare plot
fig, ax = plt.subplots()
line, = ax.plot(x, np.zeros_like(x))
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

# Set consistent y-limits
ax.set_ylim(0, 2)  # Adjust based on your data
ax.set_xlim(x[0], x[-1])
ax.set_xlabel('x')
ax.set_ylabel('|ψ|²')

# Animation (read and plot every 100th line)
with open(filename, newline='') as csvfile:
    reader = csv.reader(csvfile)
    for i, row in enumerate(reader):
        if i % stride != 0:
            continue  # Skip 99/100 rows
        row = [float(val) for val in row]
        time = row[0]
        real = row[1::2]
        imag = row[2::2]
        psi_squared = np.array(real)**2 + np.array(imag)**2
        line.set_ydata(psi_squared)
        time_text.set_text(f'time = {time:.2f}')
        plt.waitforbuttonpress()

plt.show()

