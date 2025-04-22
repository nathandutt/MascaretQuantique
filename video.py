import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.animation import FuncAnimation, FFMpegWriter

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
x_range = config.get('x_range', 20)  # Default to 20 if not in config
dx = config.get('dx', 0.1)  # Default to 0.1 if not in config
x = np.arange(-x_range, x_range, dx)
x = x[:int(2 * x_range / dx)]

# Prepare plot
fig, ax = plt.subplots()
line, = ax.plot(x, np.zeros_like(x))
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

# Set consistent y-limits
ax.set_ylim(0, 2)  # Adjust based on your data
ax.set_xlim(x[0], x[-1])
ax.set_xlabel('x')
ax.set_ylabel('|ψ|²')

# Prepare the writer for saving video
writer = FFMpegWriter(fps=30, metadata=dict(artist="Your Name", title="Evolution of Psi"))
filename = 'video_NLS.mp4'  # Name of the output video file

# Animation function
def update(frame):
    # Read the CSV file and plot each frame
    with open('evolution.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        for i, row in enumerate(reader):
            if i != frame:
                continue  # Skip until we get the correct frame
            
            row = [float(val) for val in row]
            time = row[0]
            real = row[1::2]
            imag = row[2::2]
            psi_squared = np.array(real)**2 + np.array(imag)**2
            
            # Update the plot with new data
            line.set_ydata(psi_squared)
            time_text.set_text(f'time = {time:.2f}')
            break

    return line, time_text

# Get the total number of time steps (frames) in the CSV
total_frames = sum(1 for row in open('evolution.csv')) - 1  # Subtract 1 for header

# Create the animation
ani = FuncAnimation(fig, update, frames=range(1, total_frames), interval=50, blit=True)

# Save the animation as a video without displaying it
ani.save(filename, writer=writer)

