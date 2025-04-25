import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
from matplotlib.animation import FuncAnimation, FFMpegWriter
from csv_chomper import lambda_plus, lambda_minus

def read_config(filename):
    config = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) == 2:
                key, value = parts
                try:
                    config[key] = float(value)
                except ValueError:
                    print(f"Warning: Skipping invalid line (non-numeric value) in config: {line}")
            else:
                print(f"Warning: Skipping invalid line (expected key-value pair) in config: {line}")
    return config

# --- Command-line flags ---
parser = argparse.ArgumentParser()
parser.add_argument('--lambda', dest='use_lambda', action='store_true', help='Plot lambda+ and lambda- instead of u and |psi|^2')
parser.add_argument('--out', default='output.mp4', help='Output video filename')
args = parser.parse_args()

# --- Read config ---
config = read_config('config.txt')
x_range = config.get('x_range', 20)
dx      = config.get('dx',    0.1)
gamma   = config.get('exponent', 2.0)  # required for lambda±

x = np.arange(-x_range, x_range, dx)
x = x[:int(2 * x_range / dx)]

# --- Load all data first ---
filename = "evolution.csv"
frames = []

with open(filename, newline='') as csvfile:
    reader = csv.reader(csvfile)
    rows = list(reader)

for i in range(0, len(rows) - 1, 2):  # every pair of lines is a frame
    row1 = [float(v) for v in rows[i]]
    row2 = [float(v) for v in rows[i + 1]]

    time = row1[0]
    real = row1[1::2]
    imag = row1[2::2]
    rho = np.array(real)**2 + np.array(imag)**2

    u = np.array(row2[1:])

    if args.use_lambda:
        y1 = lambda_plus(rho, u, gamma)
        y2 = lambda_minus(rho, u, gamma)
    else:
        y1 = rho
        y2 = u

    frames.append((time, y1, y2))

# --- Set up figure ---
fig, ax = plt.subplots()
line1, = ax.plot(x, np.zeros_like(x), label='|ψ|²' if not args.use_lambda else 'λ⁺')
line2, = ax.plot(x, np.zeros_like(x), label='u'     if not args.use_lambda else 'λ⁻', linestyle='--')
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

ax.set_xlim(x[0], x[-1])
ax.set_ylim(-3, 3)  # adjust if needed
ax.set_xlabel('x')
ax.set_ylabel('Field values')
ax.legend()

# --- Animation update function ---
def update(frame):
    time, y1, y2 = frame
    line1.set_ydata(y1)
    line2.set_ydata(y2)
    time_text.set_text(f'time = {time:.2f}')
    return line1, line2, time_text

# --- Create animation ---
ani = FuncAnimation(fig, update, frames=frames, blit=True)

# --- Save to video ---
writer = FFMpegWriter(fps=10, bitrate=1800)
ani.save(args.out, writer=writer)
print(f"Saved animation to {args.out}")

