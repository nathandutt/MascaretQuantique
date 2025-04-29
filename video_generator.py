import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
from matplotlib.animation import FuncAnimation, FFMpegWriter
from csv_chomper import lambda_plus, lambda_minus, lambda_plus_theo, lambda_moins_theo, lambda_to_urho
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
L = config.get('L', 20.0)

x = np.arange(-x_range, x_range, dx)
x = x[:int(2 * x_range / dx)]

def rho_theoretical(x, t):
    lambda_plus_val = lambda_plus_theo(L, gamma, t, x)
    lambda_moins_val = lambda_moins_theo(L, gamma, t, x)

    r, u = lambda_to_urho(gamma, lambda_plus_val, lambda_moins_val)
    return r

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
    y3 = rho_theoretical(x, time)
    frames.append((time, y1, y2, y3))

# --- Set up figure ---
fig, ax = plt.subplots()
line1, = ax.plot(x, np.zeros_like(x), label='|ψ|²' if not args.use_lambda else 'λ⁺', linewidth=2.0, antialiased=True)
line2, = ax.plot(x, np.zeros_like(x), label='u'     if not args.use_lambda else 'λ⁻', linestyle='--', linewidth=2.0, antialiased=True)
line3, = ax.plot(x, np.zeros_like(x), label='ρ_theoretical', linewidth=2.0, linestyle=':', color='gray', antialiased=True)

time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

ax.set_xlim(x[0], x[-1])
ax.set_ylim(-3, 3)  # adjust if needed
ax.set_xlabel('x')
ax.set_ylabel('Field values')
ax.legend()

# --- Animation update function ---
def update(frame):
    time, y1, y2, y3 = frame
    line1.set_ydata(y1)
    line2.set_ydata(y2)
    line3.set_ydata(0)
    time_text.set_text(f'time = {time:.2f}')
    if time < 105 and time > 102: fig.savefig("test_frame.png", dpi=1000)
    return line1, line2, line3, time_text

# --- Create animation ---
ani = FuncAnimation(fig, update, frames=frames, blit=True)

# --- Save to video ---
writer = FFMpegWriter(fps=10, bitrate=5000)
ani.save(args.out, writer=writer,dpi=300)
print(f"Saved animation to {args.out}")

