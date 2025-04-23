import argparse
import matplotlib.pyplot as plt
from csv_chomper import *

import csv
import numpy as np


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
args = parser.parse_args()

config = read_config('config.txt')
x_range = config.get('x_range', 20)
dx      = config.get('dx',    0.1)
gamma   = config.get('exponent', 2.0)
dt      = config.get('dt',     0.0005)
t_max   = config.get('t_max', 10.0)
L       = config.get('L',      1.0)
print("Gamma = ", gamma)
x = np.arange(-x_range, x_range, dx)
x = x[:int(2 * x_range / dx)]

# --- Read CSV data ---
filename = 'evolution.csv'
u_data, rho_data, times = read_csv_fields(filename)

# --- Plot setup ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
line_rho, = ax1.plot([], [], label='ρ (numerical)')
line_rho_theo, = ax1.plot([], [], '--', label='ρ (theoretical)')

line_u, = ax2.plot([], [], label='u (numerical)')
line_u_theo, = ax2.plot([], [], '--', label='u (theoretical)')

ax1.set_ylabel('Density ρ')
ax2.set_ylabel('Velocity u')
ax2.set_xlabel('x')
ax1.legend()
ax2.legend()

fig.tight_layout()
fig.subplots_adjust(top=0.9)
title = fig.suptitle("")

# --- Animation loop ---
for t_index in range(len(times)):
    t = times[t_index]

    lambda_p = np.array([lambda_plus_theo(L, gamma, t, xi) for xi in x])
    lambda_m = np.array([lambda_moins_theo(L, gamma, t, xi) for xi in x])
    rho_theo, u_theo = lambda_to_urho(gamma, lambda_p, lambda_m)

    line_rho.set_data(x, rho_data[t_index])
    line_rho_theo.set_data(x, rho_theo)

    line_u.set_data(x, u_data[t_index])
    line_u_theo.set_data(x, u_theo)

    ax1.set_xlim(x[0], x[-1])
    ax1.set_ylim(-0.1, np.max(rho_data[t_index]) * 1.1)
    ax2.set_ylim(np.min(u_data[t_index]) - 0.5, np.max(u_data[t_index]) + 0.5)

    title.set_text(f"Comparison at time t = {t:.3f}")
    plt.pause(0.01)

    if plt.waitforbuttonpress(timeout=-1) is None:
        break

plt.close(fig)

