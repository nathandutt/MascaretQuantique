import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
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
args = parser.parse_args()

config = read_config('config.txt')
x_range = config.get('x_range', 20)
dx      = config.get('dx',    0.1)
gamma   = config.get('exponent', 2.0)  # required for lambda±
L = config.get('L', 10.0)

x = np.arange(-x_range, x_range, dx)
x = x[:int(2 * x_range / dx)]

filename = "evolution.csv"
stride   = 1

fig, ax = plt.subplots()
line1, = ax.plot(x, np.zeros_like(x), label='|ψ|²' if not args.use_lambda else 'λ⁺')
line2, = ax.plot(x, np.zeros_like(x), label='u'     if not args.use_lambda else 'λ⁻', linestyle='--')
#line3, = ax.plot(x, np.zeros_like(x), label='rho_theoretical', linewidth = 1.0, color = 'blue')

time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

ax.set_xlim(x[0], x[-1])
ax.set_ylim(-3, 3)  # Fixed for consistent scale
ax.set_xlabel('x')
ax.set_ylabel('Field values')
ax.legend()

# Define rho_theoretical function
def rho_theoretical(x, t):
    lambda_plus_val = lambda_plus_theo(L, gamma, t, x)
    lambda_moins_val = lambda_moins_theo(L, gamma, t, x)

    r, u = lambda_to_urho(gamma, lambda_plus_val, lambda_moins_val)
    return r

with open(filename, newline='') as csvfile:
    reader   = csv.reader(csvfile)
    row_iter = iter(reader)
    i = 0
    for row in row_iter:
        if not plt.fignum_exists(fig.number):
            break

        try:
            row2 = next(row_iter)
        except StopIteration:
            break

        if i % stride != 0:
            i += 1
            continue

        row = [float(v) for v in row]
        time = row[0]
        real = row[1::2]
        imag = row[2::2]

        rho = np.array(real)**2 + np.array(imag)**2

        row2 = [float(v) for v in row2]
        u = np.array(row2[1:])

        if args.use_lambda:
            y1 = lambda_plus(rho, u, gamma)
            y2 = lambda_minus(rho, u, gamma)
        else:
            y1 = rho
            y2 = u
        
        # Update line3 with rho_theoretical
        #y3 = rho_theoretical(x, time)
        #line3.set_ydata(y3)
        line1.set_ydata(y1)
        line2.set_ydata(y2)
        time_text.set_text(f'time = {time:.2f}')

        if plt.waitforbuttonpress(timeout=-1) is None:
            break

        i += 1

plt.close(fig)

