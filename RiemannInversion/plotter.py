import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib

# Explicitly set backend (use 'TkAgg' for interactive windows)
matplotlib.use('TkAgg')

# === CONFIG ===
l_0 = np.sqrt(2)

# === HELPERS ===
def read_csv(filename):
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        data = []
        for row in reader:
            if not row:
                continue
            data.append([float(cell) if cell != 'NaN' else np.nan for cell in row])
        return data

def clip_lambdas(lambda_data):
    return [np.clip(np.array(row), -l_0, l_0) for row in lambda_data]

# === LOAD DATA ===
lambda_plus_raw = read_csv('inverted_plus.csv')
lambda_minus_raw = read_csv('inverted_minus.csv')

times = [row[0] for row in lambda_plus_raw]
lambda_plus = clip_lambdas([row[1:] for row in lambda_plus_raw])
lambda_minus = clip_lambdas([row[1:] for row in lambda_minus_raw])

# === PROCESS COMBINATIONS ===
lambda_sum = [lp + lm for lp, lm in zip(lambda_plus, lambda_minus)]
lambda_diff = [lp - lm for lp, lm in zip(lambda_plus, lambda_minus)]

x = np.arange(len(lambda_sum[0]))  # spatial axis

# === DEBUGGING: CHECK DATA ===
print(f"Times: {len(times)}")
print(f"lambda_sum: {len(lambda_sum)}")
print(f"lambda_diff: {len(lambda_diff)}")

# === SETUP PLOT ===
fig, ax = plt.subplots(figsize=(10, 4))
line_sum, = ax.plot([], [], label='λ₊ + λ₋', lw=2)
line_diff, = ax.plot([], [], label='λ₊ − λ₋', lw=2)
title = ax.set_title("")
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-2 * l_0, 2 * l_0)  # adjust if needed
ax.set_xlabel('x')
ax.set_ylabel('value')
ax.grid(True)
ax.legend()

# === ANIMATION ===
def update(frame):
    y_sum = lambda_sum[frame]
    y_diff = lambda_diff[frame]
    mask = ~np.isnan(y_sum) & ~np.isnan(y_diff)

    line_sum.set_data(x[mask], y_sum[mask])
    line_diff.set_data(x[mask], y_diff[mask])
    title.set_text(f"Time = {times[frame]:.2f}")
    return line_sum, line_diff, title

ani = animation.FuncAnimation(fig, update, frames=len(times), interval=100, blit=True)

# Show the plot
plt.tight_layout()
plt.show()

