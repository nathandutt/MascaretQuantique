import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def load_csv(filename):
    """Loads a CSV file without a header. Returns times and values."""
    data = np.genfromtxt(filename, delimiter=',', skip_header=0)
    t = data[:, 0]    # first column = time
    y = data[:, 1:]   # rest = values as a function of x
    return t, y

# Load data
t_plus, y_plus = load_csv('inverted_plus.csv')
t_minus, y_minus = load_csv('inverted_minus.csv')

x = np.arange(y_plus.shape[1])  # x-axis is just indices: 0,1,2,...
x = 0.1 * x - 10.
# Set up figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
line_plus, = ax1.plot([], [], 'b-')
line_minus, = ax2.plot([], [], 'r-')

ax1.set_xlim(-10, 10)
ax1.set_ylim(np.nanmin(y_plus), np.nanmax(y_plus))
ax1.set_title('Lambda Plus')
ax1.set_xlabel('x')
ax1.set_ylabel('Value')

ax2.set_xlim(-10, 10)
ax2.set_ylim(np.nanmin(y_minus), np.nanmax(y_minus))
ax2.set_title('Lambda Minus')
ax2.set_xlabel('x')
ax2.set_ylabel('Value')

# Update function for animation
def update(frame_idx):
    line_plus.set_data(x, y_plus[frame_idx])
    line_minus.set_data(x, y_minus[frame_idx])
    fig.suptitle(f'Time = {t_plus[frame_idx]:.3f}', fontsize=16)
    return line_plus, line_minus

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(t_plus), interval=50, blit=True)

plt.tight_layout()
plt.show()

