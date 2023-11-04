import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def calculate_next_generation_vector_field(x_c_initial, x_n_initial, r, p, delta_c, delta_n, epsilon):
    delta_x_c = r*(1-p)*(x_n_initial - x_c_initial) + r*delta_c*np.sqrt(p*(1-p))*(x_n_initial - x_c_initial) + r*delta_c*epsilon
    delta_x_n = r*(-p)*(x_n_initial - x_c_initial) + r*delta_n*np.sqrt(p*(1-p))*(x_n_initial - x_c_initial) + r*delta_n*epsilon
    
    # Calculate next positions
    x_c_next = x_c_initial + delta_x_c
    x_n_next = x_n_initial + delta_x_n
    
    # Normalize the positions
    total = x_c_next + x_n_next
    x_c_next /= total
    x_n_next /= total

    return x_c_next - x_c_initial, x_n_next - x_n_initial

def iterate_system(x_c, x_n, r, p, delta_c, delta_n, epsilon, steps=10):
    for _ in range(steps):
        dx_c, dx_n = calculate_next_generation_vector_field(x_c, x_n, r, p, delta_c, delta_n, epsilon)
        x_c += dx_c
        x_n += dx_n
        
        # Normalizing
        total = x_c + x_n
        x_c /= total
        x_n /= total

    return x_c, x_n

r_val = 0.5
p_val = 0.2
delta_n_val = 0.1
epsilon_val = 0.99

x = np.linspace(0, 1, 20)
y = np.linspace(0, 1, 20)
X, Y = np.meshgrid(x, y)

fig, ax = plt.subplots(figsize=(8, 6))
Q = ax.quiver(X, Y, np.zeros(X.shape), np.zeros(Y.shape), angles='xy')

def update(delta_c_val):
    u, v = np.zeros(X.shape), np.zeros(Y.shape)
    converge_points_x, converge_points_y = [], []
    NI, NJ = X.shape
    
    for i in range(NI):
        for j in range(NJ):
            x_c = X[i, j]
            x_n = Y[i, j]
            dx_c, dx_n = calculate_next_generation_vector_field(x_c, x_n, r_val, p_val, delta_c_val, delta_n_val, epsilon_val)
            u[i,j], v[i,j] = dx_c, dx_n

            # Check where the system converges after several iterations
            x_c_end, x_n_end = iterate_system(x_c, x_n, r_val, p_val, delta_c_val, delta_n_val, epsilon_val)
            converge_points_x.append(x_c_end)
            converge_points_y.append(x_n_end)

    Q.set_UVC(u, v)
    ax.scatter(converge_points_x, converge_points_y, color='red', s=5)  # Highlight where points converge in red
    ax.set_title(f"Phase Portrait for r={r_val}, p={p_val}, delta_c={delta_c_val:.3f}, delta_n={delta_n_val}, epsilon={epsilon_val}")
    return Q,

ani = FuncAnimation(fig, update, frames=np.linspace(0.089, 0.102, 10), blit=True)
plt.xlabel('x_c')
plt.ylabel('x_n')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.show()