import numpy as np
import matplotlib.pyplot as plt

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

r_val = 0.5
p_val = 0.2
delta_c_val = 0.095
delta_n_val = 0.1
epsilon_val = 0.99

# Generate a grid of points
x = np.linspace(0, 1, 20)
y = np.linspace(0, 1, 20)

X, Y = np.meshgrid(x, y)
u, v = np.zeros(X.shape), np.zeros(Y.shape)

NI, NJ = X.shape
for i in range(NI):
    for j in range(NJ):
        x_c = X[i, j]
        x_n = Y[i, j]
        dx_c, dx_n = calculate_next_generation_vector_field(x_c, x_n, r_val, p_val, delta_c_val, delta_n_val, epsilon_val)
        u[i,j], v[i,j] = dx_c, dx_n

plt.quiver(X, Y, u, v, angles='xy')
plt.xlabel('x_c')
plt.ylabel('x_n')
plt.title(f"Phase Portrait for r={r_val}, p={p_val}, delta_c={delta_c_val}, delta_n={delta_n_val}, epsilon={epsilon_val}")
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.show()

