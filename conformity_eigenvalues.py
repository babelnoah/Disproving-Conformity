import numpy as np
import matplotlib.pyplot as plt

def jacobian_matrix(r, p, delta_c, delta_n, epsilon):
    A11 = r*(1-2*p) + r*delta_c*np.sqrt(p*(1-p)) + r*delta_c*epsilon
    A12 = r*(p-1) - r*delta_c*np.sqrt(p*(1-p))
    A21 = r*(p-1) - r*delta_n*np.sqrt(p*(1-p))
    A22 = r*(1-2*p) + r*delta_n*np.sqrt(p*(1-p)) + r*delta_n*epsilon
    
    A = np.array([[A11, A12], [A21, A22]])
    
    return A

# Define the function to compute the eigenvalues of the Jacobian matrix
def compute_eigenvalues(r, p, delta_c, delta_n, epsilon):
    J = jacobian_matrix(r, p, delta_c, delta_n, epsilon)
    return np.linalg.eigvals(J)
# Initialize some fixed parameter values for reference
r_fixed = 0.5
p_fixed = 0.2
delta_c_fixed = 0.61
delta_n_fixed = 0.1
epsilon_fixed = 0.99


# Define a range for each parameter
p_values = np.linspace(0.1, 0.9, 100)
delta_c_values = np.linspace(0, 1, 100)
delta_n_values = np.linspace(0, 1, 100)
epsilon_values = np.linspace(0, 1.5, 100)

parameters = [
    ("p", p_values, p_fixed),
    ("delta_c", delta_c_values, delta_c_fixed),
    ("delta_n", delta_n_values, delta_n_fixed),
    ("epsilon", epsilon_values, epsilon_fixed)
]

# Create plots for each parameter
for param_name, param_values, param_fixed in parameters:
    eigen1_list = []
    eigen2_list = []

    for value in param_values:
        if param_name == "p":
            eigenvalues = compute_eigenvalues(r_fixed, value, delta_c_fixed, delta_n_fixed, epsilon_fixed)
        elif param_name == "delta_c":
            eigenvalues = compute_eigenvalues(r_fixed, p_fixed, value, delta_n_fixed, epsilon_fixed)
        elif param_name == "delta_n":
            eigenvalues = compute_eigenvalues(r_fixed, p_fixed, delta_c_fixed, value, epsilon_fixed)
        else:  # epsilon
            eigenvalues = compute_eigenvalues(r_fixed, p_fixed, delta_c_fixed, delta_n_fixed, value)

        eigen1_list.append(eigenvalues[0])
        eigen2_list.append(eigenvalues[1])

    plt.figure(figsize=(10, 6))
    plt.plot(param_values, eigen1_list, label='Eigenvalue 1', color='blue')
    plt.plot(param_values, eigen2_list, label='Eigenvalue 2', color='red')
    plt.xlabel(param_name)
    plt.ylabel('Eigenvalue')
    plt.legend()
    plt.title(f'Eigenvalues as a function of {param_name}')
    plt.grid(True)
    plt.show()