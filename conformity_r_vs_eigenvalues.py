import numpy as np
import matplotlib.pyplot as plt

# Given parameters
p = 0.5
delta_c = 0
delta_n = 0.4
epsilon = 0.5

def compute_matrix(r):
    A11 = r*(1-2*p) + r*delta_c*np.sqrt(p*(1-p)) + r*delta_c*epsilon
    A12 = r*(p-1) - r*delta_c*np.sqrt(p*(1-p))
    A21 = r*(p-1) - r*delta_n*np.sqrt(p*(1-p))
    A22 = r*(1-2*p) + r*delta_n*np.sqrt(p*(1-p)) + r*delta_n*epsilon
    
    A = np.array([[A11, A12], [A21, A22]])
    
    return A

def compute_eigenvalues(A):
    eigenvalues, _ = np.linalg.eig(A)
    return eigenvalues

r_values = np.linspace(0.1,0.9, 100)  # increased resolution for better plotting
eigenvalue1_list = []
eigenvalue2_list = []

for r in r_values:
    A = compute_matrix(r)
    eigenvalues = compute_eigenvalues(A)
    eigenvalue1_list.append(eigenvalues[0].real)  # only taking the real part
    eigenvalue2_list.append(eigenvalues[1].real)

# Plotting
plt.figure(figsize=(10,6))
plt.plot(r_values, eigenvalue1_list, label="Eigenvalue 1", color='blue')
plt.plot(r_values, eigenvalue2_list, label="Eigenvalue 2", color='red')
plt.axhline(y=0, color='gray', linestyle='--')  # Stability threshold
plt.xlabel("r values")
plt.ylabel("Real part of Eigenvalues")
plt.title("Eigenvalues vs r")
plt.legend()
plt.grid(True)
plt.show()