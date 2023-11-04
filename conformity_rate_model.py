import numpy as np
import matplotlib.pyplot as plt

# Create a range of r values
r_values = [0.5]

# initial conditions
p = 1/2
delta_c = 0
delta_n = 0.4
epsilon = 0.5
x_c = 0.2
x_n = 0.8

def calculate_next_generation(x_c_initial, x_n_initial, r, p, delta_n, delta_c, epsilon):
    delta_x_c = r*(1-p)*(x_n_initial - x_c_initial) + r*(delta_c)*(np.sqrt(p*(1-p)))*(x_n_initial - x_c_initial) + r*(delta_c)*(epsilon)
    delta_x_n = r*(-p)*(x_n_initial - x_c_initial) + r*(delta_n)*(np.sqrt(p*(1-p)))*(x_n_initial - x_c_initial) + r*(delta_n)*(epsilon)
    x_c_next = x_c_initial + delta_x_c
    x_n_next = x_n_initial + delta_x_n
    return x_c_next, x_n_next, delta_x_c, delta_x_n

def check_convergence(p, delta_c, delta_n, epsilon, x_c, x_n, r):
    stable_generations = 0
    change_threshold = 10**-15
    max_generations = 10**2
    num_generations = 0

    rate_of_change_x_c = [] # List to store rate of changes for x_c
    rate_of_change_x_n = [] # List to store rate of changes for x_n

    delta_D_old = 0
    d_old = 0

    while num_generations < max_generations:
        try:
            num_generations += 1
            x_c_old, x_n_old = x_c, x_n  # Store the old values before calculating new ones
            x_c, x_n, delta_x_c, delta_x_n = calculate_next_generation(x_c, x_n, r, p, delta_n, delta_c, epsilon)
            delta_D = delta_x_n - delta_x_c
            d = x_n - x_c

            if np.all(np.abs(np.array([x_c, x_n]) - np.array([x_c_old, x_n_old])) < change_threshold) and np.abs(delta_D - delta_D_old) and np.abs(d - d_old) < change_threshold:
                stable_generations += 1
            else:
                stable_generations = 0

            # Save the absolute change at each generation for x_c and x_n
            rate_of_change_x_c.append(np.abs(x_c - x_c_old))
            rate_of_change_x_n.append(np.abs(x_n - x_n_old))

            if stable_generations >= 100:
                break
            delta_D_old = delta_D
            d_old = d
        except OverflowError:
            return np.inf, np.inf, rate_of_change_x_c, rate_of_change_x_n
    if num_generations >= max_generations:
        return np.inf, np.inf, rate_of_change_x_c, rate_of_change_x_n
    return x_c, x_n, rate_of_change_x_c, rate_of_change_x_n

# Calculate the status of the function for every r value and plot the rate of change
plt.figure(figsize=(15, 10))
for r in r_values:
    x_c, x_n = 0.2, 0.8  # Reset initial conditions
    x_c_final, x_n_final, rate_of_change_x_c, rate_of_change_x_n = check_convergence(p, delta_c, delta_n, epsilon, x_c, x_n, r)

    # plot rate of change for each r
    plt.plot(range(len(rate_of_change_x_c)), rate_of_change_x_c, label='x_c, r = ' + str(r))
    plt.plot(range(len(rate_of_change_x_n)), rate_of_change_x_n, label='x_n, r = ' + str(r))

plt.xlabel('Generation')
plt.ylabel('Rate of change')
plt.legend()
plt.show()
