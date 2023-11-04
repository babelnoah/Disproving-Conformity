import numpy as np
import matplotlib.pyplot as plt

# Create a range of r values
r_values = np.linspace(0,1, 10)
#print(r_values)

#Create list of convergent values
final_x_one = []
final_x_two = []
final_delta_D = []
final_d = []

#initial conditions
x = np.zeros(2)
x[0] = 0.2
x[1] = 0.4
k = 0.5
delta = 0.5 
epsilon = 0.7

def calculate_next_generation(x, k, delta, epsilon):
    x_next = np.zeros(2)
    x_next[0] = x[0] * (1 - k/2 + k*delta/2) + x[1] * (k/2 + k*delta/2) + k*delta*epsilon
    x_next[1] = x[0] * (k/2 + k*delta/2) + x[1] * (1 - k/2 + k*delta/2) - k*delta*epsilon
    
    delta_x_one = x_next[0] - x[0]
    delta_x_two = x_next[1] - x[1]
    return x_next, delta_x_one, delta_x_two

def check_convergence(delta, k, epsilon, x, r):
    stable_generations = 0  # Initialize counter for stable generations
    change_threshold = 10**-15  # Set change threshold
    max_generations = 1.3*10**5  # Maximum number of generations to avoid infinite loop

    num_generations = 0
    delta_difference = []
    difference = []

    delta_D_old = 0  # Initial delta_D value
    d_old = 0

    while num_generations < max_generations:
        try:
            num_generations +=1
            x_old = x
            x, delta_x_one, delta_x_two = calculate_next_generation(x,k,delta,epsilon)
            delta_D = delta_x_one - delta_x_two
            d = x[1] - x[0]


            delta_difference.append(delta_D)
            difference.append(d)

            # Check if absolute change in all components of x and delta_D is less than threshold
            if np.all(np.abs(x - x_old) < change_threshold) and np.abs(delta_D - delta_D_old) and np.abs(d-d_old) < change_threshold:
                stable_generations += 1
            else:
                # Reset counter if change is above threshold
                stable_generations = 0

            if stable_generations >= 100:
                break
            delta_D_old = delta_D
            d_old = d
        except OverflowError:
            return x, delta_D,d # Return infinity for both x_c and x_n if OverflowError is encountered

    if num_generations >= max_generations:
        print(f'OverflowError at r = {r}, num_generations = {num_generations}, x_one = {x[0]}, x_two = {x[1]}, delta_D = {delta_D},d = {d}')
        return x, delta_D,d # Return infinity if max_generations is reached
    return x, delta_D,d

# Calculate the status of the function for every r value
for r in r_values:
    x = [0.2, 0.8]  # Reset initial conditions
    k = 0.5
    delta = 0.5 
    epsilon = 0.7
    x, delta_D_final,d = check_convergence(delta, k, epsilon, x, r)
    print(x)
    if np.isinf(x).any():  # If x contains any infinite values
        final_x_one.append('inf')
        final_x_two.append('inf')
    else:
        final_x_one.append(x[0])
        final_x_two.append(x[1])
    final_delta_D.append(delta_D_final)
    final_d.append(d)


print(final_delta_D)
# print(final_x_c)
# print(final_x_n)

# Plotting the graphs
plt.figure(figsize=(15, 5))

# plt.subplot(1,3,1)
# plt.plot(r_values, final_x_c, label='x_c')
# plt.ylim([0, 100])  # set the ylim to bottom, top
# plt.xlabel('r')
# plt.ylabel('Convergent x_c')
# plt.title('Convergent x_c as a function of r')

# plt.subplot(1,3,2)
# plt.plot(r_values, final_x_n, label='x_n')
# plt.ylim([0, 10])  # set the ylim to bottom, top
# plt.xlabel('r')
# plt.ylabel('Convergent x_n')
# plt.title('Convergent x_n as a function of r')

# # New subplot for delta D
# plt.subplot(1,3,3)
# plt.plot(r_values, final_delta_D, label='delta_D')
# plt.xlabel('r')
# plt.ylabel('Convergent delta_D')
# plt.title('Convergent delta_D as a function of r')

# New subplot for d
plt.plot(r_values, final_d, label='delta_D')
plt.xlabel('r')
plt.ylabel('Convergent D')
plt.title('Convergent D as a function of r')

plt.tight_layout()
print("")
print(k*delta*epsilon/(1-delta))
plt.show()

