import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Create a range of r values
r_values = np.linspace(0.0001,0.99999, 5)

#initial conditions
p = 1/2
delta_c = 0
delta_n = 0.4
epsilon = 0.5
x_c = 0.2
x_n = 0.8

def calculate_next_generation(x_c_initial, x_n_initial, r, p,delta_n,delta_c,epsilon):
    delta_x_c = r*(1-p)*(x_n_initial - x_c_initial) + r*(delta_c)*(np.sqrt(p*(1-p)))*(x_n_initial - x_c_initial) + r*(delta_c)*(epsilon)
    delta_x_n = r*(-p)*(x_n_initial - x_c_initial) + r*(delta_n)*(np.sqrt(p*(1-p)))*(x_n_initial - x_c_initial) + r*(delta_n)*(epsilon)
    x_c_next = x_c_initial + delta_x_c
    x_n_next = x_n_initial + delta_x_n
    return x_c_next, x_n_next

def check_convergence(p, delta_c, delta_n, epsilon, x_c, x_n, r):
    max_generations = 10**3
    num_generations = 0
    x_c_values = []
    x_n_values = []
    while num_generations < max_generations:
        try:
            num_generations +=1
            x_c, x_n = calculate_next_generation(x_c, x_n, r, p, delta_n, delta_c, epsilon)
            x_c_values.append(x_c)
            x_n_values.append(x_n)
        except OverflowError:
            break
    return x_c_values, x_n_values

fig, ax = plt.subplots(2, 1, figsize=(15, 10))

color_norm = plt.Normalize(r_values.min(), r_values.max())
scalar_map = cm.ScalarMappable(norm=color_norm, cmap='viridis')

# Calculate the status of the function for every r value
var = 0
for r in r_values:
    var +=1
    print(str(100*var/len(r_values)) + "% Done")
    x_c, x_n = 0.2, 0.8  # Reset initial conditions
    x_c_values, x_n_values = check_convergence(p,delta_c,delta_n,epsilon,x_c,x_n,r)
    ax[0].plot(x_c_values, color=scalar_map.to_rgba(r))
    ax[1].plot(x_n_values, color=scalar_map.to_rgba(r))

# Add colorbar and labels
color_bar = fig.colorbar(scalar_map, ax=ax, orientation='vertical')
color_bar.set_label('r values')

ax[0].set_title('x_c across generations')
ax[0].set_xlabel('Generations')
ax[0].set_ylabel('x_c')

ax[1].set_title('x_n across generations')
ax[1].set_xlabel('Generations')
ax[1].set_ylabel('x_n')

plt.tight_layout()
plt.show()




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