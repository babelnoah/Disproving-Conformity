import numpy as np

def create_matrix(k, p, delta_c, delta_n):
    W = np.zeros((2, 2))

    W[0][0] = 1 - k * (1-p) - k * delta_c * np.sqrt(p * (1 - p))
    W[0][1] = k * (1 - p) + k * delta_c * np.sqrt(p * (1 - p))
    W[1][0] = k * p - k * delta_n * np.sqrt(p * (1 - p))
    W[1][1] = 1 - k * p + k * delta_n * np.sqrt(p * (1 - p))

    return W

def compute_lyapunov_exponent(W, iterations=1000):
    # Start with a random vector of norm 1
    v = np.random.rand(2)
    v /= np.linalg.norm(v)
    
    for t in range(1, iterations+1):
        v = np.dot(W, v)
        print(v)
        v /= np.linalg.norm(v)  # renormalize the vector at each step
    print(v)
    lambda_ = np.log(np.linalg.norm(v))
    return lambda_ / iterations

# Setting up parameters
k = 0.5
p = 0.5
delta_c = 0
delta_n = 0.4

W = create_matrix(k, p, delta_c, delta_n)
lyapunov_exp = compute_lyapunov_exponent(W)
print("Maximal Lyapunov Exponent:", lyapunov_exp)
