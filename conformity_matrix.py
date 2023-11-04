import numpy as np

#initial conditions
k = 0.5
p = 0.5
delta_c = 0
delta_n = 0.4
epsilon = 0.5

#delta_c = 0.2

def create_matrix(k,p,delta_c,delta_n):
    W = np.zeros((2, 2))

    W[0][0] = 1 - k * (1-p) - k * delta_c * np.sqrt(p * (1 - p))
    W[0][1] = k * (1 - p) + k * delta_c * np.sqrt(p * (1 - p))
    W[1][0] = k * p - k * delta_n * np.sqrt(p * (1 - p))
    W[1][1] = 1 - k * p + k * delta_n * np.sqrt(p * (1 - p))

    return W

#Matrix
W = create_matrix(k,p,delta_c,delta_n)
print("The matrix is:")
print(W)
print("")

# Eigenvalues and eigenvectors
eigenvalues, eigenvectors = np.linalg.eig(W)
for i in range(len(eigenvalues)):
    print("λ" + str(i) + " = " + str(eigenvalues[i]))
    print("Corresponding eigenvector: " + str(eigenvectors[:, i]))
    print("Difference = " + str(eigenvectors[:, i][1]-eigenvectors[:, i][0]))
print("")

#Inverse
W_inverse = np.linalg.inv(W)
print("The inverse of the matrix is:")
print(W_inverse)
print("")

#Compute (I - W)^-1 (k*delta_c*epsilon, k*delta_n*epsilon)
I = np.eye(2)
vector = np.array([k*delta_c*epsilon, k*delta_n*epsilon])
IW = I - W
IW_inverse = np.linalg.inv(IW)
result = np.dot(IW_inverse, vector)
print("The result is:")
print(result)