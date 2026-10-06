import sys
import tensorflow as tf
import os.path
import numpy as np
import scipy.io
from scipy.linalg import svd
import time
from tensorflow import keras

random_seed = 1234
np.random.seed(random_seed)
tf.random.set_seed(random_seed)
import random
random.seed(random_seed)
import math
    
dataFile = 'Snapshots.mat'
data = scipy.io.loadmat(dataFile)
U = data['U']
Uddot = data['Uddot']

# Constants
T = 5

# Assuming U is already defined as a 2D numpy array
N_T = U.shape[1]  # Number of columns in U

# Time vector
t = np.linspace(0, T, N_T)
t = tf.cast(t, tf.float32)

dt = t[1] - t[0]

dt_value = dt.numpy()

dt_value = math.floor(dt_value * 1000) / 1000

nSkip = int(np.round(0.1/dt));

print('nSkip is ' + str(nSkip))

# Load scaling values
dataScaling = scipy.io.loadmat("ScalingVals.mat")

s_U = dataScaling['s_U']
s_Uddot = dataScaling['s_Uddot']

s_U = tf.convert_to_tensor(s_U, dtype=tf.float32)
s_Uddot = tf.convert_to_tensor(s_Uddot, dtype=tf.float32)

# Skip first new snapshots in which the excitiation is taking place without the dynamics being started
U = U[:, nSkip:]
Uddot = Uddot[:, nSkip:]

# Normalize matrices
U = U / s_U
Uddot = Uddot / s_Uddot

r = 100

Basis = scipy.io.loadmat("Phi_POD.mat")
Utilde = Basis['Phi']
Phi = Utilde[:, :r]  # Reduced basis

scipy.io.savemat("Phi" + str(r) + ".mat", {"Phi": Phi})

out_reduced = np.float32((Phi.T @ np.float32(U)).T)
out_ddot_reduced = np.float32((Phi.T @ np.float32(Uddot)).T)
NoutputTotal = out_reduced.shape[1]

# Initialize trainable variables for M and K
initializer = tf.keras.initializers.GlorotUniform()
L = tf.Variable(initializer((r, r)), dtype=tf.float32, trainable=True)
W = tf.Variable(initializer((r, r)), dtype=tf.float32, trainable=True)
def get_M(): return tf.matmul(L, tf.transpose(L))
def get_K(): return tf.matmul(W, tf.transpose(W))

# Placeholder second derivative of U (assuming precomputed)
d2Udt2 = out_ddot_reduced  # Replace with actual data
U = out_reduced  # Replace with actual data

# Convert NumPy arrays to TensorFlow tensors
d2Udt2 = tf.convert_to_tensor(d2Udt2, dtype=tf.float32)  # (Nt, r)
U = tf.convert_to_tensor(U, dtype=tf.float32)            # (Nt, r)

lambda_1 = 1e-7

# Compute residual with correct transpositions
def compute_residual():
    return tf.matmul(get_M(), tf.transpose(d2Udt2)) + tf.matmul(get_K(), tf.transpose(U))


# Define loss function with regularization term
@tf.function
def loss_fn():
    residual = compute_residual()
    # L2 loss term
    l2_reg = lambda_1 * (tf.reduce_sum(tf.square(get_M())) + tf.reduce_sum(tf.square(get_K()))) 
    return tf.reduce_sum(tf.square(residual)) + l2_reg  # Mean squared residual + regularization


lr_schedule = keras.optimizers.schedules.ExponentialDecay(
    initial_learning_rate=1e-3,
    decay_steps=12000,
    decay_rate=0.97)

# Optimizer
optimizer = tf.keras.optimizers.Adam(learning_rate=lr_schedule)

# Training step
@tf.function
def train_step():
    with tf.GradientTape() as tape:
        loss = loss_fn()
    gradients = tape.gradient(loss, [L, W])
    optimizer.apply_gradients(zip(gradients, [L, W]))
    return loss
    
if os.path.isdir('OpInf_Matrices/') is False:
    os.mkdir('OpInf_Matrices/')

# Training loop
num_epochs = 200000
save_freq = 50000
sys.stdout.flush()
startTime = time.time()

tol = 1e-3
patience = 1000

M_old = get_M().numpy()
K_old = get_K().numpy()

counter = 0


for epoch in range(num_epochs):

    loss_value = train_step()

    if (epoch+1) % 100 == 0:

        M_new = get_M().numpy()
        K_new = get_K().numpy()

        change_M = (
            np.linalg.norm(M_new-M_old)
            /
            np.linalg.norm(M_old)
        )

        change_K = (
            np.linalg.norm(K_new-K_old)
            /
            np.linalg.norm(K_old)
        )

        operator_change = max(
            change_M,
            change_K
        )

        print(
            f"Epoch {epoch+1}, "
            f"loss={loss_value.numpy():.4e}, "
            f"operator change={operator_change:.3e}"
        )


        if operator_change < tol:
            counter += 100
        else:
            counter = 0


        if counter >= patience:
            print("Operators converged")
            break


        M_old = M_new
        K_old = K_new
    
        
    if (epoch + 1) % save_freq == 0:
        M_numpy, K_numpy = get_M().numpy(), get_K().numpy()
        scipy.io.savemat(f"OpInf_Matrices/learned_matrices_dt{dt_value}_r{r}_{lambda_1}_epochs{epoch+1}.mat", {"M": M_numpy, "K": K_numpy})

endTime = time.time()
elapsed = (endTime - startTime) / 60.0
print("The training took {:.4} minutes".format(elapsed))

residual = compute_residual()
loss = loss_fn()


print('Final Loss Value is ' + str(loss))

# Convert M and K to NumPy arrays
M_numpy = get_M().numpy()
K_numpy = get_K().numpy()

scipy.io.savemat(f"OpInf_Matrices/learned_matrices_dt{dt_value}_r{r}_{lambda_1}.mat", {"M": M_numpy, "K": K_numpy})
np.save(f"OpInf_Matrices/M_learned_dt{dt_value}_r{r}_{lambda_1}.npy", M_numpy)
np.save(f"OpInf_Matrices/K_learned_dt{dt_value}_r{r}_{lambda_1}.npy", K_numpy)
print("Saved M and K matrices.")

