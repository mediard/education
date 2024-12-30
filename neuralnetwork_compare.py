import numpy as np
import os
import cv2  # E
#from sklearn.model_selection import train_test_split
import random
import os
#import chatGPT
from dnn_app_utils_v3 import load_data
from NNstuff import *
#os.chdir('/Users/mehdiardavan/Documents/lML/projects/NN')

import time
import numpy as np
import h5py
import matplotlib.pyplot as plt
import scipy
from PIL import Image
from scipy import ndimage
#from dnn_app_utils_v3 import *


# Import Same Data as in Andrew Ng's course
train_x_orig, train_y, test_x_orig, test_y, classes = load_data()
print(f"shape of train_x_orig is {train_x_orig.shape}")


# Reshape image tensors into one-dimensional vectors.
# The vector for each image will occupy a column in the matrix.
train_x_flatten = train_x_orig.reshape(train_x_orig.shape[0], -1).T
test_x_flatten = test_x_orig.reshape(test_x_orig.shape[0], -1).T



# Defining the architecture.
unverified_layers = [20, 7, 5 ,  1]
unverified_normalization = "z-score"
unverified_activations = ["ReLU", "ReLU", "ReLU",    "sigmoid"]
unverified_cost = "BCE" #binary cross entropy
unverified_initilizations = [ "random-normal",  "random-normal",   "random-normal",  "random-normal"]
unverified_regularization = 0.00
unverified_reg_type = "L1"
unverified_learning_rates = [ 0.0075, 0.0075,  0.0075, 0.0075,]

# Verifying the architecture.
layers, normalization, activations, cost, initializations, regularization, reg_type, learning_rates = defineArchitecture(unverified_layers, unverified_normalization, unverified_activations, unverified_cost, unverified_initilizations, unverified_regularization, unverified_reg_type ,unverified_learning_rates)

# Skipping my own normalization in order to get similar starting results.
print("starting normalization")
#X_train_n = normalizeInput(X_train, normalization)
X_train_n = train_x_flatten/255.
X_test_n = test_x_flatten/255.
m = X_train_n.shape[1]
n_x = X_train_n.shape[0]


print("initializing parameters")
#W,B = initializeParameters(n_x, layers, initializations)
W,B = initializeParameters_He_M(n_x, layers) # Almost a He Initializer

print("Doing forward Propagation for the first time")
ZZ,AA = forwardPropagation(X_train_n,W,B, activations)

print("Computing cost for the first time")
cost_rep = computeCost(AA[-1],train_y, W,cost, regularization, reg_type)
print("first cost computed")

rep_index=[]
J=[]
rep_index.append(0)
J.append(cost_rep)
print(f"rep no. {0}:\t\tcost:\t {cost_rep}")

repetitions = 3000

print("entering the Training loop")
for rep in range(1, repetitions+1):
    dW, db, dZ = computeGradients(X_train_n,train_y,W,B,AA,ZZ,activations, regularization, reg_type)
    W, B = updatedParameters(dW,db,W,B, learning_rates)
    ZZ,AA = forwardPropagation(X_train_n,W,B, activations)
    cost_rep = computeCost(AA[-1],train_y, W,cost, regularization, reg_type)
    J.append(cost_rep)
    if rep%100==0:
        print(f"rep no. {rep}:\t\tcost:\t {cost_rep}")
y_predictions = (AA[-1]>0.5).astype(int)
accuracy = np.mean(train_y==y_predictions)*100
print(f"Training accuracy is {accuracy} percent")

# Starting testing
ZZ,AA = forwardPropagation(X_test_n,W,B, activations)
t_nx, t_m = X_test_n.shape
Y_hat = AA[-1]
y_predictions = (AA[-1]>0.5).astype(int)
accuracy = np.mean(test_y==y_predictions)*100
print(f"Test accuracy is {accuracy} percent")
