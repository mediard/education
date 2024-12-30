
import numpy as np
import os
import cv2  # E
from sklearn.model_selection import train_test_split
import random
import os
from dnn_app_utils_v3 import load_data
from NNstuff import *

def defineArchitecture(layers, normalization, activations, cost, initilizations, regularization,reg_type, learning_rates):
    """
    This function ensure the correct values are chosen of the architecture of the neural network.
    layers          :   a list of integers defining the number of neurons in each layer. The length of layers gives L, the number of layers.
    normalization   :   a string defining the normalization of the input values. Currently only z-score is supported.
    activations     :   a list of strings: Provides the activations for each layer. Its lenght is L. Possible values are "sigmoid", "tanh" and "ReLU"
    cost            :   a string defining the type of cost function. The only option for now is "BCE" which is Binay Cross Entropy. In future the following will be added:
                        MSE, MAE, HuberLoss, CCE, HingeLoss,
    regularizations :   a list of small real numbers defining the regularization parameter for each layer.
    learning_rates  :   a list of small real numbers defining the learning rates at each layer.
    """

    assert isinstance(layers, list), "layers must be a list."
    assert all(isinstance(x, int) for x in layers), "All elements in the list must be integers."

    assert normalization == "z-score", "Currently only z-score is supported"

    assert isinstance(activations, list), "activations must be a list."
    assert all(isinstance(a, str) for a in activations), "All elements in activations must be a string"
    allowed_activations = {"sigmoid", "tanh", "ReLU"}
    assert all(a in allowed_activations for a in activations), (
        f"Each element in activations must be one of f{allowed_activations}"
    )

    assert cost == "BCE", "The cost must be equal to 'BCE'."

    assert isinstance(initilizations, list), "initializations must be a list"
    assert all(isinstance(init, str) for init in initilizations), "All initializations be a strings"
    allowed_initializations = {"random-uniform", "random-normal"}
    assert all(init in allowed_initializations for init in initilizations), (
        f"Each element in initializations must be one of f{allowed_initializations}"
    )

    assert isinstance(regularization, (int, float)), "regularizations must be a real number"
    assert regularization>=0, "Regularizations must be zero or positive"
    assert reg_type in ["L1", "L2"], "reg_type must be either 'L1' or 'L2'"

    assert isinstance(learning_rates, list), "learning_rates must be a list"
    assert all(isinstance(lr, (int, float)) for lr in learning_rates ), "All elements in learning_rates must be real numbers"
    assert all(lr>0 for lr in learning_rates), "All elements in learning_rates must be positive"

    return layers, normalization, activations, cost, initilizations, regularization, reg_type, learning_rates





def load_images_and_labels(dataset_path, target_size=(256, 256)):
    categories = ["cats", "dogs", "snakes"]
    X = []
    Y = []

    for label, category in enumerate(categories):
        category_path = os.path.join(dataset_path, category)
        if not os.path.exists(category_path):
            print(f"Warning: Folder {category_path} does not exist.")
            continue

        for file_name in os.listdir(category_path):
            file_path = os.path.join(category_path, file_name)

            # Read the image
            image = cv2.imread(file_path)
            if image is None:
                print(f"Warning: Unable to read {file_path}. Skipping.")
                continue

            # Convert to RGB color space (if needed)
            image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

            # Resize the image to the target size
            image_resized = cv2.resize(image_rgb, target_size)

            # Normalize the image data
            image_normalized = image_resized / 255.0

            # Append to X and Y
            X.append(image_normalized)
            Y.append(1 if category == "cats" else 0)
        print(f"Uploading category {category} finished.")
    X = np.array(X, dtype=np.float32)
    Y = np.array(Y, dtype=np.int32)
    return X, Y




def normalizeInput(X, normalization="z-score"):
    if normalization=="z-score":
        mean = np.mean(X, axis=1, keepdims=True)
        std = np.std(X, axis=1, keepdims=True)+1e-8
        normalizedX = (X-mean)/std
    elif normalization is None:
        normalizedX = X
    else:
        raise ValueError(f"Unsupported normalization: {normalization}")
    return normalizedX

# Initialize parameters
def initializeParameters_He_M(n_x, layers):
    W = []
    B = []
    np.random.seed(1)
    prev_layer_size = n_x
    for layer_size in layers:
        W.append(np.random.randn(layer_size, prev_layer_size) / np.sqrt( prev_layer_size))  # He initialization ... Almost
        B.append(np.zeros((layer_size, 1)))
        prev_layer_size = layer_size
    return W, B


def initializeParameters(n_x, layers, initializations):
    W = []
    B = []
    np.random.seed(1)
    extended_layers = [n_x] + layers
    counter = 0
    for layer in layers:
        counter += 1
        if initializations[counter-1]=="random-uniform":
            w = np.random.rand( extended_layers[counter], extended_layers[counter-1] )/100
            b = np.random.rand( extended_layers[counter],1 )/100
        elif initializations[counter-1]=="random-normal":
            w = np.random.randn( extended_layers[counter], extended_layers[counter-1] )/np.sqrt(extended_layers[counter-1])
            b = np.zeros( extended_layers[counter], 1 )
        W.append(w)
        B.append(b)
    return W,B


#def sigmoid(z):
#    return np.where(
#        z>=0,
#        1 / (1+np.exp(-z)),
#        np.exp(z) / (1+np.exp(z))
#    )

def sigmoid(z):
    z = np.clip(z, -500, 500) # prevents attempting compute exp of large numbers
    return np.where( #exp(z) unstable when z>>1. A trick to avoid exp(z) when z>0.
        z >= 0,
        1 / (1 + np.exp(-z)), # Original formula which is stable for large positive z
        np.exp(z) / (1 + np.exp(z)) # Multiply numerator and denominator by exp(z) then it's stable for large negative z
    )


def sigmoid_p(z):
    temp = sigmoid(z)
    return temp*(1-temp)

def tanh(x):
    return np.tanh(x)

def tanh_p(x):
    a = tanh(x)
    return 1-a*a

def ReLU(x):
    return np.maximum(0,x)

def ReLU_p(x):
    return (x>0).astype(float)

def forwardPropagation(X,W,B, activations):
    AA=[]
    ZZ=[]
    L = len(W)

    for l in range(L):
        if l==0:
            A_prev = X
        else:
            A_prev = AA[l-1]
        Z = np.dot(W[l],A_prev) + B[l]
        command = activations[l]+"(Z)"
        A = eval(command)
        AA.append(A)
        ZZ.append(Z)
    return ZZ,AA

def computeCost(A,Y, W,cost, regularization, reg_type):
    m = Y.shape[1]
    if cost == 'BCE':
        epsilon = 1e-8
        y_predictions = np.clip(A.flatten(), epsilon, 1-epsilon)
        y_originals = Y.flatten()
        bce_loss = -np.mean(y_originals * np.log(y_predictions) + (1 - y_originals) * np.log(1 - y_predictions))
        penalty = 0
        for w in W:
            wf = w.flatten()
            if reg_type == "L1":
                penalty += regularization * np.sum(np.abs(wf))/m
            elif reg_type == "L2":
                penalty += regularization * np.sum(np.square(wf))/(2*m)
            else:
                raise ValueError("Error: Invalid regularization type!")
        total_loss = bce_loss + penalty
    else:
        raise ValueError("Error: Invalid cost type provided!")

    return total_loss


def computeGradients(X,Y,W,B,AA,ZZ,activations, regularization, reg_type):
    # returns dW, dB, dZ all list of lenth L
    # start from layer L
    dW = []
    db = []
    dZ = []
    m = X.shape[1]
    L = len(AA)
    dZ_L = AA[-1]-Y
    if len(AA)<2:
        A_prev = X
    else:
        A_prev = AA[-2]
    dw_L = 1.0/m * np.dot(A_prev, dZ_L.T).T
    if reg_type=="L1":
        dw_L = dw_L + np.sign(W[-1]) * regularization / m
    elif reg_type=="L2":
        dw_L = dw_L + (W[-1]) * regularization / m
    db_L = 1.0/m * np.sum(dZ_L, axis=1, keepdims=True)
    # Add dZ_L, dw_L and db_L to the pool of data
    dW.append(dw_L)
    db.append(db_L)
    dZ.append(dZ_L)
    dzl = dZ_L

    for l in range(L - 1, 0, -1): # l will indicate the actual layer number not the index in the list AA or W or B. For index use l-1
        dzl = np.dot(W[l+1-1].T, dzl) * eval(activations[l-1]+'_p(ZZ[l-1])')
        dbl = 1/m*np.sum(dzl, axis=1, keepdims=True)
        dwl = 1/m* np.dot(dzl, AA[l-1-1].T if l>1 else X.T)
        if reg_type=="L1":
            dwl = dwl + regularization * np.sign(W[l-1]) / m
        else:
            dwl = dwl + regularization * (W[l-1]) / m
        dZ.insert(0,dzl) #= dzl + dZ
        db.insert(0,dbl) # = dbl + db
        dW.insert(0,dwl) #= dwl + dW

    return dW, db, dZ


def updatedParameters(dW,db,W,B, learning_rates):

    L = len(learning_rates)
    for l in range(1-1,L+1-1):
        w = W[l]
        assert w.shape==dW[l].shape, "w and dW must be of the same shape"


        w = w - learning_rates[l] * dW[l]
        W[l] = w
        b = B[l]
        assert b.shape==db[l].shape, "w and dW must be of the same shape"
        b = b - learning_rates[l] * db[l]
        B[l] = b

    return W, B
