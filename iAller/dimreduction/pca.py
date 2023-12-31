import numpy as np
from sklearn.decomposition import PCA

def pca(encodings, n_components = 2):
    encodings = np.array(encodings)[1:]
    data = encodings[:, 1:]
    shape = data.shape
    data = np.reshape(data, shape[0] * shape[1])
    data = np.reshape([float(i) for i in data], shape)
    newData = PCA(n_components = n_components).fit_transform(data)
    pca = []
    for i in range(len(data)):
        pca.append([encodings[i][0]] + list(newData[i]))
    return pca

