# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 21:02:40 2017

@author: geqiao
"""

import numpy as np
import pandas as pd
import scipy.spatial.distance as sdis
import itertools as itools
import easygui
import os

def nchoosek(A, num):
    return np.array(tuple(itools.combinations(A, num)))


# def zeros(a, b):
#     return np.zeros([a, b], dtype=np.int)


# def tril(A):
#     return np.tril(A)


# def ones(a, b):
#     return np.ones([a, b], dtype=np.int)


# def eye(k):
#     return np.identity(k, dtype=np.int)


def quasiOT_sample():
    p = int(input("Number of levels ?\n"))
    r = int(input("Number of sampled random trajectories ?\n"))
    M = int(input("Number of generated random trajectories ?\n"))
    # p= 4
    # r= 50
    # M= 200
    datafile = easygui.fileopenbox(msg="Choose the input file (*.csv)", title="Open File", default="*.csv")

    while datafile is None:
        easygui.msgbox("No file is selected, please try again !")
        datafile = easygui.fileopenbox(msg="Choose the input file (*.csv)", title="Open File", default="*.csv")

    data = pd.read_csv(datafile, sep=',')

    datamatrix = data[['Min', 'Max']].values
    ColName = data['Parameter'].values
    k = len(ColName)
    datamatrix_diff = datamatrix[:, 1] - datamatrix[:, 0]
    B = np.append(np.zeros([1, k]), np.tril(np.ones([k,k])), axis=0)
    delta = p / (2 * (p - 1))
    J_1 = np.ones([k + 1, k])
    J_k = np.ones([k + 1, k])
    T = [0 for i in range(M)]
    i = 1
    A = np.eye(k)

    while i <= M:
        x_Star = np.random.randint(0, p - 1, [1, k]) / (p - 1)
        d = np.ones([1, k])
        d[np.random.rand(1, k) <= 0.5] = -1
        D_Star = np.eye(k)
        D_Star[np.diag_indices(k)] = d
        P_Star = A[np.random.permutation(k), :]
        B_Star = (J_1 * x_Star + (delta / 2) * ((2 * B - J_k).dot(D_Star) + J_k)).dot(P_Star)

        if k > 10:
            for cp in np.nditer(np.arange(p / 2)):
                B_Star[B_Star == (cp - 0.5 * p) / (p - 1)] = (cp + 0.5 * p) / (p - 1)
            for cp in np.nditer(np.arange(p / 2) + p / 2):
                B_Star[B_Star == (cp + 0.5 * p) / (p - 1)] = (cp - 0.5 * p) / (p - 1)

        if np.min(B_Star) >= 0 and np.max(B_Star) <= 1:
            T[i - 1] = B_Star
            i = i + 1

    DIST = np.zeros([M, M])
    for i in range(2, M + 1):
        for j in range(1, i):
            DIST[i - 1, j - 1] = np.sum(sdis.cdist(T[i - 1], T[j - 1], 'euclidean'))
            DIST[j - 1, i - 1] = DIST[i - 1, j - 1]

    vector = np.arange(1, M + 1)
    CombDist = np.zeros([M, 1])
    CombMatrix0 = nchoosek(vector, M - 1)

    for i in range(M, r, -1):
        if i == M:
            r_comb_matrix = nchoosek(vector, 2)
            CombDist_total = 0.
            for index in range(np.size(r_comb_matrix, 0)):
                CombDist_total = CombDist_total + (DIST[r_comb_matrix[index, 1] - 1, r_comb_matrix[index, 0] - 1]) ** 2
            discard_element = np.arange(M, 0, -1)
            for j in range(i):
                CombDist[j] = CombDist_total - np.sum(DIST[CombMatrix0[j, :] - 1, discard_element[j] - 1] ** 2)
        else:
            for j in range(i):
                CombDist[j] = CombDist[j] - np.sum(DIST[CombMatrix0[j, :] - 1, discard_element - 1] ** 2)

        index = np.argmax(CombDist)
        old_vector = vector[:]
        vector = CombMatrix0[index, :]
        discard_element = np.setdiff1d(old_vector, vector)
        CombDist = np.delete(CombDist, index)
        CombMatrix0 = np.delete(CombMatrix0, index, 0)
        CombMatrix0 = CombMatrix0[CombMatrix0 != discard_element]
        CombMatrix0 = np.reshape(CombMatrix0, (i - 1, i - 2))

    best_comb = np.sort(vector)
    TrajectorySet = [T[i] for i in best_comb - 1]
    ParameterSet = TrajectorySet[:]

    datamatrix_diff_transver = np.matlib.repmat(datamatrix_diff, k + 1, 1)
    datamatrix_transver = np.matlib.repmat(np.transpose(datamatrix[:, 0]), k + 1, 1)

    for i in range(r):
        trajectory = TrajectorySet[i]
        ParameterSet[i] = datamatrix_transver + trajectory * datamatrix_diff_transver

    t = pd.DataFrame([], columns=ColName)

    for i in ParameterSet:
        t = t.append(pd.DataFrame(i, columns=ColName))

    t.to_csv(os.path.join(os.path.dirname(datafile), "quasiOT_sample.csv"),
             index=False,
             float_format='%.5f')


if __name__ == '__main__':
    quasiOT_sample()
