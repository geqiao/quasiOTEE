import numpy as np
import pandas as pd


def g_func(X):
    a = np.array([0.02, 0.03, 0.05, 11, 12.5, 13, 34, 35, 37])
    p = 1
    for i in range(len(a)):
        p *= (np.abs(4 * X[:, i] - 2) + a[i]) / (1 + a[i])

    return p


def model(X):
    # user's model which needs to run the SA
    y1 = np.sum(X, axis=1)
    y2 = g_func(X)
    y = np.array([y1, y2]).T
    return y


def save_model_output(y, path):
    df = pd.DataFrame(y, columns=['Output-{}'.format(i + 1) for i in range(y.shape[1])])
    df.to_csv(path, index=False)


def run_model():
    df_X = pd.read_csv('data/quasiOT_sample.csv')
    y = model(df_X.values)
    save_model_output(y, 'data/model_output.csv')


if __name__ == '__main__':
    run_model()
