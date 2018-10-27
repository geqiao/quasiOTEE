# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 10:17:51 2017

@author: geqiao
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plotEEresults(miu, std, abs_miu, r, filename):
    f, (ax1, ax2) = plt.subplots(2,1)
    f.suptitle('QuasiOTEE Result')
    ax1.plot(miu, std, 'ro')
    for i, (x, y) in enumerate(zip(miu, std)):
        ax1.annotate('x{}'.format(i+1), (x, y))
    ax1.plot([0, 2 * max(std / np.sqrt(r))], [0, max(std)], 'b-')
    ax1.plot([0, -2 * max(std / np.sqrt(r))], [0, max(std)], 'b-')
    ax1.set_xlabel('$\mu$')
    ax1.set_ylabel('$\sigma$')

    ax2.plot(abs_miu, std, 'ro')
    for i, (x, y) in enumerate(zip(abs_miu, std)):
        ax2.annotate('x{}'.format(i+1), (x, y))
    ax2.set_xlim(0)
    ax2.set_xlabel('$\mu^*$')
    ax2.set_ylabel('$\sigma$')

    fig = plt.gcf()
    fig.set_dpi(100)
    fig.set_size_inches(15, 10)
    fig.savefig(filename, format='pdf')
    return



in_datafile = 'data/quasiOT_sample.csv'
res_datafile = 'data/model_output.csv'
input_range_file = 'data/model_input_range.csv'

input_data_cell = pd.read_csv(in_datafile, sep=',')
input_data = input_data_cell.values
result_data_cell = pd.read_csv(res_datafile, sep=',')
result_data = result_data_cell.values
input_range = pd.read_csv(input_range_file, sep=',')


k = np.size(input_data, 1)
r = np.size(result_data, 0) // (k + 1)
num_outputs = np.size(result_data, 1)

delta = abs(sum(input_data[1:(k + 1), 0] - input_data[0:(k + 1) - 1, 0])) / \
        (input_range.loc[0,'Max']- input_range.loc[0,'Min'])

EE = np.zeros([r, num_outputs, k])
for i in range(r):
    diff_res = result_data[2 + (k + 1) * i - 1:(k + 1) * (i + 1), :] - result_data[(k + 1) * i:(k + 1) * (i + 1) - 1, :]
    diff_inp = (input_data[2 + (k + 1) * i - 1:(k + 1) * (i + 1), :] - input_data[(k + 1) * i:(k + 1) * (i + 1) - 1, :]) != 0
    delta_matrix = np.zeros([k, k])
    delta_matrix[(input_data[2 + (k + 1) * i - 1:(k + 1) * (i + 1), :] - input_data[(k + 1) * i:(k + 1) * (i + 1) - 1, :]) > 0] = delta
    delta_matrix[(input_data[2 + (k + 1) * i - 1:(k + 1) * (i + 1), :] - input_data[(k + 1) * i:(k + 1) * (i + 1) - 1, :]) < 0] = -delta
    for j in range(k):
        EE[i, :, diff_inp[j, :]] = diff_res[j, :] / delta_matrix[j, diff_inp[j, :]]

norm_mean = np.empty([k, num_outputs])
abs_mean = np.empty([k, num_outputs])
stdeviation = np.empty([k, num_outputs])

for i in range(num_outputs):
    norm_mean[:, i] = np.mean(EE[:, i, :], 0)
    abs_mean[:, i] = np.mean(abs(EE[:, i, :]), 0)
    stdeviation[:, i] = np.std(EE[:, i, :], 0)

for i in range(num_outputs):
    filename = 'data/EEResult-{}.pdf'.format(i+1)
    plotEEresults(norm_mean[:,i],
                  stdeviation[:,i],
                  abs_mean[:,i], r, filename)

# ==============================================================================
# input_data_cell = importdata(in_datafile);
#   input_data=input_data_cell.data;
#   result_data_cell=importdata(res_datafile);
#   result_data=result_data_cell.data;
#   
#     r= size(result_data,1)/(k+1);
#   num_traveltimemeasurement=size(result_data,2);
#   EE=zeros(r,num_traveltimemeasurement,k);
#    for i = 1:r
#     diff_res = result_data(2+(k+1)*(i-1):(k+1)*i,:)...
#       -result_data(1+(k+1)*(i-1):(k+1)*i-1,:);
#     diff_inp = (input_data(2+(k+1)*(i-1):(k+1)*i,:)...
#       -input_data(1+(k+1)*(i-1):(k+1)*i-1,:))~=0;
#     delta_matrix=zeros(k,k);
#     delta_matrix((input_data(2+(k+1)*(i-1):(k+1)*i,:)...
#       -input_data(1+(k+1)*(i-1):(k+1)*i-1,:))> 0)= delta;
#     delta_matrix((input_data(2+(k+1)*(i-1):(k+1)*i,:)...
#       -input_data(1+(k+1)*(i-1):(k+1)*i-1,:))<0)=-delta;
#     for j=1:k
#       EE(i,:,diff_inp(j,:))=diff_res(j,:)./delta_matrix(j,diff_inp(j,:));
#     end
#   end
#   
#   % the following codes are used to calculate the mean, absolutemean and
#   % standard deviation
#   
# 
#   norm_mean=(squeeze(mean(EE,1)))';%% mean of EE
#   abs_mean=(squeeze(mean(abs(EE),1)))';%% absolute mean of EE
#   stdeviation=(squeeze(std(EE,1)))';%% standard deviation of EE
# ==============================================================================
