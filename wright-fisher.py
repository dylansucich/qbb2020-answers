#!/Users/dylansucich/miniconda3/bin/python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import seaborn as sns
import random
import matplotlib as mpl
import seaborn as sns 
from numpy.random import binomial

# def simulation(p, N, s, repetitions = 1000):
#     N = int(2 * N) # why wouldn't this be an integer by default? 
#     n1 = np.ones(repetitions) * (N * p)
#     T = np.empty_like(n1)
#     update = (n1 > 0) & (n1 < N)
#     t = 0
#     while update.any():
#         t += 1
#         p = n1 * (1 + s) / (N + n1 * s)
#         print(p)
#         n1[update] = np.random.binomial(N, p[update])
#         T[update] = t
#         update = (n1 > 0) & (n1 < N)
#     return n1 == N, T
# #print(simulation(N=1000, s = 0.1)) 
# fixations, times = simulation(p=0.5, N=100, s=0, repetitions=1000)
# fixation_prob = fixations.mean()
# fixation_time = times[fixations].mean()
# w, h = plt.rcParams['figure.figsize'] # customize window size 
# #fig, ax = plt.subplots(1, 2, figsize = (2 * w, h))
# fig, ax = plt.subplots()
# sns.distplot(times[fixations], ax = ax)
# ax.axvline(times[fixations].mean(), color = 'k', ls = '--')
# ax.set(xlabel='Fixation time', ylabel = 'Frequency')
# ax.set_title("Time to Fixation")
# fig.savefig("Histogram_part1")

# p=0.5
# repetitions = 1000
# s = 0
# Nrange = np.logspace(2, 3, dtype=np.uint64)

# def fix_time_simulation(N):
#     fixations, times = simulation(p=p, N=N, s=s, repetitions=repetitions)
#     fixation_time_mean = times[fixations].mean()
#     fixation_time_std =  times[fixations].std(ddof=1) / np.sqrt(repetitions)
#     return fixation_time_mean, fixation_time_std

# fix_time_sim = np.array([
#     fix_time_simulation(N=N)
#     for N in Nrange
# ])

# def fixation_time_plot(N, mean, sem):
#     fig, ax = plt.subplots(1, 1)
#     ax.errorbar(x=N, y=mean, yerr=sem, 
#                 fmt='o', capsize=5, label='Simulation')
#     ax.set(
#         xlabel='Population size (N)',
#         ylabel='Fixation time',
#         xscale='log', 
#         xlim=(0.5 * Nrange.min(), 1.5 * Nrange.max()),
#     )
#     return fig, ax

# fixation_time_plot(Nrange, fix_time_sim[:,0], fix_time_sim[:,1]);
# plt.show()
# fig.savefig("Histogram_part2.png")


def simulation_3(p, N, s, repetitions = 1000):
    N = int(2 * N) # why wouldn't this be an integer by default? 
    n1 = np.ones(repetitions) * (N * p)
    T = np.empty_like(n1)
    update = (n1 > 0) & (n1 < N)
    t = 0
    while update.any():
        t += 1
        p = n1 * (1 + s) / (N + n1 * s)
        n1[update] = np.random.binomial(N, p[update])
        T[update] = t
        update = (n1 > 0) & (n1 < N)
    return n1 == N, T
#print(simulation(N=1000, s = 0.1)) 
#fixations, times = simulation_3(p=allele_freq, N=100, s=0, repetitions=10000)
allele_freqs = [0.2, 0.4, 0.6, 0.8]
fixation_times3 = {}
for allele_freq in allele_freqs:
    fixation_times3[allele_freq] = []
    for i in range(100):
        fixation_times3[allele_freq].append(simulation_3(p=allele_freq, N=100, s=0, repetitions=1000))
mean3 = []
std_dev3 = []
for allele_freq in allele_freqs:
    mean3.append(np.mean(fixation_times3[allele_freq]))
    std_dev3.append(np.std(fixation_times3[allele_freq]))
fig,ax = plt.subplots()
ax.bar([x for x in range(1,len(allele_freqs)+1)], mean3, yerr= std_dev3)
ax.set_xlabel("Starting Allele Frequency")
ax.set_xticks([x for x in range(1,len(allele_freqs)+1)])
ax.set_xticklabels(["0.2", "0.4", "0.6", "0.8"])
ax.set_ylabel("Generations to Fixation")
plt.title("Varying Allele Frequency")
fig.savefig("wright-fisherPlot3.png")
plt.show()

def simulation_4(p, N, s, repetitions = 100):
    N = int(2 * N) # why wouldn't this be an integer by default? 
    n1 = np.ones(repetitions) * (N * p)
    T = np.empty_like(n1)
    update = (n1 > 0) & (n1 < N)
    t = 0
    while update.any():
        t += 1
        p = n1 * (1 + s) / (N + n1 * s)
        n1[update] = np.random.binomial(N, p[update])
        T[update] = t
        update = (n1 > 0) & (n1 < N)
    return n1 == N, T
#print(simulation(N=1000, s = 0.1)) 
#fixations, times = simulation_3(p=allele_freq, N=100, s=0, repetitions=10000)
sel_coeffs = [0.2, 0.4, 0.6, 0.8]
fixation_times4 = {}
for sel_coeff in sel_coeffs:
    fixation_times4[sel_coeff] = []
    for i in range(100):
        fixation_times4[sel_coeff].append(simulation_4(p=0.5, N=100, s=sel_coeff, repetitions=1000))
mean3 = []
std_dev3 = []
for sel_coeff in sel_coeffs:
    mean3.append(np.mean(fixation_times4[sel_coeff]))
    std_dev3.append(np.std(fixation_times4[sel_coeff]))
fig,ax = plt.subplots()
ax.bar([x for x in range(1,len(sel_coeffs)+1)], mean3, yerr= std_dev3)
ax.set_xlabel("Selection Coefficient")
ax.set_xticks([x for x in range(1,len(sel_coeffs)+1)])
ax.set_xticklabels(["0.2", "0.4", "0.6", "0.8"])
ax.set_ylabel("Generations to Fixation")
plt.title("Varying Selection Coefficient")
fig.savefig("wright-fisherPlot4.png")
plt.show()    
# def simulation(p, N, s, repetitions=1000):
#     N = int(2*N) # very important! if this is a float than the while loop will be ~infinite
#     n1 = np.ones(repetitions) * (N*p)
#     T = np.empty_like(n1)
#     update = (n1 > 0) & (n1 < N)
#     t = 0
    
#     while update.any():
#         t += 1
#         p = n1 * (1 + s) / (N + n1 * s) 
#         n1[update] = np.random.binomial(N, p[update])
#         T[update] = t
#         update = (n1 > 0) & (n1 < N)

#     return n1 == N, T


# print("Fixate %", simulation(p=0.5, N=100, s=0)[0].mean())


# fixations, times = simulation(p=0.5, N=100, s=0, repetitions=1000)
# fixation_prob = fixations.mean()
# fixation_time = times[fixations].mean()


# w, h = plt.rcParams['figure.figsize']
# # fig, ax = plt(1, 2, figsize=(2 * w, h))
# fig, ax = plt.subplots()

# sns.distplot(times[fixations], ax=ax)
# ax.set_title( "Time to fixation" )
# ax.axvline(times[fixations].mean(), color='k', ls='--')
# ax.set(xlabel='Fixation time', ylabel='Frequency')

# fig.savefig("wright-fisherPlot1.png")


# repetitions = 1000
# s = 0
# Nrange = np.logspace(1, 6, 10, dtype=np.uint64)

# def fix_time_simulation(N):
#     fixations, times = simulation(p=0.5, N=N, s=s, repetitions=repetitions)
#     fixation_time_mean = times[fixations].mean()
#     fixation_time_std =  times[fixations].std(ddof=1) / np.sqrt(repetitions)
#     return fixation_time_mean, fixation_time_std

# fix_time_sim = np.array([
#     fix_time_simulation(N=N)
#     for N in Nrange
# ])

# def fixation_time_plot(N, mean, sem):
#     fig, ax = plt.subplots(1, 1)

#     ax.errorbar(x=N, y=mean, yerr=sem, 
#                 fmt='o', capsize=5, label='Simulation')

#     ax.set(
#         xlabel='Population size (N)',
#         ylabel='Fixation time',
#         xscale='log', 
#         xlim=(0.5 * Nrange.min(), 1.5 * Nrange.max()),
#     )
#     return fig, ax

# fixation_time_plot(Nrange, fix_time_sim[:,0], fix_time_sim[:,1]);










