#!/usr/bin/env python3
from string import *
#import os, commands, getopt, sys, exceptions
import numpy as np
import matplotlib.pyplot as plt

# read time, convert to days since start
print("reading file...")
data = np.loadtxt("tbot.out",usecols=(1,))

# statefreq=24, timestep=8.33s
delta = 24*8.3333333333333333333/(24*3600.)
t = np.arange(0, delta*data.size, delta)


print("constructing plot...")
fig, axs = plt.subplots()
axs.plot(t,data,label='min tbot',color='b')
axs.grid(True)
axs.set(xlabel='days', ylabel='degrees K',
       title='TBOT min')

#plt.show()
print("writing plot...")
plt.savefig("temp.png")
