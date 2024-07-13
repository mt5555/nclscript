#!/usr/bin/env python3
from string import *
#import os, commands, getopt, sys, exceptions
import numpy as np
import matplotlib.pyplot as plt
#
# simple script to plot x-y curve from text file
#
#
# read in data (y axis data):
#filename = "swtc2-ne30.linf.errors"
filename = "swtc2-CA.linf.errors"
title = "SWTC2"
ylabel = "max error"

#filename = "TBOT.out"
#title = "min tbot"
#ylabel = "degrees K'




print("reading file...")
data = np.loadtxt(filename,usecols=(1,))

t  = np.loadtxt(filename,usecols=(0,))
# or compute time based on knowledge of timestep between samples:
# statefreq=24, timestep=8.33s
#delta = 24*8.3333333333333333333/(24*3600.)
#t = np.arange(0, delta*data.size, delta)


print("constructing plot...")
fig, axs = plt.subplots()
axs.plot(t,data,label=title,color='b')
axs.grid(True)
axs.set(xlabel='days', ylabel=ylabel,title=title)

#plt.show()
print("writing plot...")
plt.savefig("temp.png")
