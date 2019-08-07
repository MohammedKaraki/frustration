#! python3.7

import matplotlib.pyplot as plt
import numpy as np

L = input("Enter L: ")

filename = '../data/thermalization-L' + L + '.txt'

with open(filename, 'r') as file:
    lines = []
    for line in file:
        lines.append(line)
    lines = lines[1:-4]
    l = len(lines);
    every = 1
    x = [(float(line.split()[0])) for line in lines[:l-1:every]]
    ene = [(float(line.split()[1])) for line in lines[:l-1:every]]
    mag = [np.abs(float(line.split()[2])) for line in lines[:l-1:every]]
    tim = [(float(line.split()[3])) for line in lines[:l-1:every]]


plt.figure(figsize=(12,7))
plt.subplot(311)
# plt.ylim([-500,0])
plt.plot(x, ene)
plt.subplot(312)
plt.plot(x, mag)
plt.subplot(313)
plt.plot(x, tim)

plt.show()
