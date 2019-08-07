#! python3.7

import matplotlib.pyplot as plt
import numpy as np

L = input("Enter L: ")

filename = '../data/collecting-L' + L + '.hx.txt'


with open(filename, 'r') as file:
    lines = []
    for line in file:
        lines.append(line)
    lines = lines[1:]
    l = len(lines);
    every = 1 #int(np.ceil(l / 100))
    field = []
    ene = []
    ene_err = []
    mag = []
    binder = []
    for i in range(l-1):
        line = lines[i]
        sp = line.split();
        field.append(float(sp[0]))
        ene.append(float(sp[1]))
        ene_err.append(float(sp[2])*100)
        mag.append(float(sp[3]))
        binder.append(float(sp[5]))

print (ene)
print (ene_err)

plt.figure(figsize=(12,7))
plt.subplot(211)
plt.errorbar(field, ene, yerr=ene_err)
plt.ylabel("Energy")
plt.title("Err Bars 1000x smaller")
plt.subplot(212)
plt.plot(field, binder)
plt.ylabel("Energy Binder Cumulant")
plt.xlabel("Field")


plt.show()
