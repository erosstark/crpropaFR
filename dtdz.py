from crpropa import *
import numpy as np

z = np.linspace(0, 2.5, 100)
dzdt = [1/((1+i)*hubbleRate(i)) for i in z]

with open("dtdzL.txt", "w") as f:
    for i in range(len(dzdt)):
        f.write(str(dzdt[i]) + " " + str(z[i]) + "\n")