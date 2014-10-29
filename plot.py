import numpy as np
import matplotlib.pyplot as pl

k, l, m, t = np.loadtxt("perf.txt").transpose()
pl.plot(k, t, "x", label="current")

k, l, m, t = np.loadtxt("data/tutor.txt").transpose()
pl.plot(k, t, "x", label="tutor")

k, l, m, t = np.loadtxt("data/blas.txt").transpose()
pl.plot(k, t, "x", label="blas")

pl.title("Square Matrix Matrix Multiplication")
pl.xlabel("matrix size")
pl.ylabel("time [s]")

pl.xlim(16, 4096)
pl.loglog(basex = 2)
pl.legend(loc = "upper left")

pl.savefig("current.pdf")
