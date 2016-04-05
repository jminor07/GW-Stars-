import matplotlib.pyplot as plot
import numpy as np
import math
from pylab import *

wave,flux=np.loadtxt('t05000_g+3.5_m01p04.spec',skiprows=7,usecols=(0,1),unpack=True)

plot(wave,flux)
show()
