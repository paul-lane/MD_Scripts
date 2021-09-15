#####################################################################
# Inputs
#####################################################################

rollpts=10			# Number of points for rolling average 
threshold = 10			# Threshold for d_Energy/d_time to determine convergence

#####################################################################
# import numpy as np
import numpy as np
# import pandas as pd
import pandas as pd

# read data file
inputfile = "grepEnergy.txt"
# Unpack as an array
data = np.genfromtxt(inputfile,unpack=True).T

# Time is in the first column (0), energy is in second column (1)
time = data[:,0]
energy = data[:,1]

# Find the derivative dy/dx (d_Energy/d_time)
from numpy import diff
x = time
y = energy
# set up dydx as a pandas series
dydx = pd.Series(diff(y)/diff(x))

# Find the rolling average of the derivative
avedydx = dydx.rolling(rollpts, win_type ='boxcar').mean()
# Pull out the last value of the array, we take the absolute value as this easier to work with
lastvalue=abs(avedydx.tail(1).item())

print(lastvalue)


if lastvalue < threshold:
	print("Converged")
else:
	print("Not converged")

