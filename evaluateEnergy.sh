# Script to run gmx energy on a file, grep the columns of numbers and uses energyConverged.py to determine if the energy has flattened off sufficiently
# IN order to do this it reads a .xvg file and greps the two columns of data into a file called 'grepEnergy.txt' energyConverged.py reads this in
# and calculates a rolling average of dE/dt, the last value in this array (i.e. the one at the end of the run is then compared to some convergence
# threshold set in the code, if the convergence criteria is met the the code outputs "converged" otherwise it outputs "not converged". This can be built
# into more complex scripts such as allowing the next stage to be run automatically...

#  Runs  gmx energy with parameter 9 (potential)
echo "9" | gmx energy -f C8mim_Mixture_EM.edr -o energy_test.xvg
# finds everything from the start of the columns to the end of file and  writes it to file 'grepEnergy.txt'
sed -n '/0.000000/,$p' energy_test.xvg > grepEnergy.txt
python energyConverged.py 

