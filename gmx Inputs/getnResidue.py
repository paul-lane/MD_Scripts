import sys
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC

# filename is the first arguement passed when this python script is called 
filename =  str(sys.argv[1])

#print(filename)


##############################################################
# Set up the .gro file Universe
##############################################################

# Create Universe from .gro file
grofile = mda.Universe(filename)

# Get number of residues in Universe
nresidue = len(grofile.residues)
#print(nresidue)

# write number of residues to file
f=open("residues.txt", "w")
f.write(str(nresidue))
f.close()
