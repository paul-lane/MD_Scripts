import sys
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC

filename =  str(sys.argv[1])

#print(filename)


##############################################################
# Set up the .gro file Universe
##############################################################

# Create Universe from .gro file
grofile = mda.Universe(filename)
# Create a selection containing all residues (there is only one)
selection=grofile.residues
# Create a variable to hold the name of the residue
residueName = selection.resnames
print(residueName)

# write number of residues to file
f=open("residueName.txt", "w")
f.write(str(residueName))
f.close()
