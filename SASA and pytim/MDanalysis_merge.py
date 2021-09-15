import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC

# Create Universe from .gro file
grofile = mda.Universe("pdb200.gro")	#"C8mim25_s_17.gro"

# Get number of atoms in Universe
natoms = len(grofile.atoms)
print(natoms)


import numpy as np

# setup reading the SASA file
xvgfile = "pdb200_A.xvg"		#"C8mim25_s_17A_22.xvg"
sasa = np.genfromtxt(xvgfile,unpack=True).T

# setup reading the pytim file
#pytimfile="layers.pdb"
#pytim = np.genfromtxt(pytimfile,unpack=True).T



# Atom selection
all =  grofile.atoms					# Create a selection of all atoms
#all.write('all.gro')					# Write a .gro file (should be the same as the input file - this is just a check)

# Convert position arrays
pos = all.positions/10.0				# Positions are a factor 10 different to input
vel = all.velocities/10.0				# ... as are speeds, so we correct both

# Get values from grofile
residueNo = all.resnums
residueName = all.resnames
atomName = all.names

# Seperate positions and velocities into x, y, z components
xpos = pos[:,0]
ypos = pos[:,1]
zpos = pos[:,2]

xvel = vel[:,0]
yvel = vel[:,1]
zvel = vel[:,2]

# SASA output create named arrays
atomNo = sasa[:,0]
SASAarea = sasa[:,1]
SASAerror = sasa[:,2]

# At this point I need all the info from the grofile as variables stored above as I am creating a new Universe

pdbfile = mda.Universe("layers.pdb")    #"C8mim25_s_17.gro"

allpytim =  pdbfile.atoms


pytimsurf=allpytim.tempfactors
pos=allpytim.positions/10.0
xpos2 = pos[:,0]
ypos2 = pos[:,1]
zpos2 = pos[:,2]

residueNo2 = allpytim.resnums
residueName2 = allpytim.resnames
atomName2 = allpytim.names


SASAcount = 0						# Initialise SASA count
pytimCount = 0						# INitialise pytim count


# create two output files
outfile1 = open("output1.txt", 'w')
outfile2 = open("output2.txt", 'w')

# Write Output file 1 - This contains all the information about atoms that SASA has identified as being at the surface
for i in range(natoms):					# For all atoms in selection...
	if (SASAarea[i] != 0.0):			# If they have a non-zero SASA area write out the gro file and SASA information

# Write Output file 1 - This contains all the information about atoms that SASA has identified as being at the surface
		outfile1.write(str(i+1) +"\t" +		# Write to file:	# atom no (derived from counter)
		str(residueNo[i]) + "\t" + 					# residue number
		str(residueName[i])+ "\t" + 					# residue name
		str(atomName[i]) + "\t" + 					# atom name
		str('{:.3f}'.format(xpos[i])) + "\t" +				# x position (formatted as in .gro file)
		str('{:.3f}'.format(ypos[i])) + "\t" +                          # y position (formatted as in .gro file)
		str('{:.3f}'.format(zpos[i])) + "\t" +                          # z position (formatted as in .gro file)
		str(int(atomNo[i])) + "\t" +					# atom number (from SASA file)
		str(SASAarea[i]) + "\t" +					# SASA area
                str(residueNo2[i]) + "\t" +					# residue number (from pytim)
		str(residueName2[i]) + "\t" +					# residue name (from pytim)
		str(atomName2[i]) + "\t" +					# atom name (from pytim)
		str('{:.3f}'.format(xpos2[i])) + "\t" +                         # x position (from pytim)
                str('{:.3f}'.format(ypos2[i])) + "\t" +                         # y position (from pytim)
                str('{:.3f}'.format(zpos2[i])) + "\t" +                         # z position (from pytim)
		str(pytimsurf[i]) + "\n")					# Pytim surface atom identifier

		SASAcount = SASAcount + 1		# Count the number of SASA surface atoms

	if (pytimsurf[i] == 1.0):			# If they have been identified by pytim as being at the surface

# Write Output file 2 - This contains all the information about atoms that pytim has identified as being at the surface
		outfile2.write(str(i+1) +"\t" +         # Write to file:        # atom no (derived from counter)
		str(residueNo[i]) + "\t" +                                      # residue number
		str(residueName[i])+ "\t" +                                     # residue name
		str(atomName[i]) + "\t" +                                       # atom name
		str('{:.3f}'.format(xpos[i])) + "\t" +                          # x position (formatted as in .gro file)
		str('{:.3f}'.format(ypos[i])) + "\t" +                          # y position (formatted as in .gro file)
		str('{:.3f}'.format(zpos[i])) + "\t" +                          # z position (formatted as in .gro file)
		str(int(atomNo[i])) + "\t" +                                    # atom number (from SASA file)
		str(SASAarea[i]) + "\t" +                                       # SASA area
                str(residueNo2[i]) + "\t" +                                     # residue number (from pytim)
                str(residueName2[i]) + "\t" +                                   # residue name (from pytim)
                str(atomName2[i]) + "\t" +                                      # atom name (from pytim)
		str('{:.3f}'.format(xpos2[i])) + "\t" +                         # x position (from pytim)
		str('{:.3f}'.format(ypos2[i])) + "\t" +                         # y position (from pytim)
		str('{:.3f}'.format(zpos2[i])) + "\t" +                         # z position (from pytim)
		str(pytimsurf[i]) + "\n")                                       # Pytim surface atom identifier

		pytimCount = pytimCount + 1

outfile1.close()
outfile2.close()
print(SASAcount)
print(pytimCount)
