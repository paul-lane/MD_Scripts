###########################################################################
# Inputs
###########################################################################
filelist=("C8mim.gro" "C8Fmim.gro" "Tf2N.gro")  # Create array of filenames
nmolecule=("800" "0" "800")			# Create array of nmolecules
outfile="Mixture.gro"				# Name of output file

#final parameters
try=3000 					# try is the no of tries to position molecule in box

# test parameters
testtry=100					# test try is the number of tries when calculating the volume
testSide=4.0					# testSide is the initial guess at box side

###########################################################################
# start of code
###########################################################################

nfiles="${#filelist[@]}"                        # Get size of file array
nnmol="${#nmolecule[@]}"			# Get size of nmolecules array

if(( $nfiles != $nnmol ))			# If array sizes are different
then 						# we have an error so stop
	echo "ERROR: Number of file and number of molecule arrays mismatched" 
else						# otherwise continue
	echo "continue with code"

	i=0
	completed=()				# create an empty array called completes
	volume=()				# create an empty array for volumes
	totalnmol=()				# create an empty array for total no of molecules
	while [ $i -lt $nfiles ]		# for every item in array...
	do
	volume[$i]=0				# set initial volume value to zero
	if [ ${nmolecule[$i]} -gt 0 ]
	then
		completed[$i]="false"		# initialize array completed to false
	else					# for all files with a nmolecules greater than zero
		completed[$i]="true"		# (we want to skip any with no molecules in as they don't contribute)
	fi
	totalnmol=$(($totalnmol+${nmolecule[$i]}))
#	totalnmol=$(bc -l <<< "${nmolecule[$i]}+$totalnmol")  
	echo ${nmolecule[$i]}
	echo ${completed[$i]}
	i=$((i+1))				# iterate counter
	done


	echo "no of files in array " $nfiles

	i=0      				# Initialise counter
	while [ $i -lt $nfiles ]		# While we are below the number of files...
	do
		while [ ${completed[$i]} == "false" ]
		do
			gmx insert-molecules -box $testSide $testSide $testSide -ci ${filelist[$i]} -nmol ${nmolecule[$i]} -try $testtry -rot xyz -o test.gro
			# We see how many mol1 molecules are in the box using MDAnalysis
			python getnResidue.py test.gro
			#Read the file 'residues.txt' this is produced by checkmResidue.py to output the number of residues in the output file
			read -r nres<residues.txt
			# Factor is the value you would need to scale the starting volume by to get all the molecules in
			factor=$(bc -l <<< "(1.0*${nmolecule[$i]}/$nres)+0.05")
			# We are only interested in the factor to 1 decimal place (the 0.05 above means we round up)
			# fact is the rounded version of factor
			printf -v fact "%.1f \n" $factor
			#The volume of each component is dependent on the factor multiplied by the test volume
			volume[$i]=$(bc -l <<< $testSide*$testSide*$testSide*$fact)
			echo ${volume[$i]}
			# remove temporary files
			rm residues.txt
			rm test.gro
			# Change completed status to indicate we are complete
			completed[$i]="true"
		done
        	i=$((i+1))                 	 # Iterate counter
	done

	i=0
	totalvol=0
	while [ $i -lt $nfiles ]
	do
		echo $i
		totalvol=$(bc -l <<< "${volume[$i]}+$totalvol")		# Calculate total volume
		i=$((i+1))                       # Iterate counter


	done

	echo "total volume" $totalvol
	sidelen=$(bc -l <<< "e(l($totalvol)/3) +0.05")			# Calculate length of box side 0.05 added for rounding later 
	echo $sidelen
        printf -v side "%.1f \n" $sidelen				# Convert side length to a 1dp value
	echo $side
	echo $totalnmol
	echo "starting production run"

############################################################################################################
# Actual production run
############################################################################################################
	filegenerated=false
	echo $filegenerated
	echo $nfiles
#	nsteps=$nfiles-1

	while [ $filegenerated == false ]
	do

        	i=0                                     # Initialise counter
        	while [ $i -lt $nfiles ]                # While we are below the number of files...
        	do

			if [ $i == 0 ]
			then
				# Gromacs command to add molecule for first iteration add to an empty box
				echo "Creating box" $i
        			gmx insert-molecules -box $side $side $side -ci ${filelist[$i]} -nmol ${nmolecule[$i]} -try $try -rot xyz -o $outfile
			else
		        	# For subsequent iterations add to the box created in first  iteration
                        	echo "adding to box" $i
        			gmx insert-molecules -f $outfile -ci ${filelist[$i]} -nmol ${nmolecule[$i]} -try $try -rot xyz -o $outfile
			fi
	        	i=$((i+1))                       # Iterate counter

		done

		# Call getnResidue to find out how many residues have been added to the output file
		python getnResidue.py $outfile
		# Read the file 'residues.txt' this is produced by checkmResidue.py to output the number of residues in the output file
		read -r firstline<residues.txt
		echo $firstline
		# if less molecules than desired rerun with bigger box
		if(($firstline == $totalnmol))
		then
        		echo "Input file generated"
			filegenerated=true
		# otherwise redo with box side length increase
		else
        		echo "Box is too small rerunning with bigger box"
			side=$(bc -l <<< $side+0.1)
			echo $side "New box size"
		fi

	done

fi


