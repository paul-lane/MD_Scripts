###########################################################################
# Inputs
###########################################################################
filelist=("C8mim_working.gro" "C8Fmim_working.gro" "Tf2N_working.gro")  # Create array of filenames
toplist=("C8mim_working.top" "C8Fmim_working.top" "Tf2N_working.top")
nmolecule=("400" "400" "800")			# Create array of nmolecules
outfile="Mixture.gro"				# Name of output file
topout="Mixture2.top"

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
ntop="${#toplist[@]}" 


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

###################################################################
# Create top file
###################################################################

residues=()                            # create an empty array called residues

# For all the individual input files correct the residue names...
	i=0
	while [ $i -lt $ntop ]
	do
		# Run script insert residue name
		./insertResidueNameToTop.sh ${filelist[$i]} ${toplist[$i]} 
		read -r ${residues[i]}<residueName.txt
		i=$((i+1))
	done

	# Copy the first file in the list of topology files to the output file
	cp ${toplist[0]} $topout
	# In specified file remove everything from "Include Position restraint file" t$
	sed -i '/Include Position restraint file/,$d' $topout

	# At this point we have all the info that we need for the first molecule now w$

	# Initialise counter (not we set it to 1 rather than 0) as the file that was a$
	i=1
	#While we are below the number of files...
	while [ $i -lt $ntop ]
	do
        	if [ ${nmolecule[$i]} != 0 ]
        	then
                	# In first file {toplist[$i]} match from line containing "mole$
                	sed -n '/moleculetype ]/,/Include/p' ${toplist[$i]} >> $topout
                	# In file $topout remove everything from "Include Position res$
                	sed -i '/Include Position restraint file/,$d' $topout
	        fi
        	#iterate counter
        	i=$((i+1))
	done

	# Add the end part - copy everything from "Include Postion restraint file" to $
        sed -n '/Include Position restraint file/,$p' ${toplist[0]} >> $topout
	# find number of lines in file count back ntop, replace '1' with nmolecules fo$
        sed -i "$(( $(wc -l < $topout) - 1)),\$s/1/${nmolecule[0]}/g" $topout


	# Add the other residues and the no of atoms
	# these are the last lines of the appropriate files (containing residue name a$
	i=1
	while [ $i -lt $ntop ]
	do
        	if [ ${nmolecule[$i]} != 0 ]
	        then 
                	# Add the last line of each of the input files to the end of t$
                	sed -n '$p' ${toplist[$i]} >> $topout
                	# Replace the "1" in the last line with the number of molecules $
                	sed -i "$(( $(wc -l < $topout) - 1)),\$s/1/${nmolecule[$i]}/g" $topout
        	fi
        	i=$((i+1))
	done

fi


