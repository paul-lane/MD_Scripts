# Script is called with command "./residueName.sh input.gro input.top"
#	where input.gro is the filename of a .gro file containing a single molecule which we want to put multiple into a box for an MD sim
#	where input.top is the filename of a .top file containing of the single molecule
# This code edits the input file by changing the [ moleculetype ] in the .gro file

# filename is the 1st arguement passed when calling the code
grofile=$1
topfile=$2
#echo $grofile

#####################################################################################
# Get residue name
#####################################################################################
#Call the python script to work on the file passed to it when this script was called
python getResidueName.py $grofile

#The python script creates a test file called 'residueName.txt' which contains a list of the residues in the file.
#As we are running this for a single molecule there should only be one residue in there
#Read the file 'residueName.txt' and get whatever is between the single quotation marks - this should be the residue name
residue=$(cat residueName.txt | cut -d"'" -f2 | cut -d"'" -f1)
echo $residue
#Writes residue name to file
echo $residue > residueName.txt

# Overwrite all instances of "Other" in the .top file with the residue name, this should just be in the [ moleculetype ] and [ molecules ] sections
# but care should be taken to ensure that it isn't doing something else that may cause problems! 
sed -i "s/Other/$residue/g" $topfile


