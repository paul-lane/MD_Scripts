The files in this folder are used to populate boxes with a mixture of molecules or residues and are produced from the .gro and .top files of the individual molecules
you specify a lists of files of .gro files and .top files and the number of molecules of each type to be inserted.
populateBox.sh is a script that generates the .gro file of the mixture, it calculates the size of box required and populates it accordingly
populateBox2.sh does the same as populateBox.sh but also generates the .top files based
The other files are scripts of python codes which are called by these main scripts
