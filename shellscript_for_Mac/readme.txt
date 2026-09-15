The shellscript example in this directory, ml3_183610.sh, is for a Mac with a zsh shell. The shellscript 
is put into the user's bin file along with the executable for the Marylie code (in this example the executable is named
ml3_183610.x)

Marylie uses default Fortran file names for input/output.  The default names
have the format fort.N, where N is the logical unit number of a 
Fortran read or write statement.  For input files, the shellscript copies
files with a user-supplied suffix to files with the default name fort.N.

The logical unit numbers for various Marylie I/O read/write commands are specified by the 
Marylie Input File, itself commonly given by a user-supplied name followed by the suffix .mip
Therefore, in the following the Marylie Input File will be called the mip file.

The logical unit number 11 is reserved for the Marylie Input File or mip file. Therefore,
the first command in ml3_18.sh is

cp $1.mip fort.11

In the above command, $1 stands for the user-supplied job name, which will appear in all of the
output file names specified in the shellscript. 

Other user-supplied input-file names and associated logical unit numbers are specified in the mip file.
 
For output files, the shellscript takes files with the default name fort.N and renames them. 
The name assigned by the shellscript is the job name $1 followed by a user-supplied suffix.

Example:

In the example in the directory type6test, the job name is type6 and the mip file type6.mip contains the lines
mapsav   tmo     
16.0000000000000

To run Marylie with the type 6 example, one types in from any directory containing type6.mip and fort.155
ml3_183610.sh type6

-----------------------------------------------------------------------------------------------------------
Note: The menu in this example and the other examples contains elements that are not used; only the 
elements in #lines and #labor are actually used.
-----------------------------------------------------------------------------------------------------------
mapsav is the user-supplied name and tmo is a standard Marylie command name (short for Transfer Map Out).
The 16.00000000000 specifies the logical unit number, which is 16, so the map file (if written) will be fort.16
The map file is actually written if mapsav  appears in the #labor command list, as it does in this example.
Since the user-supplied job name is type6, the shellscript renames fort.16 and gives it the name type6.outmap.

The type 6 magnet is described by a user-supplied shape function (see the coil stack manual, coil_man.pdf). The logical unit number for the input data file file describing the shape function is 155 in this example, so the file fort.155 is read in during initialization of the type 6 magnet.
 
