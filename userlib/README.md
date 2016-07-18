## before compiling
mkdir userlib/lib
mkdir userlib/obj

# to recompile the root dictionary, in case something is changed:	
make dictionary	 

#to make the library
make

###################################
Info on existing Executables:
###################################

######################
## digitizer.cpp
# Submit script is submitDigi.py.
# Run in the same way as submitProd~
for i in `seq 0 5`; do python submitDigi.py  -r ${i} -o /afs/cern.ch/work/o/ocolegro/test/ -n 1000 ; done


