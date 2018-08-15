# COMPILER
f90comp = gfortran

# FLAGS
p_flags = -fopenmp
flags = -O3 -g
###-fbounds-check  -fdefault-real-8

# THINGS
%.o : %.f90
	$(f90comp) $(flags) -c -o $@ $<
%.o : %.f
	$(f90comp) $(flags) -c -o $@ $<

# MIGRATION
mk3d : mk3d.o
	$(f90comp) $(flags) -o mk3d mk3d.o

# DATA GENERATION
gendata : GenData.o
	$(f90comp) $(flags) -o gendata GenData.o

# SLICE
slice : Slice.o
	$(f90comp) $(flags) -o slice Slice.o

# 2D INTERFACES
mint : make_interfaces.o idefs.o
	$(f90comp) -o mint make_interfaces.o idefs.o

# 2D VELOCITIES
mvel : make_vgrids.o vdefs.o
	$(f90comp) -o mvel make_vgrids.o vdefs.o

# MOVE FILES TO RUN
copy :
	cp mk3d gendata slice ../Run/.
move :
	mv mk3d gendata slice ../Run/.

# CLEAN
clean :
	rm *.o
