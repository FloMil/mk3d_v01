fm3d_objects = 3dfmlib.o propagate.o rays.o frechet.o visual.o 3dfm_main.o matchref.o teleseismic.o
fm3d_modules = mod_3dfm.o
tt_objects = libsun.o libtau.o ellip.o sphdist.o
nn_objects = nn_subsf.o stack.o
svd_objects = svdlib.o

f90comp = gfortran

%.o : %.f90
	$(f90comp) -c  $< -O3

%.o : %.f
	$(f90comp) -c  $< -O3

fm3d : $(fm3d_objects) $(fm3d_modules) $(svd_objects) $(nn_objects) $(tt_objects)
	$(f90comp) -o fm3d  $(fm3d_modules) $(fm3d_objects) $(svd_objects) $(nn_objects) $(tt_objects) -O3

$(fm3d_objects) : $(fm3d_modules)

clean:
	rm $(fm3d_objects) $(fm3d_modules) $(svd_objects) $(nn_objects) $(tt_objects)

mint : make_interfaces.o idefs.o
	$(f90comp) -o mint make_interfaces.o idefs.o

mvel : make_vgrids.o vdefs.o
	$(f90comp) -o mvel make_vgrids.o vdefs.o

copy :
	cp mint mvel fm3d ../../Run/fm3d/.
move :
	mv mint mvel fm3d ../../Run/fm3d/.
