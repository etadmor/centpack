libs:
	(cd src; make libs)

samples:
	(cd src; make samples)
	
clean:
	(cd src; make clean)
	
clean_samples:
	(cd src; make clean_samples)
	
burgers_1d_SD3:
	(cd src/1d/samples/burgers; make samples)
	
euler_1d_SD2:
	(cd src/1d/samples/euler; make samples)

MHD_1d_FD2:
	(cd src/1d/samples/MHD_BW; make samples)	
	
scalar2d_2d_SD3:
	(cd src/2d/samples/scalar2d; make samples)
	
euler_2d_FD2:
	(cd src/2d/samples/euler; make samples)

MHD_2d_SD2:
	(cd src/2d/samples/MHD_OT; make samples)
	
clean_burgers_1d_SD3:
	(cd src/1d/samples/burgers; make clean)
	
clean_euler_1d_SD2:
	(cd src/1d/samples/euler; make clean)

clean_MHD_1d_FD2:
	(cd src/1d/samples/MHD_BW; make clean)	
	
clean_scalar2d_2d_SD3:
	(cd src/2d/samples/scalar2d; make clean)
	
clean_euler_2d_FD2:
	(cd src/2d/samples/euler; make clean)

clean_MHD_2d_SD2:
	(cd src/2d/samples/MHD_OT; make clean)

