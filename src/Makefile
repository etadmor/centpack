CPROOT = ..

COMMON = disclaimer.o minmod3.o
ARRAY  = doublearray1d.o doublearray2d.o doublearray3d.o
1D     = read_in_1d.o run_info_1d.o time_step_1d.o
1DFD2  = centpack_main_1d_FD2.o corrector_1d_FD2.o  end_of_step_1d_FD2.o \
	 evolution_1d_FD2.o predictor_1d_FD2.o reconstruction_1d_FD2.o
1DSD2  = centpack_main_1d_SD2.o C_flux_1d_SD2.o end_of_step_1d_SD2.o \
	 evolution_1d_SD2.o H_flux_1d_SD2.o reconstruction_1d_SD2.o
1DSD3  = centpack_main_1d_SD3.o C_flux_1d_SD3.o end_of_step_1d_SD3.o \
	 evolution_1d_SD3.o H_flux_1d_SD3.o indicators_1d_SD3.o \
         norm_1d_SD3.o reconstruction_1d_SD3.o sum_1d_SD3.o
2D     = read_in_2d.o run_info_2d.o time_step_2d.o
2DFD2  = centpack_main_2d_FD2.o corrector_2d_FD2.o  end_of_step_2d_FD2.o \
	 evolution_2d_FD2.o predictor_2d_FD2.o reconstruction_2d_FD2.o
2DSD2  = centpack_main_2d_SD2.o C_flux_2d_SD2.o end_of_step_2d_SD2.o \
	 evolution_2d_SD2.o Hx_flux_2d_SD2.o Hy_flux_2d_SD2.o \
	 reconstruction_2d_SD2.o
2DSD3  = centpack_main_2d_SD3.o C_flux_2d_SD3.o end_of_step_2d_SD3.o \
	 evolution_2d_SD3.o Hx_flux_2d_SD3.o Hy_flux_2d_SD3.o \
	 indicators_2d_SD3.o indicators_diag_2d_SD3.o \
         norm_2d_SD3.o reconstruction_2d_SD3.o reconstruction_diag_2d_SD3.o \
         sum_2d_SD3.o
CP     = libcentpack_1d_FD2.a libcentpack_1d_SD2.a libcentpack_1d_SD3.a \
         libcentpack_2d_FD2.a libcentpack_2d_SD2.a libcentpack_2d_SD3.a \
         libarray.a

COMMON := $(addprefix common/,    $(COMMON))
ARRAY  := $(addprefix common/,    $(ARRAY))
1D     := $(addprefix 1d/common/, $(1D))
1DFD2  := $(addprefix 1d/FD2/,    $(1DFD2))
1DSD2  := $(addprefix 1d/SD2/,    $(1DSD2))
1DSD3  := $(addprefix 1d/SD3/,    $(1DSD3))
2D     := $(addprefix 2d/common/, $(2D))
2DFD2  := $(addprefix 2d/FD2/,    $(2DFD2))
2DSD2  := $(addprefix 2d/SD2/,    $(2DSD2))
2DSD3  := $(addprefix 2d/SD3/,    $(2DSD3))
CP     := $(addprefix $(CPROOT)/lib/, $(CP))

OBJS =  $(COMMON) $(1D) $(1DFD2) $(1DSD2) $(1DSD3) \
	$(2D) $(2DFD2) $(2DSD2) $(2DSD3) $(ARRAY)

libs: $(CP)

samples:
	(cd 1d/samples; make samples)
	(cd 2d/samples; make samples)

clean:
	rm $(ARRAY)
	rm $(COMMON)
	rm $(1D) $(1DFD2) $(1DSD2) $(1DSD3)
	rm $(2D) $(2DFD2) $(2DSD2) $(2DSD3)
	
clean_samples:
	(cd 1d/samples; make clean)
	(cd 2d/samples; make clean)

$(OBJS): %.o : %.cc
	g++ -c -fPIC -I../include $< -o $@

$(CPROOT)/lib/libarray.a: $(ARRAY)
	ar cr $@ $(ARRAY)

$(CPROOT)/lib/libcentpack_1d_FD2.a: $(COMMON) $(1D) $(1DFD2)
	ar cr $@ $(COMMON) $(1D) $(1DFD2)

$(CPROOT)/lib/libcentpack_1d_SD2.a: $(COMMON) $(1D) $(1DSD2)
	ar cr $@ $(COMMON) $(1D) $(1DSD2)

$(CPROOT)/lib/libcentpack_1d_SD3.a: $(COMMON) $(1D) $(1DSD3)
	ar cr $@ $(COMMON) $(1D) $(1DSD3)

$(CPROOT)/lib/libcentpack_2d_FD2.a: $(COMMON) $(2D) $(2DFD2)
	ar cr $@ $(COMMON) $(2D) $(2DFD2)

$(CPROOT)/lib/libcentpack_2d_SD2.a: $(COMMON) $(2D) $(2DSD2)
	ar cr $@ $(COMMON) $(2D) $(2DSD2)

$(CPROOT)/lib/libcentpack_2d_SD3.a: $(COMMON) $(2D) $(2DSD3)
	ar cr $@ $(COMMON) $(2D) $(2DSD3)

