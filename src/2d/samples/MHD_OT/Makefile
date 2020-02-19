APPNAME   = MHD_2d_SD2

CPROOT    = ../../../..
SAMPLEDIR = $(CPROOT)/samples
APPDIR    = $(SAMPLEDIR)/$(APPNAME)

OBJS = boundary_conditions.o flux_x.o flux_y.o initial_conditions.o \
       spectral_radii.o writeout.o ../sample_mains/main_2d_SD2.o

DIRS = $(SAMPLEDIR) $(APPDIR) $(APPDIR)/rho_files $(APPDIR)/u1_files \
       $(APPDIR)/u2_files $(APPDIR)/u3_files $(APPDIR)/b1_files \
	   $(APPDIR)/b2_files $(APPDIR)/b3_files $(APPDIR)/p_files \
	   $(APPDIR)/t_files $(APPDIR)/rho_frames $(APPDIR)/u_frames \
       $(APPDIR)/b_frames $(APPDIR)/p_frames

samples all: $(APPNAME)

clean:
	rm $(OBJS)

$(DIRS):
	mkdir $@

$(OBJS): %.o : %.cc
	g++ -c -fPIC -I$(CPROOT)/include $< -o $@

$(APPNAME): $(DIRS) $(OBJS)
	g++ $(OBJS) -I$(CPROOT)/include -L$(CPROOT)/lib -lcentpack_2d_SD2 -larray -o $@ -O3

	mv $(APPNAME) $(APPDIR)
	cp input $(APPDIR)
	cp MHD_2d_frames.m $(APPDIR)
