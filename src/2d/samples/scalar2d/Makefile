APPNAME   = scalar2d_2d_SD3

CPROOT    = ../../../..
SAMPLEDIR = $(CPROOT)/samples
APPDIR    = $(SAMPLEDIR)/$(APPNAME)

OBJS = boundary_conditions.o flux_x.o flux_y.o initial_conditions.o \
       spectral_radii.o writeout.o ../sample_mains/main_2d_SD3.o

DIRS = $(SAMPLEDIR) $(APPDIR) $(APPDIR)/u_files $(APPDIR)/t_files \
$(APPDIR)/u_frames

samples all: $(APPNAME)

clean:
	rm $(OBJS)
	
$(DIRS):
	mkdir $@

$(OBJS): %.o : %.cc
	g++ -c -fPIC -I$(CPROOT)/include $< -o $@

$(APPNAME): $(DIRS) $(OBJS)
	g++ $(OBJS) -I$(CPROOT)/include -L$(CPROOT)/lib -lcentpack_2d_SD3 -larray -o $@ -O3

	mv $(APPNAME) $(APPDIR)
	cp input $(APPDIR)
	cp scalar2d_2d_frames.m $(APPDIR)
