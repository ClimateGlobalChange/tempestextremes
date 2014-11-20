##############################################################################
# Compiler and flags
CC= g++
CCOMP= gcc
CFLAGS= -O3

USEBLAS= True

# NETCDF library directories
NETCDF_INCLUDEDIR=/opt/local/include
NETCDF_LIBDIR=/opt/local/lib

# Library files to include
LDFILES= -lnetcdf -lnetcdf_c++ -framework accelerate

##############################################################################
# DO NOT MODIFY BELOW THIS LINE
##############################################################################

# Local files
STITCHNODES_FILES= StitchNodes.cpp Announce.cpp
STITCHNODES_CFILES= kdtree.c

STITCHBLOBS_FILES= StitchBlobs.cpp Announce.cpp
STITCHBLOBS_CFILES= kdtree.c

DETECTCYCLONES_FILES= DetectCyclones.cpp Announce.cpp TimeObj.cpp
DETECTCYCLONES_CFILES= kdtree.c

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

include Make.defs

##
## Build instructions
##
all: StitchNodes StitchBlobs DetectCyclones

StitchNodes: $(STITCHNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHNODES_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(STITCHNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHNODES_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

StitchBlobs: $(STITCHBLOBS_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHBLOBS_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(STITCHBLOBS_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHBLOBS_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

DetectCyclones: $(DETECTCYCLONES_FILES:%.cpp=$(BUILDDIR)/%.o) $(DETECTCYCLONES_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(DETECTCYCLONES_FILES:%.cpp=$(BUILDDIR)/%.o) $(DETECTCYCLONES_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f StitchNodes StitchBlobs DetectCyclones *.o
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)
include $(FILES:%.c=$(DEPDIR)/%.d)

# DO NOT DELETE

