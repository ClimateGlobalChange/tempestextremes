##############################################################################
# Compiler and flags
CC= g++
CCOMP= gcc
CFLAGS= -O3

USEBLAS= True

# NETCDF library directories
NETCDF_INCLUDEDIR=$(NETCDF_DIR)/include
NETCDF_LIBDIR=$(NETCDF_DIR)/lib

# Library files to include
LDFILES= -lnetcdf -lnetcdf_c++ -L. 

##############################################################################
# DO NOT MODIFY BELOW THIS LINE
##############################################################################

# Local files
######################################
#ADDED MY FILES HERE!
######################################
CLIVAR_FILES = blockingPV.cpp Announce.cpp NetCDFUtilities.cpp blockingUtilities.cpp interpolate.cpp TimeObj.cpp
CLIVAR_CFILES =

BLOCK_FILES = blockingAvg.cpp Announce.cpp NetCDFUtilities.cpp TimeObj.cpp  blockingUtilities.cpp 
BLOCK_CFILES = 

DEV_FILES = blockingDevs.cpp Announce.cpp NetCDFUtilities.cpp blockingUtilities.cpp TimeObj.cpp
DEV_CFILES = 

DENS_FILES = densityCalculations.cpp Announce.cpp NetCDFUtilities.cpp blockingUtilities.cpp TimeObj.cpp
DENS_CFILES = 

GH_FILES = blockingGH.cpp Announce.cpp NetCDFUtilities.cpp interp_z500.cpp blockingUtilities.cpp TimeObj.cpp
GH_CFILES =  
######################################

STITCHNODES_FILES= StitchNodes.cpp Announce.cpp
STITCHNODES_CFILES= kdtree.c

STITCHBLOBS_FILES= StitchBlobs.cpp Announce.cpp NetCDFUtilities.cpp
STITCHBLOBS_CFILES= kdtree.c

DETECTCYCLONES_FILES= DetectCyclones.cpp Announce.cpp TimeObj.cpp
DETECTCYCLONES_CFILES= kdtree.c

DENSITYNODES_FILES= DensityNodes.cpp Announce.cpp
DENSITYNODES_CFILES=

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

include Make.defs

##
## Build instructions
##
all: StitchNodes StitchBlobs DetectCyclones DensityNodes blockingPV blockingAvg blockingDevs blockingDensity blockingGH

StitchNodes: $(STITCHNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHNODES_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(STITCHNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHNODES_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

StitchBlobs: $(STITCHBLOBS_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHBLOBS_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(STITCHBLOBS_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHBLOBS_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

DetectCyclones: $(DETECTCYCLONES_FILES:%.cpp=$(BUILDDIR)/%.o) $(DETECTCYCLONES_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(DETECTCYCLONES_FILES:%.cpp=$(BUILDDIR)/%.o) $(DETECTCYCLONES_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

DensityNodes: $(DENSITYNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(DENSITYNODES_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(DENSITYNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(DENSITYNODES_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

blockingPV: $(CLIVAR_FILES:%.cpp=$(BUILDDIR)/%.o) $(CLIVAR_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(CLIVAR_FILES:%.cpp=$(BUILDDIR)/%.o) $(CLIVAR_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

blockingAvg: $(BLOCK_FILES:%.cpp=$(BUILDDIR)/%.o) $(BLOCK_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(BLOCK_FILES:%.cpp=$(BUILDDIR)/%.o) $(BLOCK_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)
blockingDevs: $(DEV_FILES:%.cpp=$(BUILDDIR)/%.o) $(DEV_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(DEV_FILES:%.cpp=$(BUILDDIR)/%.o) $(DEV_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)
blockingDensity: $(DENS_FILES:%.cpp=$(BUILDDIR)/%.o) $(DENS_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(DENS_FILES:%.cpp=$(BUILDDIR)/%.o) $(DENS_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

blockingGH: $(GH_FILES:%.cpp=$(BUILDDIR)/%.o) $(GH_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GH_FILES:%.cpp=$(BUILDDIR)/%.o) $(GH_CFILES:%/c=$(BUILDDIR)/%.o) $(LDFILES)

## Clean
##
clean:
	rm -f StitchNodes StitchBlobs DetectCyclones DensityNodes blockingPV blockingAvg blockingDevs blockingDensity blockingGH
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)
include $(FILES:%.c=$(DEPDIR)/%.d)

# DO NOT DELETE

