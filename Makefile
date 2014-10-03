##############################################################################
# Compiler and flags
CC= g++
CCOMP= gcc
CFLAGS= -O3

USEBLAS= True

# NETCDF library directories
NETCDF_INCLUDEDIR=/usr/local/include
NETCDF_LIBDIR=/usr/local/lib

# Library files to include
LDFILES= -framework accelerate

##############################################################################
# DO NOT MODIFY BELOW THIS LINE
##############################################################################

# Local files
STITCHNODES_FILES= StitchNodes.cpp Announce.cpp
STITCHNODES_CFILES= kdtree.c

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

include Make.defs

##
## Build instructions
##
all: StitchNodes

StitchNodes: $(STITCHNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHNODES_CFILES:%.c=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(STITCHNODES_FILES:%.cpp=$(BUILDDIR)/%.o) $(STITCHNODES_CFILES:%.c=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f StitchNodes *.o
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)
include $(FILES:%.c=$(DEPDIR)/%.d)

# DO NOT DELETE

