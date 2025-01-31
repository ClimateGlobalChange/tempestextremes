# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

DEFAULT_BUILD_TARGETS= src/netcdf-cxx-4.2 src/base src/blobs src/nodes src/util
ALL_BUILD_TARGETS= $(DEFAULT_BUILD_TARGETS) src/blocking src/sandbox
CLEAN_TARGETS= $(addsuffix .clean,$(ALL_BUILD_TARGETS))

.PHONY: all clean $(ALL_BUILD_TARGETS) $(CLEAN_TARGETS)

# Build rules.
all: $(DEFAULT_BUILD_TARGETS)

$(ALL_BUILD_TARGETS): %:
	cd $*; $(MAKE)

# Clean rules.
clean: $(CLEAN_TARGETS)
	rm -f bin/*

$(CLEAN_TARGETS): %.clean:
	cd $*; $(MAKE) clean

# DO NOT DELETE
