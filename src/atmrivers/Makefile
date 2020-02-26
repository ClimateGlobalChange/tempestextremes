# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTEXTREMESDIR= ../..

# Load Makefile framework. 
include $(TEMPESTEXTREMESDIR)/mk/framework.make

TEMPESTEXTREMESBASEDIR= $(TEMPESTEXTREMESDIR)/src/base

TEMPESTEXTREMESBASELIB= $(TEMPESTEXTREMESBASEDIR)/libextremesbase.a

EXEC_FILES= SpineARs.cpp

EXEC_TARGETS= $(EXEC_FILES:%.cpp=%)

FILES= $(UTIL_FILES) $(EXEC_FILES)

.PHONY: all clean

# Build rules. 
all: $(EXEC_TARGETS)

$(EXEC_TARGETS): %: $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/%.o $(TEMPESTEXTREMESBASELIB)
	$(CXX) $(LDFLAGS) -o $@ $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/$*.o $(TEMPESTEXTREMESBASELIB) $(LIBRARIES)
	mv $@ $(TEMPESTEXTREMESDIR)/bin

$(TEMPESTEXTREMESBASELIB): %:
	cd $(TEMPESTEXTREMESDIR)/src/base; $(MAKE)

# Clean rules.
clean:
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
