# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTEXTREMESDIR= ../..

# Load Makefile framework. 
include $(TEMPESTEXTREMESDIR)/mk/framework.make

FILES= Announce.cpp \
       NetCDFUtilities.cpp \
       TimeObj.cpp \
       Variable.cpp \
	   DataOp.cpp \
       kdtree.cpp \
	   SimpleGridUtilities.cpp \
	   AutoCurator.cpp \
	   ArgumentTree.cpp \
	   NodeFileUtilities.cpp \
	   RLLPolygonArray.cpp \
	   SimpleGrid.cpp

LIB_TARGET= libextremesbase.a

.PHONY: all clean

# Build rules. 
all: $(LIB_TARGET)

$(LIB_TARGET): $(FILES:%.cpp=$(BUILDDIR)/%.o)
	rm -f $(LIB_TARGET)
	ar -cqs $(LIB_TARGET) build/*.o

# Clean rules.
clean:
	rm -f $(LIB_TARGET) 
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
