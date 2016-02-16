# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Rules for compiling object files

# Compilation directories
DEPDIR= $(CURDIR)/depend
BUILDDIR= $(CURDIR)/build

# Dependency file construction
MAKEDEPENDCPP= \
    mkdir -p $(DEPDIR); \
    echo "-- Generating dependencies for $<"; \
    $(CXX) -M $(CXXFLAGS) $(CURDIR)/$< > $(DEPDIR)/$*.P; \
    sed -e 's~.*:~$(BUILDDIR)/$*.o $(DEPDIR)/$*.d:~' < $(DEPDIR)/$*.P > $(DEPDIR)/$*.d; \
    sed -e 's/.*://' -e 's/\\$$//' < $(DEPDIR)/$*.P | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(DEPDIR)/$*.d; \
    rm -f $(DEPDIR)/$*.P

# Compilation rules
$(BUILDDIR)/%.o : %.cpp
	@mkdir -p $(@D)
	@$(MAKEDEPENDCPP)
	$(CXX) $(CXXFLAGS) -c -o $@ $(CURDIR)/$<

$(BUILDDIR)/%.o : %.f90
	@mkdir -p $(@D)
	$(F90) $(F90FLAGS) -c -o $@ $(CURDIR)/$<

# DO NOT DELETE
