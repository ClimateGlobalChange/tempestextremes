# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Detects the system and includes the system-specific makefile.

UNAME := $(shell uname)

ifeq ($(UNAME),Darwin)
  SYSTEM= MACOSX
  SYSTEM_MAKEFILE= macosx.make
else ifeq ($(UNAME),Linux)
  ifeq ($(NERSC_HOST),babbage)
    SYSTEM= BABBAGE
    SYSTEM_MAKEFILE= babbage.make
  endif
  ifeq ($(HOSTNAME),yslogin1)
    SYSTEM= YELLOWSTONE
    SYSTEM_MAKEFILE= yellowstone.make
  endif
  ifeq ($(HOSTNAME),yslogin2)
    SYSTEM= YELLOWSTONE
    SYSTEM_MAKEFILE= yellowstone.make
  endif
  ifeq ($(HOSTNAME),yslogin3)
    SYSTEM= YELLOWSTONE
    SYSTEM_MAKEFILE= yellowstone.make
  endif
  ifeq ($(HOSTNAME),yslogin4)
    SYSTEM= YELLOWSTONE
    SYSTEM_MAKEFILE= yellowstone.make
  endif
  ifeq ($(HOSTNAME),yslogin5)
    SYSTEM= YELLOWSTONE
    SYSTEM_MAKEFILE= yellowstone.make
  endif
  ifeq ($(HOSTNAME),yslogin6)
    SYSTEM= YELLOWSTONE
    SYSTEM_MAKEFILE= yellowstone.make
  endif
  ifeq ($(SYSTEM),)
    SYSTEM= AGRI
    SYSTEM_MAKEFILE= agri.make
  endif   
endif

include $(TEMPESTEXTREMESDIR)/mk/system/$(SYSTEM_MAKEFILE)

# Build identifier
BUILDID:= $(SYSTEM)

ifeq ($(OPT),TRUE)
  BUILDID:=$(BUILDID).OPT
endif

ifeq ($(DEBUG),TRUE)
  BUILDID:=$(BUILDID).DEBUG
endif

ifeq ($(PARALLEL),MPIOMP)
  BUILDID:=$(BUILDID).MPIOMP
else ifeq ($(PARALLEL),HPX)
  BUILDID:=$(BUILDID).HPX
endif

# DO NOT DELETE
