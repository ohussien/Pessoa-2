# $Id$
#
# Makefile for the CUDD distribution kit
#---------------------------------------------------------------------------

# Beginning of the configuration section. These symbol definitions can
# be overridden from the command line.

# C++ compiler
CXX	= g++

# Specific options for compilation of C++ files.
CXXFLAGS = 
# Stricter standard conformance for g++.
#CXXFLAGS = -std=c++98
# For Sun CC version 5, this invokes compatibility mode.
#CXXFLAGS = -compat
# On some versions of UP-UX, it is necessary to pass the option +a1
# to CC for the C++ test program to compile successfully.
#CXXFLAGS = +a1

# C compiler used for all targets except optimize_dec, which always uses cc.
CC = gcc

# On some machines ranlib is either non-existent or redundant.
# Use the following definition if your machine has ranlib and you think
# it is needed.
RANLIB	= ranlib
# Use the following definition if your machine either does not have
# ranlib (e.g., SUN running solaris) or can do without it (e.g., DEC Alpha).
#RANLIB	= :

# Use ICFLAGS to specify machine-independent compilation flags.
# These two are typical settings for optimized code with gcc.
#ICFLAGS	= -g -O3 -Wall
ICFLAGS = -g -fPIC


# Use XCFLAGS to specify machine-dependent compilation flags.
# For some platforms no special flags are needed.
#XCFLAGS	= -DHAVE_IEEE_754 -DBSD
XCFLAGS = -fPIC -DPIC -DHAVE_IEEE_754 -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 -DHAVE_IEEE_754 -fno-common

# Define the level of self-checking and verbosity of the CUDD package.
#DDDEBUG = -DDD_DEBUG -DDD_VERBOSE -DDD_STATS -DDD_CACHE_PROFILE -DDD_UNIQUE_PROFILE -DDD_COUNT
DDDEBUG	=

# Define the level of self-checking and verbosity of the MTR package.
#MTRDEBUG = -DMTR_DEBUG
MTRDEBUG =

# Loader options.
LDFLAGS	=
# This may produce faster code on the DECstations.
#LDFLAGS	= -jmpopt -Olimit 1000
# This may be necessary under some old versions of Linux.
#LDFLAGS	= -static
# This normally makes the program faster on the DEC Alphas.
#LDFLAGS	= -non_shared -om
# This is for 32-bit pointers on the DEC Alphas.
#LDFLAGS	= -non_shared -om -taso
#LDFLAGS	= -non_shared -taso

# Define PURE as purify to link with purify.
# Define PURE as quantify to link with quantify.
# Remember to compile with -g if you want line-by-line info with quantify.
PURE =
#PURE	= purify
#PURE	= quantify

# Define EXE as .exe for MS-DOS and derivatives.  Not required by recent
# versions of cygwin.
EXE	=
#EXE	= .exe

# End of the configuration section.
#---------------------------------------------------------------------------

MFLAG   = -DMNEMOSYNE
MNEMLIB	= ../mnemosyne/libmnem.a

DDWDIR	= .
IDIR	= $(DDWDIR)/include
INCLUDE = -I$(IDIR)

BDIRS	= cudd dddmp mtr st util epd
DIRS	= $(BDIRS) nanotrav

#------------------------------------------------------------------------

.PHONY : build
.PHONY : nanotrav
.PHONY : check_leaks
.PHONY : optimize_dec
.PHONY : testcudd
.PHONY : libobj
.PHONY : testobj
.PHONY : testdddmp
.PHONY : testmtr
.PHONY : lint
.PHONY : all
.PHONY : clean
.PHONY : distclean


build:
	sh ./setup.sh
	@for dir in $(DIRS); do \
		(cd $$dir; \
		echo Making $$dir ...; \
		make CC=$(CC) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" )\
	done

nanotrav: build

check_leaks:
	sh ./setup.sh
	@for dir in mnemosyne $(DIRS); do \
		(cd $$dir; \
		echo Making $$dir ...; \
		make CC=$(CC) RANLIB=$(RANLIB) MFLAG=$(MFLAG) MNEMLIB=$(MNEMLIB) ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" EXE="$(EXE)" )\
	done

optimize_dec:
	sh ./setup.sh
	@for dir in $(DIRS); do \
		(cd $$dir; \
		echo Making $$dir ...; \
		make CC=$(CC) RANLIB=$(RANLIB) XCFLAGS="$(XCFLAGS)" LDFLAGS="$(LDFLAGS)" optimize_dec )\
	done

lint:
	sh ./setup.sh
	@for dir in $(DIRS) obj; do \
		(cd $$dir; \
		echo Making lint in $$dir ...; \
		make CC=$(CC) lint )\
	done

tags:
	sh ./setup.sh
	@for dir in $(DIRS) obj; do \
		(cd $$dir; \
		echo Making tags in $$dir ...; \
		make CC=$(CC) tags )\
	done

all:
	sh ./setup.sh
	@for dir in $(DIRS); do \
		(cd $$dir; \
		echo Making all in $$dir ...; \
		make CC=$(CC) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" all )\
	done

testcudd:
	sh ./setup.sh
	@for dir in util st mtr epd; do \
		(cd $$dir; \
		echo Making $$dir ...; \
		make CC=$(CC) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" )\
	done
	@(cd cudd; \
	echo Making testcudd ...; \
	make CC=$(CC) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" testcudd$(EXE) )

objlib:
	sh ./setup.sh
	@for dir in $(BDIRS); do \
		(cd $$dir; \
		echo Making $$dir ...; \
		make CC=$(CC) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" )\
	done
	@(cd obj; \
	echo Making obj ...; \
	make CXX=$(CXX) CXXFLAGS=$(CXXFLAGS) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" )

testobj: objlib
	@(cd obj; \
	echo Making testobj ...; \
	make CXX=$(CXX) CXXFLAGS=$(CXXFLAGS) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" testobj$(EXE) )

testdddmp: build
	@(cd dddmp; \
	echo Making testdddmp ...; \
	make CC=$(CC) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" testdddmp$(EXE) )

testmtr: build
	@(cd mtr; \
	echo Making testmtr ...; \
	make CC=$(CC) RANLIB=$(RANLIB) MFLAG= MNEMLIB= ICFLAGS="$(ICFLAGS)" XCFLAGS="$(XCFLAGS)" DDDEBUG="$(DDDEBUG)" MTRDEBUG="$(MTRDEBUG)" LDFLAGS="$(LDFLAGS)" PURE="$(PURE)" EXE="$(EXE)" testmtr$(EXE) )

clean:
	@for dir in mnemosyne $(DIRS) obj; do	\
	    (cd $$dir;	\
	     echo Cleaning $$dir ...; \
	     make -s clean	) \
	done

distclean:
	@for dir in mnemosyne $(DIRS) obj; do	\
	    (cd $$dir;	\
	     echo Cleaning $$dir ...; \
	     make -s EXE="$(EXE)" distclean	) \
	done
	sh ./shutdown.sh
