#---------------------------------------------------------------------------#
# Makefile for Pessoa 							     #
# University of California, Los Angeles                                      #
# Version 2.0, November, 2015	                                     	     #
#----------------------------------------------------------------------------#

CUDD_DIR      = cudd-2.5.0
SRC_DIRS      = library Pessoa_abstract Pessoa_control Pessoa_design utils
MAKEARG       = 'CUDD_DIR=../$(CUDD_DIR)'



.PHONY: all clean cudd

all:   cudd pessoa2.out 

clean:
	$(MAKE) $@ -C $(CUDD_DIR)
	

distclean: clean
	$(MAKE) $@ -C $(CUDD_DIR)

cudd:
	$(MAKE) -C $(CUDD_DIR) 

pessoa2.out: cudd
	$(MAKE) $(MAKEARG) -C ./src

  


##



