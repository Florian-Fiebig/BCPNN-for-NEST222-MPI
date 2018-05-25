# Automake file for external dynamic modules for NEST
#
# Hans Ekkehard Plesser, April 2008
# Automake file for the Developer Module
# 
# libpt_module is built as a normal, installable library.
# It will be installed to $prefix/lib by make install.
# 
# Headers from this directory are not to be installed upon
# make install. They are therefore included in _SOURCES.


# 1. Exchange "my" in "pt_module" with the name of your model below
#    (ten locations).
# 2. Add all .cpp and .h files from your code as *_SOURCES. Header files
#    are given only so that they will be included in the tarball if you
#    run "make dist" on your module.
# 3. The ml_module* stuff creates a module that can be loaded at runtime.
#    It is called pt_module.so.
# 4. The libpt_module* stuff creates a library against which NEST can be
#    linked.

libdir= @libdir@/nest

lib_LTLIBRARIES=      pt_module.la libpt_module.la

pt_module_la_CXXFLAGS= @AM_CXXFLAGS@
pt_module_la_SOURCES =  pt_module.cpp                   pt_module.h         	\
                        bcpnn_connection.cpp            bcpnn_connection.h      \
                        iaf_cond_alpha_bias.cpp         iaf_cond_alpha_bias.h	\
                        aeif_cond_exp_multisynapse.cpp         aeif_cond_exp_multisynapse.h	\
                        iaf_cond_exp_bias.cpp         iaf_cond_exp_bias.h
pt_module_la_LDFLAGS=  -module

libpt_module_la_CXXFLAGS= $(pt_module_la_CXXFLAGS) -DLINKED_MODULE
libpt_module_la_SOURCES=  $(pt_module_la_SOURCES)

MAKEFLAGS= @MAKE_FLAGS@

AM_CPPFLAGS= @NEST_CPPFLAGS@ \
             @INCLTDL@      

.PHONY: install-slidoc

nobase_pkgdata_DATA=\
	sli/pt_module.sli

install-slidoc:
	NESTRCFILENAME=/dev/null $(DESTDIR)$(NEST_PREFIX)/bin/sli --userargs="@HELPDIRS@" $(NEST_PREFIX)/share/nest/sli/install-help.sli

install-data-hook: install-exec install-slidoc

EXTRA_DIST= sli
