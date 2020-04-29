include Make.inc

all: library 

library: libdir srcd

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi)
	(if test ! -d modules ; then mkdir modules; fi)
	($(INSTALL_DATA) Make.inc  include/Make.inc.mld-ainv)
srcd:
	$(MAKE) -C src

install:
	(mkdir -p $(INSTALL_DIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_DIR))
	(mkdir -p $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR))
	(mkdir -p $(INSTALL_MODULESDIR) && \
	   $(INSTALL_DATA) modules/*$(.mod) $(INSTALL_MODULESDIR))
	(mkdir -p $(INSTALL_INCLUDEDIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_INCLUDEDIR)/Make.inc.mld-ainv)
#	(mkdir -p $(INSTALL_DOCSDIR) && \
#	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR))

veryclean: 
	$(MAKE) -C src veryclean
	(cd lib; /bin/rm -fr *.a *$(.mod))
	(cd tests/pdegen; make clean)

clean:
	(cd src; make clean)
