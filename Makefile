include Make.inc

all: library 

library: libdir srcd

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi)
	($(INSTALL_DATA) Make.inc  include/Make.inc.mldainv)
srcd:
	cd src && $(MAKE)

install:
	(./mkdir.sh  $(INSTALL_DIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_DIR))
	(./mkdir.sh $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR))
	(./mkdir.sh $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) include/*$(.mod) $(INSTALL_INCLUDEDIR))
	(./mkdir.sh  $(INSTALL_INCLUDEDIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_INCLUDEDIR)/Make.inc.mld-ainv)
#	(./mkdir.sh  $(INSTALL_DOCSDIR) && \
#	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR))

veryclean: 
	cd src && $(MAKE) veryclean
	(cd lib; /bin/rm -fr *.a *$(.mod))
	(cd tests/pdegen; make clean)

clean:
	(cd src; make clean)
