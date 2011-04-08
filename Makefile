include Make.inc

all: library 

library: libdir srcd

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
srcd:
	(cd src; make lib)

install:
	(./mkdir.sh  $(INSTALL_DIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_DIR))
	(./mkdir.sh $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR))
	(./mkdir.sh $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) lib/*$(.mod) $(INSTALL_INCLUDEDIR))
	(./mkdir.sh  $(INSTALL_DOCSDIR) && \
	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR))
veryclean: 
	(cd src; make veryclean)
	(cd lib; /bin/rm -fr *.a *$(.mod))
	(cd tests/pdegen; make clean)

clean:
	(cd src; make clean)
