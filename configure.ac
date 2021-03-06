dnl $Id$

dnl
dnl Michele Martone's proposal for a transition towards the GNU build system for the PSBLAS library
dnl  
 
dnl Process this file with autoconf to produce a configure script.
dnl
dnl usage : aclocal -I config/ && autoconf && ./configure && make 

dnl then  : VAR=VAL ./configure

dnl In some configurations (AIX) the next line is needed:
dnl MPIFC=mpxlf95 ./configure

dnl then  : ./configure VAR=VAL 
dnl then  : ./configure --help=short
dnl then  : ./configure --help

dnl the PSBLAS modules get this task difficult to accomplish!
dnl SEE  : --module-path --include-path

dnl NOTE : There is no cross compilation support.



###############################################################################
###############################################################################
# 
#  This script is used by the PSBLAS to determine the compilers, linkers, and
# libraries to build its libraries executable code.
#  Its behaviour is driven on the compiler it finds or it is dictated to work
# with.
#
###############################################################################
###############################################################################

# NOTE: the literal for version (the second argument to AC_INIT should be a literal!)
AC_INIT([MLD_AINV],1.0, salvatore.filippone@uniroma2.it)

# VERSION is the file containing the PSBLAS version code
# FIXME
mld_ainv_cv_version="1.0"
# A sample source file
dnl AC_CONFIG_AUX_DIR(config)

# A sample source file
AC_CONFIG_SRCDIR([src/mld_base_ainv_mod.F90])
# Our custom M4 macros are in  the 'config' directory
AC_CONFIG_MACRO_DIR([config]) 
AC_MSG_NOTICE([                                             
 --------------------------------------------------------------------------------
     	Welcome to the $PACKAGE_NAME $mld_ainv_cv_version configure Script.
             
 	This  creates Make.inc, but if you read carefully the
	docum entation, you can make your own by hand for your needs.

        ./configure  --with-mld2p4=/path/to/mld2p4 
	See ./configure --help=short fore more info.
--------------------------------------------------------------------------------
		])

###############################################################################
# 			FLAGS and LIBS user customization
###############################################################################

dnl NOTE : no spaces before the comma, and no brackets before the second argument!
PAC_ARG_WITH_MLD2P4
MLD2P4_DIR="$pac_cv_mld2p4_dir";
AC_MSG_CHECKING([for MLD2P4 install dir])
case $MLD2P4_DIR in 
     /*) ;; 
     *) AC_MSG_ERROR([The MLD2P4 installation dir must be an absolute pathname
 specified with --with-mld2p4=/path/to/mld2p4])
esac
if test ! -d "$MLD2P4_DIR" ; then 
   AC_MSG_ERROR([Could not find MLD2P4 build dir $MLD2P4_DIR!])	
fi
AC_MSG_RESULT([$MLD2P4_DIR])
pac_cv_status_file="$MLD2P4_DIR/include/Make.inc.mld2p4"
if test ! -f "$pac_cv_status_file" ; then 
   AC_MSG_ERROR([Could not find an installation  in $MLD2P4_DIR.])		
fi
dnl. $pac_cv_status_file
dnl AC_MSG_NOTICE([Loaded $pac_cv_status_file $FC $MPIFC $BLACS_LIBS])
dnl AM_INIT_AUTOMAKE
dnl Specify required version of autoconf.
AC_PREREQ(2.59)
#
# Installation. 
#
#
AC_PROG_INSTALL

AC_MSG_CHECKING([where to install])
case $prefix in  
   \/* )   eval "INSTALL_DIR=$prefix";;
   * ) eval "INSTALL_DIR=/usr/local/mld_ainv";;
esac
case $libdir in 
   \/* )   eval "INSTALL_LIBDIR=$libdir";;
   * ) eval "INSTALL_LIBDIR=$INSTALL_DIR/lib";;
esac
case $includedir in 
   \/* )   eval "INSTALL_INCLUDEDIR=$includedir";;
   * ) eval "INSTALL_INCLUDEDIR=$INSTALL_DIR/include";;
esac
INSTALL_MODULESDIR=$INSTALL_DIR/modules
case $docsdir in 
   \/* )   eval "INSTALL_DOCSDIR=$docsdir";;
   * ) eval "INSTALL_DOCSDIR=$INSTALL_DIR/docs";;
esac
AC_MSG_RESULT([$INSTALL_DIR $INSTALL_INCLUDEDIR $INSTALL_LIBDIR $INSTALL_DOCSDIR])

###############################################################################
# Compilers detection: FC,F77,CC should be set, if found.
###############################################################################

AC_PROG_FC([ftn xlf2003 xlf95 xlf90 xlf pgf95 pgf90 ifort ifc nagfor  gfortran])
AC_PROG_CC([xlc pgcc icc gcc cc])
dnl AC_PROG_CXX

if test "X$CC" == "X" ; then
	AC_MSG_ERROR([Problem : No C compiler specified nor found!])
fi
if eval "$FC -qversion 2>&1 | grep XL 2>/dev/null" ; then
	# Some configurations of the XLF want "-WF," prepended to -D.. flags.
	# TODO : discover the exact conditions when the usage of -WF is needed.
	mld_cv_define_prepend="-WF,"
        if eval "$MPIFC -qversion 2>&1 | grep -e\"Version: 10\.\" 2>/dev/null"; then
        	FDEFINES="$mld_cv_define_prepend-DXLF_10 $FDEFINES"		  
        fi

 # Note : there coule be problems with old xlf compiler versions ( <10.1 )
 # since (as far as it is known to us) -WF, is not used in earlier versions.
 # More problems could be undocumented yet.
fi 
PAC_ARG_WITH_LIBS

AC_LANG([C])
###############################################################################
# First check whether the user required our serial (fake) mpi.
PAC_ARG_SERIAL_MPI

#Note : we miss the name of the Intel C compiler
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   FAKEMPI="fakempi.o";
   MPIFC="$FC";
   MPIF77="$F77";
   MPICC="$CC";
else 
if test "X$MPICC" = "X" ; then
    # This is our MPICC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPICC],[mpxlc mpcc mpicc cc])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for C]])])
AC_LANG(Fortran 77)
if test "X$MPIF77" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIF77],[mpxlf mpf77 mpif77 ftn])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran 77]])])

AC_LANG([Fortran])

if test "X$MPIFC" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIFC],[mpxlf2003_r mpxlf2003 mpxlf95_r mpxlf90 mpf95 mpf90  mpif95 mpif90  ])
fi

ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran]])])

AC_LANG(Fortran 77)
if test "X$MPIF77" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIF77],[mpxlf mpf77 mpif77])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran 77]])])
FC="$MPIFC" ;
F77="$MPIF77";
CC="$MPICC";
fi

# We leave a default language for the next checks.
dnl AC_LANG([Fortran 77])
AC_LANG([C])

dnl Now on, MPIFC should be set, as MPIF77 and MPICC

###############################################################################
# Sanity checks, although redundant (useful when debugging this configure.ac)!
###############################################################################

if test "X$MPIFC" == "X" ; then
	AC_MSG_ERROR([Problem : No MPI Fortran compiler specified nor found!])
fi

if test "X$MPICC" == "X" ; then
	AC_MSG_ERROR([Problem : No MPI C compiler specified nor found!])
fi


###############################################################################

dnl  Library oriented Autotools facilities (we don't care about this for now)
dnl AC_PROG_LIBTOOL
dnl AM_MAINTAINER_MODE
AC_PROG_RANLIB
dnl system's ranlib will be found, too

AM_INIT_AUTOMAKE
dnl Specify required version of autoconf.
AC_PREREQ(2.59)


###############################################################################
#	Compiler identification (sadly, it is necessary)
###############################################################################
mld_ainv_cv_fc=""

dnl Do we use gfortran & co ? Compiler identification.
dnl NOTE : in /autoconf/autoconf/fortran.m4 there are plenty of better tests!

PAC_CHECK_HAVE_GFORTRAN(
	[mld_ainv_cv_fc="gcc"],
)

if test "X$mld_ainv_cv_fc" == "Xgcc" ; then
PAC_HAVE_MODERN_GFORTRAN(
	[],
        [AC_MSG_ERROR([Sorry, we require GNU Fortran 4.8 or later.])]
)
fi




##############################################################################
# Choice of our compilers, needed by Make.inc
##############################################################################
if test "X$FLINK" == "X" ; then
   FLINK=${MPF77}
fi

if test "X$F90LINK" == "X" ; then
   F90LINK=${MPF90}
fi
##############################################################################
##############################################################################
#			Include flags
##############################################################################
AX_F90_MODULE_EXTENSION
AX_F90_MODULE_FLAG
MODEXT=".$ax_cv_f90_modext"
FMFLAG="${ax_cv_f90_modflag%%[ ]*}"
FIFLAG=-I

###############################################################################
###############################################################################
# NOTE :
# Missing stuff : 
# 		In the case the detected fortran compiler is ifort, icc or gcc
#		should be valid options.
#		The same for pg (Portland Group compilers).
###############################################################################

PAC_FORTRAN_HAVE_MLD2P4([AC_MSG_RESULT([yes.])],
	[AC_MSG_ERROR([no. Could not find working version of MLD2P4.])]
)


PAC_FORTRAN_MLD2P4_VERSION()

if test "x$pac_cv_mld2p4_major" == "xunknown"; then
  AC_MSG_ERROR([MLD2P4 version major "$pac_cv_mld2p4_major".])
fi
if test "x$pac_cv_mld2p4_minor" == "xunknown"; then
   AC_MSG_ERROR([MLD2P4 version minor "$pac_cv_mld2p4_minor".])
fi
if test "x$pac_cv_mld2p4_patchlevel" == "xunknown"; then
   AC_MSG_ERROR([MLD2P4 patchlevel "$pac_cv_mld2p4_patchlevel".])
fi
if (( $pac_cv_mld2p4_major < 2 )) ||
    (  (( $pac_cv_mld2p4_major == 2 )) && (( $pac_cv_mld2p4_minor < 1 )) ) ||
      ( (( $pac_cv_mld2p4_major == 2 )) && (( $pac_cv_mld2p4_minor == 1 )) && (( $pac_cv_mld2p4_patchlevel < 1 )) )
      	then			
  AC_MSG_ERROR([I need at least MLD2P4 version 2.1.1.]) 
else
  AC_MSG_NOTICE([Am configuring with MLD2P4 version $pac_cv_mld2p4_major.$pac_cv_mld2p4_minor.$pac_cv_mld2p4_patchlevel.])
fi


PAC_ARG_WITH_SAINV

PAC_FORTRAN_HAVE_SAINV([mld_ainv_ext_sainv="yes"],
	[mld_ainv_ext_sainv="no"]
)

if test "X$mld_ainv_ext_sainv" == "Xyes" ; then
   EXT_AINV_DEFS="$psblas_cv_define_prepend-DHAVE_TUMA_SAINV";
   EXT_AINV_LIBS="$SAINV_LIBS";
fi


###############################################################################
# Variable substitutions : the Make.inc.in will have these @VARIABLES@
# substituted.

AC_SUBST(MLD2P4_DIR)
AC_SUBST(INSTALL)
AC_SUBST(INSTALL_DATA)
AC_SUBST(INSTALL_DIR)
AC_SUBST(INSTALL_LIBDIR)
AC_SUBST(INSTALL_INCLUDEDIR)
AC_SUBST(INSTALL_MODULESDIR)
AC_SUBST(INSTALL_DOCSDIR)
AC_SUBST(FDEFINES)
AC_SUBST(CDEFINES)
AC_SUBST(EXT_AINV_DEFS)
AC_SUBST(EXT_AINV_LIBS)




###############################################################################
# the following files will be created by Automake

AC_CONFIG_FILES([Make.inc])
AC_OUTPUT()
#AC_OUTPUT(Make.inc Makefile)
###############################################################################

dnl Please note that brackets around variable identifiers are absolutely needed for compatibility..
AC_MSG_NOTICE([
	${PACKAGE_NAME} ${mld_ainv_cv_version} has been configured as follows:
	EXT_AINV = ${EXT_AINV_LIBS}
	If you are satisfied, run 'make' to build ${PACKAGE_NAME} and its documentation; otherwise
	type  ./configure --help=short for a complete list of configure options specific to ${PACKAGE_NAME}.
dnl	To install the program and its documentation, run 'make install' if you are root,
dnl	or run 'su -c "make install"' if you are not root.
])

###############################################################################

