include ../../Make.inc
LIBDIR=../../lib
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
PSBMODDIR=$(PSBLASDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLDINCDIR=$(MLDDIR)/include
MLDMODDIR=$(MLDDIR)/modules
HERE=..

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(LIBDIR) $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) \
 $(FMFLAG)$(MLDMODDIR)  $(FMFLAG)$(PSBMODDIR) 


F90OBJS=mld_drwclip.o \
mld_d_sp_drop.o \
mld_d_sparsify.o \
mld_d_invk_solver_bld.o \
mld_d_invk_solver_check.o \
mld_d_invk_solver_descr.o \
mld_d_invk_solver_cseti.o \
mld_d_invk_solver_clone.o \
mld_d_invk_solver_clone_settings.o \
mld_d_invt_solver_bld.o \
mld_d_invt_solver_check.o \
mld_d_invt_solver_descr.o \
mld_d_invt_solver_cseti.o \
mld_d_invt_solver_csetr.o \
mld_d_invt_solver_clone.o \
mld_d_invt_solver_clone_settings.o \
mld_d_ainv_solver_bld.o \
mld_d_ainv_solver_check.o \
mld_d_ainv_solver_descr.o \
mld_d_ainv_solver_csetc.o \
mld_d_ainv_solver_cseti.o \
mld_d_ainv_solver_csetr.o \
mld_d_ainv_solver_clone.o \
mld_d_ainv_solver_clone_settings.o \
mld_d_ainv_bld.o \
mld_d_invk_copyin.o \
mld_d_invk_copyout.o \
mld_d_invk_inv.o \
mld_d_invk_bld.o \
mld_d_invt_copyin.o \
mld_d_invt_copyout.o \
mld_d_invt_inv.o \
mld_d_invt_bld.o \
mld_dsparse_invk.o \
mld_dsparse_invt.o \
mld_dsparse_tuma_sainv.o \
mld_dsparse_tuma_lainv.o \
mld_dsparse_biconjg_llk.o \
mld_dsparse_biconjg_llk_noth.o \
mld_dsparse_biconjg_mlk.o \
mld_dsparse_biconjg_s_ft_llk.o \
mld_dsparse_biconjg_s_llk.o \
mld_d_base_ainv_update_a.o \
mld_d_base_ainv_solver_apply.o \
mld_d_base_ainv_solver_apply_vect.o \
mld_d_base_ainv_solver_cnv.o \
mld_d_base_ainv_solver_dmp.o \
mld_d_base_ainv_solver_clear_data.o \
mld_d_base_ainv_solver_free.o


#mld_d_biconjg_impl.o


OBJS=$(F90OBJS) 

LIBNAME=libainv_prec.a

lib: $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)


veryclean: clean

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

