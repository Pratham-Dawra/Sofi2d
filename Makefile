.SUFFIXES:
.SECONDARY:
.PHONY: all clean distclean help install uninstall force doc

MAKEFLAGS += -rR
SHELL     := /bin/bash

###############################################################################
# The SOFI top-level directory as checked out from git
# Example: /path/to/sofi/ (absolute) or ../ (relative to build directory)
# Note: $(BASEDIR)/src should contain the SOFI source code
###############################################################################
BASEDIR  := ..

###############################################################################
# Installation directory: 
# Example: /path/to/inst/bin (absolute) or ../bin (relative to build directory)
###############################################################################
INSTDIR := ../bin

###############################################################################
# Name of executable
###############################################################################
FDMOD := sofi2D
MERGE := snapmerge

###############################################################################
# Name of documentation (PDF)
###############################################################################
DOC := sofi2D_manual.pdf

###############################################################################
# Compiler / MPI settings
###############################################################################
### CC = C compiler to use
### CPPFLAGS = common preprocessor flags
### EXTRA_CPPFLAGS = extra preprocessor flags for sofi2D
### CFLAGS = common compiler flags
### CDEPEND = compiler flags to produce dependency files for each object
### LD = linker to use
### LDFLAGS = common linker flags
### EXTRA_LDFLAGS = extra linker flags for sofi2D
###
### In case you have trouble with ANSI color escape sequences in your terminal,
### you can switch off color coding for logging by removing the -DLOG_COLOR
### flag.
###
### Use -DEBUG to compile code with lots of NaN checks enabled. Should only be
### used for low-level debugging of instability problems. It produces tons of
### output and slows down the program significantly.

###############################################################################
### GNU compiler with OpenMPI (standard)
CC             := mpicc
CXX            := mpicxx
CPPFLAGS        = -D_GNU_SOURCE -DLOG_COLOR
EXTRA_CPPFLAGS := -DUSE_MPI
CFLAGS          = -O3 -Wall -Wextra -g -fstack-protector -std=c99 -ftree-vectorize
CXXFLAGS        = -O3 -Wall -Wextra -g -fstack-protector -ftree-vectorize
CDEPEND         = -MMD -MP -MF .dep-$(subst /,-,$@)
LD             := $(CC)
LDFLAGS         = -lm
EXTRA_LDFLAGS  := -lmpi -lfftw3 -lstdc++
#EXTRA_LDFLAGS  := -lmpi -lfftw3
#EXTRA_LDFLAGS  := -lmpi -lstdc++

###############################################################################
### Intel compiler (OneAPI) with Intel MPI
### => GPI: module load compiler/latest mpi/latest
#CC             := mpicc -cc=icx
#CXX            := 
#CPPFLAGS        = -D_GNU_SOURCE -DLOG_COLOR
#EXTRA_CPPFLAGS := -DLOG_MPI
#CFLAGS          = -O3 -Wall -Wextra -g -fstack-protector -std=c99 --intel
#CXXFLAGS        = -O3 -Wall -Wextra -g -fstack-protector --intel
#CDEPEND         = -MMD -MP -MF .dep-$(subst /,-,$@)
#LD             := $(CC)
#LDFLAGS         = -static-intel -lm 
#EXTRA_LDFLAGS  := -lmpi -lstdc++

###############################################################################
### clang (LLVM) with OpenMPI
#CC             := clang
#CXX            := clang++
#CPPFLAGS        = -D_GNU_SOURCE -DLOG_COLOR -I/usr/lib64/mpi/gcc/openmpi3/include
#EXTRA_CPPFLAGS := -DLOG_MPI
#CFLAGS          = -O3 -Wall -Wextra -g -fstack-protector -std=c99
#CXXFLAGS        = -O3 -Wall -Wextra -g -fstack-protector
#CDEPEND         = -MMD -MP -MF .dep-$(subst /,-,$@)
#LD             := $(CC)
#LDFLAGS         = -lm 
#EXTRA_LDFLAGS  := -L/usr/lib64/mpi/gcc/openmpi3/lib64 -lmpi -lstdc++

###############################################################################
###############################################################################
# Usually, there is no need to change anything below this line 
###############################################################################
###############################################################################

CURRENT   := $(PWD)
SRCTREE   := $(strip $(realpath $(BASEDIR)/src))
CONTREE   := $(strip $(realpath $(BASEDIR)/contrib))
AFFTREE   := $(strip $(realpath $(CONTREE)/Seitosh/libaff))
SEITREE   := $(strip $(realpath $(CONTREE)/Seitosh/libseife))
FOUTREE   := $(strip $(realpath $(CONTREE)/Seitosh/libfourier))
INVTREE   := $(strip $(realpath $(CONTREE)/Seitosh/libstfinv))
FDMOD_O   := $(strip ./fdmod_obj)
MERGE_O   := $(strip ./merge_obj)
AFF_O     := $(strip ./libaff_obj)
FOU_O     := $(strip ./libfourier_obj)
SEI_O     := $(strip ./libseife_obj)
INV_O     := $(strip ./libstfinv_obj)
LIBAFF    := libaff.a
LIBSEI    := libcseife.a
LATEX_O   := $(strip ./latex)
dummy     := $(shell mkdir -p $(FDMOD_O))
dummy     := $(shell mkdir -p $(MERGE_O))
dummy     := $(shell mkdir -p $(AFF_O))
dummy     := $(shell mkdir -p $(FOU_O))
dummy     := $(shell mkdir -p $(SEI_O))
dummy     := $(shell mkdir -p $(INV_O))
dummy     := $(shell mkdir -p $(LATEX_O))
dummy     := $(shell mkdir -p $(strip $(INSTDIR)))
FDMOD_EXE := $(strip $(FDMOD))
MERGE_EXE := $(strip $(MERGE))
INSTALL   := install -D -p --backup=simple --suffix=.bck
AR        := ar
ARFLAGS   := cr
RM        := rm -f --preserve-root
LATEX     := pdflatex
BIBER     := biber
DEPS      := $(wildcard .dep*)
COMMIT    := $(strip $(shell cd $(BASEDIR) && git log --format="%h %as" 2>/dev/null | head -n1))
ifeq ($(COMMIT),)
  COMMIT  := unknown
endif
LOCHANGE  := $(strip $(shell cd $(BASEDIR) ; git diff-index --quiet HEAD -- 2>/dev/null ; echo $$?))
ifeq ($(LOCHANGE),1)
  COMMIT  += + local
endif
CPPFLAGS  += -DLOG_COMMIT="\"$(COMMIT)\""
ifdef FFTW_INC
   CPPFLAGS += -I$(FFTW_INC)
endif
LIBAFF    := libaff.a
LIBFOU    := libfourierxx.a
LIBSEI    := libcseife.a
LIBINV    := libstfinv.a
LOCALLIBS  = -L$(PWD) $(addprefix -l,$(subst lib,,$(basename $(LIBAFF) $(LIBSEI) $(LIBFOU) $(LIBINV)))) \
             -lfftw3

OS := $(shell uname)
ifeq ($(OS),Darwin)
  INSTALL := install -p
  RM      := rm -f 
endif

OS := $(shell uname)
ifeq ($(OS),Darwin)
  INSTALL := install -p
  RM      := rm -f 
endif

DOC_MAIN := main
DOC_FILES := $(wildcard $(strip $(realpath $(BASEDIR)))/doc/guide_sofi2D/latex/*.tex) \
             $(wildcard $(strip $(realpath $(BASEDIR)))/doc/guide_sofi2D/latex/figures/*) \
             $(strip $(realpath $(BASEDIR)))/doc/guide_sofi2D/latex/sofi2D.json \
             $(strip $(realpath $(BASEDIR)))/doc/guide_sofi2D/latex/biblio.bib 

MERGE_FILES := snapmerge.c \
               holbergcoeff.c \
               json_parser.c \
               logging.c \
               merge.c \
               read_par_json.c \
               read_par_json_fwi.c \
               readdsk.c \
               sources.c \
               util.c \
               writedsk.c

FDMOD_FILES := sofi2D.c \
               abs_update.c \
               absorb.c \
               acq_read.c \
               apply_workflow.c \
               av_mue.c \
               av_rho.c \
               av_tau.c \
               calc_conv.c \
               calc_envelope.c \
               calc_res.c \
               catseis.c \
               check_fs.c \
               checkfd.c \
               cpml_update.c \
               debug_buffers.c \
               eprecond.c \
               eqsource.c \
               exchange_par.c \
               exchange_s.c \
               exchange_v.c \
               filter_frequencies.c \
               freemem.c \
               freemem_model.c \
               freemem_wavefield.c \
               hh_elastic.c \
               hh_elastic_TTI.c \
               hh_elastic_VTI.c \
               hh_visco.c \
               hh_visco_TTI.c \
               hh_visco_VTI.c \
               holbergcoeff.c \
               initfd.c \
               init_grad.c \
               init_inv.c \
               initmem.c \
               initmem_fwi.c \
               initmem_model.c \
               initmem_wavefield.c \
               initproc.c \
               initsrc.c \
               inseis.c \
               inversion.c \
               json_parser.c \
               kiss_fftr.c\
               kiss_fft.c\
               lbfgs.c \
               logging.c \
               matcopy.c \
               matcopy_elastic.c \
               mergemod.c \
               outseis_glob.c \
               PML_pro.c \
               prepare_update_s_ac.c \
               prepare_update_s_el.c \
               prepare_update_s_el_4.c \
               prepare_update_s_visc.c \
               prepare_update_s_visc_4.c \
               prepare_update_s_vti.c \
               prepare_update_s_tti.c \
               prepmod.c \
               psource.c \
               read_par_json.c \
               read_par_json_fwi.c \
               read_srcsig.c \
               read_su.c \
               read_workflow.c \
               readdsk.c \
               readmod.c \
               readmod_acoustic.c \
               readmod_visco.c \
               readmod_elastic.c \
               readmod_elastic_vti.c \
               readmod_acoustic_vti.c \
               readmod_acoustic_tti.c \
               readmod_elastic_tti.c \
               readmod_visco_vti.c \
               readmod_visco_tti.c \
               receiver.c \
               saveseis.c \
               saveseis_fwi.c \
               saveseis_glob.c \
               scan_topo.c \
               seismo_ssg.c \
               seismo_shift.c\
               snap_ssg.c \
               snap_store.c \
               sources.c \
               splitrec.c \
               splitsrc.c \
               stfi.c\
               stfi_apply.c\
               stfi_merge.c\
               stfi_calc.c\
               su_gather.c \
               su_struct.c \
               subgrid_bounds.c \
               surface.c \
               surface_elastic.c \
               timedomain_filt.c \
               time_loop.c \
               timeloop.c\
               update_s_elastic_abs.c \
               update_s_acoustic_abs.c \
               update_s_elastic_VTI_abs.c \
               update_s_elastic_abs_4.c \
               update_s_elastic_interior.c \
               update_s_acoustic_interior.c \
               update_s_elastic_VTI_interior.c \
               update_s_acoustic_VTI_interior.c \
               update_s_elastic_TTI_interior.c \
               update_s_elastic_interior_4.c \
               update_s_elastic_PML.c \
               update_s_acoustic_PML.c \
               update_s_elastic_VTI_PML.c \
               update_s_elastic_TTI_PML.c \
               update_s_elastic_TTI_abs.c \
               update_s_elastic_PML_4.c \
               update_s_visc_abs.c \
               update_s_visc_interior.c \
               update_s_visc_VTI_interior.c \
               update_s_visc_TTI_interior.c \
               update_s_visc_PML.c \
               update_s_visc_VTI_PML.c \
               update_s_visc_TTI_PML.c \
               update_s_visc_VTI_abs.c \
               update_s_visc_TTI_abs.c \
               update_s_visc_abs_4.c \
               update_s_visc_interior_4.c \
               update_s_visc_PML_4.c \
               update_v_abs.c \
               update_v_abs_4.c \
               update_v_interior.c \
               update_v_interior_4.c \
               update_v_PML.c \
               update_v_PML_4.c \
               util.c \
               v_derivatives.c \
               wavefield_update_s_el.c \
               wavefield_update_s_ac.c \
               wavefield_update_s_el_vti.c \
               wavefield_update_s_ac_vti.c \
               wavefield_update_s_el_tti.c \
               wavefield_update_s_el_4.c \
               wavefield_update_s_visc.c \
               wavefield_update_s_visc_VTI.c \
               wavefield_update_s_visc_TTI.c \
               wavefield_update_s_visc_4.c \
               wavefield_update_v.c \
               wavefield_update_v_ac.c \
               wavefield_update_v_4.c \
               wavelet.c \
               write_par.c \
               write_par_fwi.c \
               write_su.c \
               writedsk.c \
               writemod.c \
               zero_wavefield.c

AFF_FILES := error.cc \
             seriesstepper.cc \
             strided.cc \
             stridedstepper.cc

SEI_FILES := cseife.c \
             cseife_deriv.c \
             cseife_gauss.c \
             cseife_rekfl.c \
             cseife_rfk.c \
             cseife_tides.c

FOU_FILES := error.cc \
             fcommand.cc \
             fftwaffar.cc \
             fftwaff.cc \
             filters.cc \
             polesnzeroes.cc

INV_FILES := error.cc \
             fouriertools.cc \
             parameterhandler.cc \
             stfinvany.cc \
             stfinvany_description_usage.cc \
             stfinvbase.cc \
             stfinvbase_description_usage.cc \
             stfinvbase_summary_usage.cc \
             stfinv.cc \
             stfinv_description_usage.cc \
             stfinvfdleastsquares.cc \
             stfinvfdleastsquares_description_usage.cc \
             stfinvfdleastsquares_summary_usage.cc \
             stfinvfinitecausal.cc \
             stfinvfixedstf.cc \
             stfinvfourier.cc \
             stfinvfourier_description_usage.cc \
             stfinvfourier_summary_usage.cc \
             stfinvidentity.cc \
             stfinvidentity_description_usage.cc \
             stfinvidentity_summary_usage.cc \
             stfinvnormalize.cc \
             stfinv_summary_usage.cc \
             tools.cc

MERGE_OBJ := $(addprefix $(MERGE_O)/, $(patsubst %.c, %.o, $(MERGE_FILES)))
FDMOD_OBJ := $(addprefix $(FDMOD_O)/, $(patsubst %.c, %.o, $(FDMOD_FILES)))
AFF_OBJ   := $(addprefix $(AFF_O)/, $(patsubst %.cc, %.o, $(AFF_FILES)))
SEI_OBJ   := $(addprefix $(SEI_O)/, $(patsubst %.c, %.o, $(SEI_FILES)))
FOU_OBJ   := $(addprefix $(FOU_O)/, $(patsubst %.cc, %.o, $(FOU_FILES)))
INV_OBJ   := $(addprefix $(INV_O)/, $(patsubst %.cc, %.o, $(INV_FILES)))

ifdef V
  ifeq ("$(origin V)", "command line")
    BUILD_VERBOSE := $(V)
  endif
endif
ifndef BUILD_VERBOSE
  BUILD_VERBOSE := 0
endif
ifeq ($(BUILD_VERBOSE),1)
  quiet := 
  Q     := 
else
  quiet := quiet_
  Q     := @
endif
ifneq ($(findstring s,$(MAKEFLAGS)),)
  quiet := silent_
endif

help:
	@echo
	@echo -e "\033[4mUsage:\033[0m"
	@echo "   $(MAKE) [V=0|1] <target>"
	@echo
	@echo -e "\033[4mTargets:\033[0m"
	@echo "   all       : compile main programs"
	@echo "   force     : force recompilation of main programs"
	@echo "   install   : install main programs"
	@echo "   doc       : compile documentation (requires LaTeX)"
	@echo "   clean     : remove object files and executables"
	@echo "   distclean : clean up entire build directory"
	@echo "   uninstall : remove previously installed programs"
	@echo 
	@echo -e "Main programs:"
	@echo "   $(FDMOD_EXE), $(MERGE_EXE)"
	@echo -e "Installation directory:"
	@echo "   $(strip $(INSTDIR))"
	@echo -e "Compiler / linker and flags:"
	@echo "   CC: $(CC)"
	@echo "   CXX: $(CXX)"
	@echo "   CPPFLAGS: $(CPPFLAGS) [$(EXTRA_CPPFLAGS)]"
	@echo "   CFLAGS: $(CFLAGS)"
	@echo "   CXXFLAGS: $(CXXFLAGS)"
	@echo "   LD: $(LD)"
	@echo "   LOCALLIB: $(LOCALLIB)"
	@echo "   LDFLAGS: $(LDFLAGS) [$(EXTRA_LDFLAGS)]"
	@echo 
	@echo "Use V=1 to see full command being executed."
	@echo 

install: all
	$(call cmd,install,$(FDMOD_EXE))
	$(call cmd,install,$(MERGE_EXE))

clean: 
	$(call cmd,rmdir,$(FDMOD_O))
	$(call cmd,rmdir,$(MERGE_O))
	$(call cmd,rmdir,$(AFF_O))
	$(call cmd,rmdir,$(SEI_O))
	$(call cmd,rmdir,$(FOU_O))
	$(call cmd,rmdir,$(INV_O))
	@echo -e "   \e[31m[RM]\e[0m    executables $(FDMOD_EXE) $(MERGE_EXE)"
	-$(Q)$(RM) $(FDMOD_EXE) $(MERGE_EXE)
	@echo -e "   \e[31m[RM]\e[0m    archives $(LIBAFF) $(LIBSEI) $(LIBFOU) $(LIBINV)"
	-$(Q)$(RM) $(LIBAFF) $(LIBSEI) $(LIBFOU) $(LIBINV)
	@$(RM) *~

distclean: clean
	@echo -e "   \e[31m[RM]\e[0m    dependency files .dep-*"
	-$(Q)find . -iname ".dep-*" | xargs $(RM)
	@echo -e "   \e[31m[RM]\e[0m    $(DOC)"
	-$(Q)$(RM) $(DOC)
	$(call cmd,rmdir,$(LATEX_O))
	-$(Q)$(RM) *~

uninstall:
	@echo -e "   \e[31m[RM]\e[0m    main programs in $(strip $(INSTDIR))"
	-$(Q)$(RM) $(strip $(INSTDIR))/$(FDMOD_EXE) $(strip $(INSTDIR))/$(MERGE_EXE) 

force:
	@$(MAKE) -B all

all: $(FDMOD_EXE) $(MERGE_EXE)

doc: $(DOC)

$(DOC): $(DOC_FILES)
	@cd $(strip $(realpath $(BASEDIR)))/doc/guide_sofi2D/latex && \
	$(LATEX) -output-directory $(CURRENT)/$(LATEX_O) $(DOC_MAIN) && \
	$(BIBER) --output-directory $(CURRENT)/$(LATEX_O) $(DOC_MAIN) && \
	$(LATEX) -output-directory $(CURRENT)/$(LATEX_O) $(DOC_MAIN) && \
	$(LATEX) -output-directory $(CURRENT)/$(LATEX_O) $(DOC_MAIN) && \
	cd $(CURRENT) && mv $(LATEX_O)/$(DOC_MAIN).pdf $(DOC) && \
	echo -e "   \e[33m[DOC]\e[0m   $(DOC)"
	@cd $(CURRENT)

$(MERGE_EXE): $(MERGE_OBJ)
	$(call cmd,link)

$(MERGE_O)/%.o: $(SRCTREE)/%.c
	$(call cmd,compile_cc)

$(LIBAFF): $(AFF_OBJ)
	$(call cmd,ar)

$(LIBSEI): $(SEI_OBJ)
	$(call cmd,ar)

$(LIBFOU): $(FOU_OBJ)
	$(call cmd,ar)

$(LIBINV): $(INV_OBJ)
	$(call cmd,ar)

LDFLAGS.$(FDMOD_EXE) = $(LOCALLIBS)
$(FDMOD_EXE): LDFLAGS += $(EXTRA_LDFLAGS)
$(FDMOD_EXE): $(LIBAFF) $(LIBSEI) $(LIBFOU) $(LIBINV) $(FDMOD_OBJ)
	$(call cmd,link)

$(FDMOD_O)/%.o: CPPFLAGS += $(EXTRA_CPPFLAGS) -I$(SEITREE) -I$(INVTREE)
$(FDMOD_O)/%.o: $(SRCTREE)/%.c
	$(call cmd,compile_cc)

$(AFF_O)/%.o: $(AFFTREE)/%.cc
	$(call cmd,compile_cxx)

$(SEI_O)/%.o: $(SEITREE)/%.c
	$(call cmd,compile_cc)

$(FOU_O)/%.o: CPPFLAGS += -I$(AFFTREE)
$(FOU_O)/%.o: $(FOUTREE)/%.cc
	$(call cmd,compile_cxx)

$(INV_O)/%.o: CPPFLAGS += -I$(AFFTREE) -I$(FOUTREE)
$(INV_O)/%.o: $(INVTREE)/%.cc
	$(call cmd,compile_cxx)

quiet_cmd_compile_cc = \e[32m[CC]\e[0m    $@
      cmd_compile_cc = $(CC) -o$@ $(CDEPEND) $(CFLAGS.$@) $(CPPFLAGS) -I$(SRCTREE) $(CFLAGS) -c $<

quiet_cmd_compile_cxx = \e[32m[CXX]\e[0m   $@
      cmd_compile_cxx = $(CXX) -o$@ $(CDEPEND) $(CXXFLAGS.$@) $(CPPFLAGS) $(CXXFLAGS) -c $<

quiet_cmd_link = \e[94m[LD]\e[0m    $@
      cmd_link = $(LD) $(CFLAGS) -o$@ $(filter-out $(LIBAFF) $(LIBSEI) $(LIBFOU) $(LIBINV),$^) \
                 $(LDFLAGS.$@) $(LDFLAGS)

quiet_cmd_ar = \e[33m[AR]\e[0m    $@
      cmd_ar = $(AR) $(ARFLAGS) $(abspath $@) $^ && ranlib $(abspath $@)

quiet_cmd_rmdir = \e[31m[RM]\e[0m    directory $(1)
      cmd_rmdir = $(RM) -r $(1)

quiet_cmd_install = \e[36m[INST]\e[0m  $(1) -> $(strip $(INSTDIR))
      cmd_install = $(INSTALL) $(1) $(strip $(INSTDIR))

cmd = @$(if $($(quiet)cmd_$(1)),echo -e '   $(call $(quiet)cmd_$(1),$(2))' &&) $(call cmd_$(1),$(2))

ifneq ($(DEPS),)
  -include $(DEPS)
endif
