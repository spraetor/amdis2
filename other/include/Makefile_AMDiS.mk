PARMETIS_DIR = $(AMDIS_DIR)/lib/ParMetis-3.1

# ============================================================================
# ===== includes pathes ======================================================
# ============================================================================

AMDIS_INCLUDE = -I$(AMDIS_DIR)/include -I$(AMDIS_DIR)/include/compositeFEM -I$(AMDIS_DIR)/include/reinit -I$(AMDIS_DIR)/include/mtl4
PARMETIS_INCLUDE = -I$(PARMETIS_DIR)
PARALLEL_INCLUDE = -I$(AMDIS_DIR)/src/parallel

INCLUDES = -I. $(AMDIS_INCLUDE) $(PARALLEL_INCLUDE)

# ============================================================================
# ===== libraries ============================================================
# ============================================================================

AMDIS_LIB = -L$(AMDIS_DIR)/lib -lamdis -lcompositeFEM -lreinit
PNG_LIB = 
LIBS =

UMFPACK_LIB = 
ifeq ($(strip $(USE_SERVER)), mars)
	UMFPACK_LIB += -lmkl -lumfpack -lamd
else
	ifeq ($(strip $(USE_SERVER)), themisto)
		UMFPACK_LIB += $(MKL_LIB) -lmkl -lguide -lpthread -lumfpack -lamd
	else	
		ifeq ($(strip $(USE_SERVER)), deimos)
#			LIBS += -L/fastfs/wir/local/lib
			MPI_DIR = /licsoft/libraries/openmpi/1.2.6/64bit
			UMFPACK_LIB += -lumfpack -lamd -L/licsoft/libraries/goto -lgoto -lpthread
		else
			UMFPACK_LIB += -lblas -lumfpack -lamd
		endif
	endif
endif

PARMETIS_LIB = -L$(PARMETIS_DIR) -lparmetis -lmetis
ZOLTAN_LIB = -L$(AMDIS_DIR)/lib/zoltan_build/lib -lzoltan

LIBS += $(AMDIS_LIB) $(PNG_LIB)
LIBS += -lboost_iostreams -lboost_filesystem -lboost_system -lboost_date_time

ifeq ($(strip $(USE_UMFPACK)), 1)
	LIBS += $(UMFPACK_LIB)
	CPPFLAGS += -DHAVE_UMFPACK -DMTL_HAS_UMFPACK
	INCLUDES += -I$(AMDIS_DIR)/include/umfpack -I$(AMDIS_DIR)/include/ufconfig -I$(AMDIS_DIR)/include/amd
endif

ifeq ($(strip $(USE_MKL)), 1)
	ifeq ($(strip $(USE_SERVER)), themisto)
		LIBS += $(MKL_LIB) 
	else
		LIBS += -L$(MKL_LIB) 
	endif

	LIBS += -lmkl -lmkl_solver -lguide -lpthread
endif

# ============================================================================
# ===== parallel or sequential ? =============================================
# ============================================================================

ifeq ($(strip $(USE_PARALLEL_AMDIS)), 1)
	include ${PETSC_DIR}/conf/variables
	CFLAGS = ${PETSC_CC_INCLUDES}
	FFLAGS = ${PETSC_FC_INCLUDES}

	ifeq ($(strip $(USE_SERVER)), mars)
		ifeq ($(strip $(USE_COMPILER)), gcc)
			COMPILE = g++
		else
			COMPILE = icpc
		endif
	else
		COMPILE = $(MPI_DIR)/bin/mpiCC
	endif

	CPPFLAGS += -DHAVE_PARALLEL_DOMAIN_AMDIS
	INCLUDES += $(PETSC_CC_INCLUDES)
	LIBS += $(PARMETIS_LIB) -lmpi $(PETSC_LIB)
	LIBS += $(ZOLTAN_LIB)
else
	ifeq ($(strip $(USE_COMPILER)), gcc)
		COMPILE = g++
	else
		COMPILE = icpc
	endif
endif

ifeq ($(strip $(USE_MARMOT)), 1)
	COMPILE = marmotcxx
endif

# ============================================================================
# ===== compile flags ========================================================
# ============================================================================

ifeq ($(strip $(DEBUG)), 0)
       CPPFLAGS += -g -O3 -DNDEBUG
else
       CPPFLAGS += -g -O0 -DDEBUG=1
endif

# ============================================================================
# ===== object directory =====================================================
# ============================================================================
ifeq ($(OBJDIR),)
	OBJDIR = .
endif
# ============================================================================
# ===== libtool linking ======================================================
# ============================================================================

LIBTOOL = $(AMDIS_DIR)/libtool

ifeq ($(strip $(USE_PARALLEL_AMDIS)), 1)
	LINK = $(LIBTOOL) --tag=mpiCC --mode=link $(COMPILE)
else
	LINK = $(LIBTOOL) --mode=link $(COMPILE)
endif

# ============================================================================
# ===== rules ================================================================
# ============================================================================

all : 
	make $(PROGRAMS)

clean: 
	-rm -rf $(OBJDIR)/*.o
	-rm -rf $(PROGRAMS)

.cc.o: $*.cc
	$(COMPILE) $(DEFS) $(INCLUDES) $(CPPFLAGS) -c -o $(OBJDIR)/$*.o $^ 



