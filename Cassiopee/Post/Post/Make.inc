#
# -- Compilation des fortrans --
#
# Definis dans le Makefile :
# SOURCEF77 : liste des source f77
# SOURCEF90 : liste des source f90
# LIBNAME : nom librairie
# Repertoire de production des .o fortrans
PRODDIR   = ../../build/$(ELSAPROD)/

# -- Ne pas modifier ci-dessous --
SRC77     = $(shell ../MakeObjs.mk "$(SOURCEF77)")
SRC90     = $(shell ../MakeObjs.mk "$(SOURCEF90)")
OBJS      = $(SRC77:%.for=%.o) $(SRC90:%.f90=%.o) 
OBJS2     = $(shell ../MakeObjsPath.mk "$(OBJS)") 
LIB2      = $(shell ../MakeObjsPath.mk "lib$(LIBNAME).a")
all: link $(LIB2)

.SUFFIXES: .for .f .f90

.for.o:
	echo "Compiling fortran 77: " $<
	@(cpp -traditional -I../../../KCore/KCore -DINTEGER_E="INTEGER*4" -DREAL_E="REAL*8" $< > $(PRODDIR)/temp.f; \
	$(FC) $(FOPT0) $(FOPT1) $(FOPT2) -c $(PRODDIR)/temp.f -o $@; rm $(PRODDIR)/temp.f)

.f90.o:
	echo "Compiling fortran 90: " $<
	@(cpp -traditional -I../../../KCore/KCore -DINTEGER_E="INTEGER*4" -DREAL_E="REAL*8" $< > $(PRODDIR)/temp.f90; \
        $(F90) $(FOPT0) $(FOPT1) $(FOPT2) -c $(PRODDIR)/temp.f90 -o $@; rm $(PRODDIR)/temp.f90)

link:	$(SRC77) $(SRC90)	
	@if [ -d ../../build/ ]; then :; else mkdir ../../build; fi
	@if [ -d $(PRODDIR) ]; then :; else mkdir $(PRODDIR); fi
	@for i in $(SRC77); do \
	../MakeLink.mk "`pwd`" $$i $(PRODDIR);\
	done
	@for i in $(SRC90); do \
	../MakeLink.mk "`pwd`" $$i $(PRODDIR);\
	done

$(LIB2) : $(OBJS2)
	@echo Making $(DIRECTORY) library $(LIB2)
	@ar clr  $@ $(OBJS2) ; ranlib $(LIB2)
	@echo done
