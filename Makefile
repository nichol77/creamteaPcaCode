#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= 
SYSLIBS         = 

SIMPLE_GEANT_SIM_DIR=/home/rjn/creamtea/simpleGeantSim


#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) -I $(SIMPLE_GEANT_SIM_DIR)/include
LDFLAGS      += -g $(ROOTLDFLAGS) 

LIBS          = $(ROOTLIBS) -lMathMore -lMinuit $(SYSLIBS) 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Now the bits we're actually compiling
ROOT_LIBRARY = libMagicDisplay.${DLLSUF}
LIB_OBJS =  AnitaCanvasMaker.o WaveformGraph.o MagicDisplay.o MagicDisplayConventions.o AnitaRFCanvasMaker.o MagicControlPanel.o FFTGraph.o magicDict.o
CLASS_HEADERS = AnitaCanvasMaker.h AnitaRFCanvasMaker.h WaveformGraph.h MagicDisplay.h MagicDisplayConventions.h MagicControlPanel.h FFTGraph.h


all : makePCAFile


#The library
makePCAFile : makePCAFile.$(OBJSUF) 
	@echo "Linking $@ ..."
	$(LD) $(LDFLAGS) $< $(LIBS) $(OutPutOpt) $@
	@echo "$@ done"

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


magicDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h

install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
	cp $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
else
	cp $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	cp  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
#############################################################################



