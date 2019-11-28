# source directories
SRC_DIR		= ./

# compiler switches 
CC              = gcc
CXX             = g++
CXXFLAGS        = -O4
LOADLIBES       = -lm
CLINKER         = g++

# target macros
EXECS		= ncltest basiccmdline makedoc
NCLOBJS		= allelesblock.o assumptionsblock.o charactersblock.o \
   datablock.o discretedatum.o discretematrix.o distancedatum.o \
   distancesblock.o nexus.o nexusblock.o nexustoken.o setreader.o \
   nxsstring.o nxsdate.o taxablock.o treesblock.o xnexus.o 
NCLTESTOBJS	= ncltest.o
BASICCMDLINEOBJS = basiccmdline.o
MAKEDOCOBJS	= makedoc.o

# implicit construction rule
%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

# targets
default : ncltest

all : $(EXECS)

ncltest : $(NCLOBJS) $(NCLTESTOBJS)
	$(CLINKER) -o ncltest $(NCLOBJS) $(NCLTESTOBJS) -lm

basiccmdline : $(NCLOBJS) $(BASICCMDLINEOBJS)
	$(CLINKER) -o basiccmdline $(NCLOBJS) $(BASICCMDLINEOBJS) -lm

makedoc : $(MAKEDOCOBJS)
	$(CLINKER) -o makedoc $(MAKEDOCOBJS) -lm
	
libncl.a : $(NCLOBJS)
	ar rv $@ $^
	ranlib $@

clean : FORCE
	rm -f *.o
	rm -f basiccmdline
	rm -f ncltest
	rm -f makedoc
	rm -f libncl.a

FORCE :

# object files and dependencies
allelesblock.o: allelesblock.cpp nexusdefs.h nxsstring.h xnexus.h \
 nexustoken.h nexus.h setreader.h taxablock.h discretedatum.h \
 discretematrix.h charactersblock.h allelesblock.h assumptionsblock.h
assumptionsblock.o: assumptionsblock.cpp nexusdefs.h nxsstring.h \
 xnexus.h nexustoken.h nexus.h setreader.h taxablock.h discretedatum.h \
 discretematrix.h charactersblock.h assumptionsblock.h
charactersblock.o: charactersblock.cpp nexusdefs.h nxsstring.h \
 xnexus.h nexustoken.h nexus.h setreader.h taxablock.h discretedatum.h \
 discretematrix.h assumptionsblock.h charactersblock.h
datablock.o: datablock.cpp nexusdefs.h nxsstring.h discretedatum.h \
 discretematrix.h nexustoken.h nexus.h taxablock.h charactersblock.h \
 datablock.h
discretedatum.o: discretedatum.cpp nexusdefs.h nxsstring.h \
 discretedatum.h
discretematrix.o: discretematrix.cpp nexusdefs.h nxsstring.h \
 discretedatum.h discretematrix.h
distancedatum.o: distancedatum.cpp distancedatum.h
distancesblock.o: distancesblock.cpp nexusdefs.h nxsstring.h xnexus.h \
 nexustoken.h nexus.h taxablock.h distancedatum.h distancesblock.h
emptyblock.o: emptyblock.cpp nexusdefs.h nxsstring.h xnexus.h \
 nexustoken.h nexus.h emptyblock.h
makedoc.o: makedoc.cpp
ncltest.o: ncltest.cpp nexusdefs.h nxsstring.h nexustoken.h nexus.h \
 taxablock.h treesblock.h discretedatum.h discretematrix.h \
 charactersblock.h allelesblock.h assumptionsblock.h datablock.h \
 distancedatum.h distancesblock.h
nexus.o: nexus.cpp nexusdefs.h nxsstring.h xnexus.h nexustoken.h \
 nexus.h
nexusblock.o: nexusblock.cpp nexusdefs.h nxsstring.h nexustoken.h \
 nexus.h
nexustoken.o: nexustoken.cpp nexusdefs.h nxsstring.h xnexus.h \
 nexustoken.h
nxsstring.o: nxsstring.cpp nxsstring.h
nxsdate.o: nxsdate.cpp nxsdate.h
setreader.o: setreader.cpp nexusdefs.h nxsstring.h xnexus.h \
 nexustoken.h nexus.h setreader.h
taxablock.o: taxablock.cpp nexusdefs.h nxsstring.h xnexus.h \
 nexustoken.h nexus.h taxablock.h
treesblock.o: treesblock.cpp nexusdefs.h nxsstring.h xnexus.h \
 nexustoken.h nexus.h taxablock.h treesblock.h
xnexus.o: xnexus.cpp nexusdefs.h nxsstring.h xnexus.h

