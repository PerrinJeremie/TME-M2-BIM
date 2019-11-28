
noinst_LIBRARIES = libtreelib.a
if ENT
INCLUDES = -Igport -I../graph
libtreelib_a_SOURCES = nodeiterator.h Parse.cpp Parse.h \
   profile.h tokeniser.cpp tokeniser.h TreeLib.cpp TreeLib.h \
   treedrawer.cpp treedrawer.h treereader.cpp treereader.h \
   treewriter.cpp treewriter.h treeorder.cpp treeorder.h   \
   ntree.cpp ntree.h  mast.cpp mast.h lcaquery.cpp lcaquery.h quartet.cpp quartet.h  \
   gport/gport.cpp
else
INCLUDES = -Igport 

libtreelib_a_SOURCES = nodeiterator.h Parse.cpp Parse.h \
   profile.h tokeniser.cpp tokeniser.h TreeLib.cpp TreeLib.h \
   treedrawer.cpp treedrawer.h treereader.cpp treereader.h \
   treewriter.cpp treewriter.h treeorder.cpp treeorder.h
endif

EXTRA_DIST = gport/gport.h gport/gdefs.h gport/gport.cpp
