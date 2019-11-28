# $Id: Makefile.am,v 1.21 2005/08/30 11:48:53 rdmp1c Exp $

# The root makefile. 

TREELIBDIR	= TreeLib
GPORTDIR	= $(TREELIBDIR)/gport

# to do: add XML support back into Unix version (cf. TreeEdit)
#EXPATPPDIR 	= $(TREELIBDIR)/expatpp
#XMLPARSEDIR = $(EXPATPPDIR)/xmlparse

NCLDIR		= ncl-2.0/src
#LIBDIRS	= -L$(NCLDIR) 
#-L$(XMLPARSEDIR)

LIBS		= @LIBS@ -L$(NCLDIR) -lncl -L$(TREELIBDIR) -ltreelib
#-lexpat
INCLUDES 	= -I$(NCLDIR) -I$(TREELIBDIR) -I$(GPORTDIR) 
#-I$(EXPATPPDIR) -I$(XMLPARSEDIR) -I$(EXPATPPDIR)/xmltok

# We want to build the tree and NEXUS parsing code as separate libraries
SUBDIRS = ncl-2.0 TreeLib

# The actual TreeView X program
bin_PROGRAMS = tv
tv_SOURCES = tview.cpp tdoc.cpp tv.cpp treeorder_dialog.cpp


# Include headers, bitmaps, icons, and IDE files in the distribution
EXTRA_DIST = tdoc.h tproject.h tv.h tview.h treeorder_dialog.h \
   bitmaps/paste.xpm bitmaps/new.xpm bitmaps/copy.xpm bitmaps/phylogram.xpm bitmaps/print_preview.xpm \
   bitmaps/slanttree.xpm bitmaps/next.xpm bitmaps/previous.xpm bitmaps/rectangletree.xpm \
   bitmaps/treeview.xpm bitmaps/open.xpm bitmaps/print.xpm bitmaps/saveas.xpm \
   bitmaps/mac/copy.xpm bitmaps/mac/phylogram.xpm bitmaps/mac/print_preview.xpm \
   bitmaps/mac/slanttree.xpm bitmaps/mac/next.xpm bitmaps/mac/previous.xpm bitmaps/mac/rectangletree.xpm \
   bitmaps/mac/open.xpm bitmaps/mac/print.xpm bitmaps/mac/saveas.xpm  \
   bitmaps/mac/internal.xpm bitmaps/mac/new.xpm bitmaps/mac/paste.xpm  \
   bitmaps/mac/zoomin.xpm bitmaps/mac/zoomout.xpm bitmaps/zoomin.xpm bitmaps/zoomout.xpm bitmaps/internal.xpm \
   tv.pbproj/project.pbxproj tv.pbproj/rpage.pbxuser \
   tv.xcode/project.pbxproj tv.xcode/rpage.pbxuser \
   bitmaps/mac/app.icns bitmaps/mac/doc.icns



