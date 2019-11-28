/*
 * TreeView
 * A phylogenetic tree viewer.
 * Copyright (C) 2001 Roderic D. M. Page <r.page@bio.gla.ac.uk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307, USA.
 */
 
 // $Id: tv.cpp,v 1.43 2007/03/09 13:58:19 rdmp1c Exp $

#ifdef __GNUG__
//	#pragma implementation "docview.h"
#endif

// hack because gcc 2.95 on my Solaris box can't handle the VERSION flag passed as -DVERSION=\"0.3\"
#ifdef VERSION
	#define TV_VERSION VERSION
#else
	#define TV_VERSION "0.5.0"
#endif

/*
 * Purpose:  Document/view architecture demo for wxWindows class library - MDI
 */
 

#include "tdoc.h"
#include "tview.h"


// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
	#pragma hdrstop
#endif

#ifndef WX_PRECOMP
	#include "wx/wx.h"
#endif

#if !wxUSE_DOC_VIEW_ARCHITECTURE
	#error You must set wxUSE_DOC_VIEW_ARCHITECTURE to 1 in setup.h!
#endif

#if !wxUSE_MDI_ARCHITECTURE
	#error You must set wxUSE_MDI_ARCHITECTURE to 1 in setup.h!
#endif



#include <wx/toolbar.h>

#if defined(__WXGTK__) || defined(__WXMOTIF__) || defined (__WXX11__)
// For Unix use 24x24 xpm files modelled on GNOME
//    #include "mondrian.xpm"
   #include "bitmaps/new.xpm"
	#include "bitmaps/open.xpm"
//    #include "bitmaps/save.xpm"
	#include "bitmaps/saveas.xpm"
	#include "bitmaps/copy.xpm"
//    #include "bitmaps/cut.xpm"
//    #include "bitmaps/paste.xpm"
	#include "bitmaps/print_preview.xpm"
	#include "bitmaps/print.xpm"
//    #include "bitmaps/help.xpm"
	#include "bitmaps/previous.xpm"
	#include "bitmaps/next.xpm"
	#include "bitmaps/slanttree.xpm"
	#include "bitmaps/rectangletree.xpm"
	#include "bitmaps/phylogram.xpm"
	#include "bitmaps/paste.xpm"
	#include "bitmaps/zoomin.xpm"
	#include "bitmaps/zoomout.xpm"
	#include "bitmaps/internal.xpm"

#elif defined(__WXMAC__)
// For MacOS we use 16x15 xpm files modelled on Windows
	#include "bitmaps/mac/new.xpm"
	#include "bitmaps/mac/open.xpm"
//    #include "bitmaps/save.xpm"
	#include "bitmaps/mac/saveas.xpm"
	#include "bitmaps/mac/copy.xpm"
//    #include "bitmaps/cut.xpm"
//    #include "bitmaps/paste.xpm"
	#include "bitmaps/mac/print_preview.xpm"
	#include "bitmaps/mac/print.xpm"
//    #include "bitmaps/help.xpm"
	#include "bitmaps/mac/previous.xpm"
	#include "bitmaps/mac/next.xpm"
	#include "bitmaps/mac/slanttree.xpm"
	#include "bitmaps/mac/rectangletree.xpm"
	#include "bitmaps/mac/phylogram.xpm"
	#include "bitmaps/mac/paste.xpm"
	#include "bitmaps/mac/zoomin.xpm"
	#include "bitmaps/mac/zoomout.xpm"
	#include "bitmaps/internal.xpm"


#endif

#include "tv.h"


#if wxUSE_DRAG_AND_DROP

// Handle dropping files
#include "wx/dnd.h"

class DnDFile : public wxFileDropTarget
{
public:
    DnDFile(wxDocManager *pOwner) { m_pOwner = pOwner; }

    virtual bool OnDropFiles(wxCoord x, wxCoord y,
                             const wxArrayString& filenames);

private:
    wxDocManager *m_pOwner;
};

bool DnDFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
    size_t nFiles = filenames.GetCount();
 
    for ( size_t n = 0; n < nFiles; n++ ) 
        m_pOwner->CreateDocument (filenames[n].c_str(), wxDOC_SILENT);

    return TRUE;
}
#endif // wxUSE_DRAG_AND_DROP

MyFrame *frame = (MyFrame *) NULL;

IMPLEMENT_APP(MyApp)

MyApp::MyApp(void)
{
	m_docManager = (wxDocManager *) NULL;
}

#ifdef __WXMAC__
void MyApp::MacOpenFile(const wxString &fileName)
{
	m_docManager->CreateDocument (fileName, wxDOC_SILENT);
}
#endif

bool MyApp::OnInit(void)
{
	// Optional command line argument is name of file to open
#if wxUSE_UNICODE
	wchar_t *InputFile = NULL;
#else
	char *InputFile = NULL;
#endif
	// Read input/output files
	if (argc > 1)
	{
		InputFile = argv[1];
	}

	//// Create a document manager                  
	m_docManager = new wxDocManager;

	//// Create a template relating drawing documents to their views
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("NEXUS tree file"), wxT("*.tre"), wxT(""), wxT("tre"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("NEXUS tree file"), wxT("*.trees"), wxT(""), wxT("trees"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("NEXUS file"), wxT("*.nex"), wxT(""), wxT("nex"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("TreeBASE"), wxT("*.acgi"), wxT(""), wxT("acgi"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("MrBayes consensus"), wxT("*.con"), wxT(""), wxT("con"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("MrBayes tree file"), wxT("*.t"), wxT(""), wxT("t"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("CLUSTALX tree"), wxT("*.ph"), wxT(""), wxT("ph"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("CLUSTALX bootstrap tree"), wxT("*.phb"), wxT(""), wxT("phb"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("CLUSTALX guide tree"), wxT("*.dnd"), wxT(""), wxT("dnd"), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("PHYLIP tree"), wxT("*.*"), wxT(""), wxT(""), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
	(void) new wxDocTemplate((wxDocManager *) m_docManager, wxT("Any document"), wxT("*"), wxT(""), wxT(""), wxT("Tree Doc"), wxT("Tree View"),
		CLASSINFO(TDocument), CLASSINFO(TView));
          
	//// Create the main frame window
#ifdef __WXMAC__
	// Create the frame window off screen (5000, 5000) so we don't see it.
	frame = new MyFrame((wxDocManager *) m_docManager, (wxFrame *) NULL, (const wxString) wxT("Toolbar"), 
	wxPoint(5000,5000), wxSize(400, 400), 
	wxDEFAULT_FRAME_STYLE);
#else
	frame = new MyFrame((wxDocManager *) m_docManager, (wxFrame *) NULL, (const wxString) wxT("TreeViewX"), 
	wxPoint(0, 0), wxSize(500, 400), wxDEFAULT_FRAME_STYLE);
#endif

	//// Give it an icon (this is ignored in MDI mode: uses resources)
#ifdef __WXMSW__
	frame->SetIcon(wxIcon("app")); // use the name of the icon in the resource file
#endif
#if defined(__WXGTK__) || defined(__WXMOTIF__)
	frame->SetIcon(wxIcon(wxT("bitmaps/treeview.xbm")));
#endif

#ifdef __WXMAC__
    wxApp::s_macAboutMenuItemId = DOCVIEW_ABOUT;
#endif

	//// Make a menubar
	wxMenu *file_menu = new wxMenu;
	wxMenu *edit_menu = (wxMenu *) NULL;

	file_menu->Append(wxID_NEW, wxT("&New...\tCtrl-N"));
	file_menu->Append(wxID_OPEN, wxT("&Open...\tCtrl-O"));

#ifndef __WXMAC__
	file_menu->AppendSeparator();
	file_menu->Append(wxID_EXIT, wxT("E&xit\tAlt-X"));
#endif

	// A nice touch: a history of files visited. Use this menu.
	m_docManager->FileHistoryUseMenu(file_menu);

	wxMenu *help_menu = new wxMenu;
	help_menu->Append(DOCVIEW_ABOUT, wxT("&About\tF1"));

	wxMenuBar *menu_bar = new wxMenuBar;

	menu_bar->Append(file_menu, wxT("&File"));
	if (edit_menu)
	menu_bar->Append(edit_menu, wxT("&Edit"));
	menu_bar->Append(help_menu, wxT("&Help"));

	//// Associate the menu bar with the frame
	frame->SetMenuBar(menu_bar);

	// Application window status bar
	frame->CreateStatusBar();

#if wxUSE_DRAG_AND_DROP
	// Frame window is a drop target
    frame->SetDropTarget(new DnDFile(m_docManager));
#endif

#ifndef __WXMAC__
	// If we don't centre the frame we don't see the it under MacOS X (which is good because
	// the MDI model on the Mac should not show a frame).
	frame->Centre(wxBOTH); 
#endif
	frame->Show(TRUE);

	SetTopWindow(frame);
	
	// Open file passed on command line
	if (InputFile)
		m_docManager->CreateDocument (InputFile, wxDOC_SILENT);
	
	return TRUE;
}

int MyApp::OnExit(void)
{
	delete m_docManager;
	return 0;
}

/*
 * Centralised code for creating a document frame.
 * Called from view.cpp, when a view is created.
 */
 
wxMDIChildFrame *MyApp::CreateChildFrame(wxDocument *doc, wxView *view, bool isCanvas)
{
	//// Make a child frame
	wxDocMDIChildFrame *subframe = new wxDocMDIChildFrame(doc, view, GetMainFrame(), -1, wxT("Child Frame"),
	wxPoint(10, 10), wxSize(300, 300), wxDEFAULT_FRAME_STYLE);

#ifdef __WXMSW__
	subframe->SetIcon(wxString("doc"));
#endif
#ifdef __X__
	subframe->SetIcon(wxIcon("bitmaps/doc.xbm"));
#endif

	//// Make a menubar
	wxMenu *file_menu = new wxMenu;

	file_menu->Append(wxID_NEW, wxT("&New...\tCtrl-N"));
	file_menu->Append(wxID_OPEN, wxT("&Open...\tCtrl-O"));
	file_menu->Append(wxID_CLOSE, wxT("&Close\tCtrl-W"));
	file_menu->Append(wxID_SAVE, wxT("&Save\tCtrl-S"));
	file_menu->Append(wxID_SAVEAS, wxT("Save &As...\tShift-Ctrl-S"));

	if (isCanvas)
	{
#ifdef __WXMSW__
		file_menu->AppendSeparator();
		file_menu->Append(SAVEAS_PICTURE_CMD, wxT("Save As Picture..."), wxT("Save picture of tree to metafile"));
#else
	#if USE_SVG
		file_menu->AppendSeparator();
		file_menu->Append(SAVEAS_PICTURE_CMD, wxT("Save As Picture..."), wxT("Save picture of tree to SVG file"));
	#endif
#endif
		file_menu->AppendSeparator();
		file_menu->Append(wxID_PRINT, wxT("&Print...\tCtrl-P"));
		file_menu->Append(wxID_PRINT_SETUP, wxT("Page &Setup...\tShift-Ctrl-P"));
		file_menu->Append(wxID_PREVIEW, wxT("Print Pre&view"));
	}

#ifndef __WXMAC__
	file_menu->AppendSeparator();
	file_menu->Append(wxID_EXIT, wxT("E&xit"));
#endif

	wxMenu *edit_menu = (wxMenu *) NULL;

    edit_menu = new wxMenu;
    edit_menu->Append(wxID_UNDO, wxT("&Undo\tCtrl-Z"));
    edit_menu->Append(wxID_REDO, wxT("&Redo\tShift-Ctrl-Z"));
    edit_menu->AppendSeparator();
    edit_menu->Append(wxID_CUT, wxT("Cu&t\tCtrl+X"));
    edit_menu->Append(wxID_COPY, wxT("&Copy\tCtrl+C"), wxT("Copy tree picture to clipboard"));
    edit_menu->Append(COPY_AS_TEXT_CMD, wxT("&Copy As Text"), wxT("Copy tree description to clipboard"));
    edit_menu->Append(wxID_PASTE, wxT("&Paste\tCtrl+V"), wxT("Paste tree from clipboard"));
    edit_menu->Append(wxID_CLEAR, wxT("&Delete"));
	
	edit_menu->Enable (wxID_UNDO, FALSE);
	edit_menu->Enable (wxID_REDO, FALSE);
	edit_menu->Enable (wxID_CUT, FALSE);
	edit_menu->Enable (wxID_COPY, FALSE);
	edit_menu->Enable (COPY_AS_TEXT_CMD, FALSE);
	edit_menu->Enable (wxID_PASTE, FALSE);
	edit_menu->Enable (wxID_CLEAR, FALSE);

  wxMenu *tree_menu = (wxMenu *) NULL;

	if (isCanvas)
	{
		tree_menu = new wxMenu;
		// These next three items are checkable
		tree_menu->Append(SLANTED_TREE_CMD, wxT("&Slanted Cladogram"), wxT("Draw tree as a slanted cladogram"), TRUE);
		tree_menu->Append(RECTANGULAR_TREE_CMD, wxT("&Rectangular Cladogram"), wxT("Draw tree as a rectangular cladogram"), TRUE);
		tree_menu->Append(PHYLOGRAM_CMD, wxT("&Phylogram"), wxT("Draw tree as a phylogram"), TRUE);
		tree_menu->AppendSeparator();
		tree_menu->Append(INTERNAL_LABELS_CMD, wxT("&Show internal node labels"), wxT("Show internal node labels"), TRUE);
		tree_menu->Append(LEAF_FONT_CMD, wxT("&Leaf Font..."), wxT("Set font for displaying leaves"));
		tree_menu->AppendSeparator();
		tree_menu->Append(CHOOSE_TREE_CMD, wxT("&Choose tree..."), wxT("Choose a tree"));
		tree_menu->AppendSeparator();
		tree_menu->Append(ORDER_TREE_CMD, wxT("&Order tree..."), wxT("Order the tree"));
		tree_menu->AppendSeparator();
		tree_menu->Append(ZOOM_IN_CMD, wxT("Zoom in\tCtrl+["), wxT("Zoom in"));
		tree_menu->Append(ZOOM_OUT_CMD, wxT("Zoom out\tCtrl+]"), wxT("Zoom out"));
		tree_menu->Append(ZOOM_TO_FIT_CMD, wxT("Zoom to fit"), wxT("View whole tree"));
	}

	wxMenu *help_menu = new wxMenu;

	help_menu->Append(DOCVIEW_ABOUT, wxT("&About"));

	wxMenuBar *menu_bar = new wxMenuBar;

	menu_bar->Append(file_menu, wxT("&File"));
	menu_bar->Append(edit_menu, wxT("&Edit"));
	if (isCanvas)
		menu_bar->Append(tree_menu, wxT("&Trees"));
	menu_bar->Append(help_menu, wxT("&Help"));

	//// Associate the menu bar with the frame
	subframe->SetMenuBar(menu_bar);

#if defined(__WXMSW__) || defined(__WXMAC__)
	// Status bar for the view window
	subframe->CreateStatusBar();
#endif

#ifdef __WXMAC__
	// Create a toolbar in the view window
	subframe->CreateToolBar(wxNO_BORDER | wxTB_FLAT | wxTB_HORIZONTAL);
	frame->InitToolBar(subframe->GetToolBar());	
#endif


	return subframe;
}

/*
 * This is the top-level window of the application.
 */
 
IMPLEMENT_CLASS(MyFrame, wxDocMDIParentFrame)
BEGIN_EVENT_TABLE(MyFrame, wxDocMDIParentFrame)
    EVT_MENU(DOCVIEW_ABOUT, MyFrame::OnAbout)
END_EVENT_TABLE()

MyFrame::MyFrame(wxDocManager *manager, wxFrame *frame, const wxString& title,
    const wxPoint& pos, const wxSize& size, long type):
  wxDocMDIParentFrame(manager, frame, -1, title, pos, size, type, wxT("myFrame"))
{
//  editMenu = (wxMenu *) NULL;

	// Toolbar
	CreateToolBar(wxNO_BORDER | wxTB_FLAT | wxTB_HORIZONTAL);
	InitToolBar(GetToolBar());


	// Accelerators
	wxAcceleratorEntry entries[3];
	entries[0].Set(wxACCEL_CTRL, (int) 'N', wxID_NEW);
	entries[1].Set(wxACCEL_CTRL, (int) 'O', wxID_OPEN);
	entries[2].Set(wxACCEL_CTRL, (int) 'X', wxID_EXIT);
	wxAcceleratorTable accel(3, entries);
	SetAcceleratorTable(accel);

}

#ifdef __WXMSW__
	#define NUM_BITMAPS	15
#else
	#define NUM_BITMAPS	15
#endif

void MyFrame::InitToolBar(wxToolBar* toolBar)
{
	wxBitmap* bitmaps[NUM_BITMAPS];

	#ifdef __WXMSW__
		bitmaps[0] = new wxBitmap("open", wxBITMAP_TYPE_RESOURCE);
		bitmaps[1] = new wxBitmap("saveas", wxBITMAP_TYPE_RESOURCE);
		bitmaps[2] = new wxBitmap("copy", wxBITMAP_TYPE_RESOURCE);
		bitmaps[3] = new wxBitmap("print_preview", wxBITMAP_TYPE_RESOURCE);
		bitmaps[4] = new wxBitmap("print", wxBITMAP_TYPE_RESOURCE);
		bitmaps[5] = new wxBitmap("previous", wxBITMAP_TYPE_RESOURCE);
		bitmaps[6] = new wxBitmap("next", wxBITMAP_TYPE_RESOURCE);
		bitmaps[7] = new wxBitmap("slanttree", wxBITMAP_TYPE_RESOURCE);
		bitmaps[8] = new wxBitmap("rectangletree", wxBITMAP_TYPE_RESOURCE);
		bitmaps[9] = new wxBitmap("phylogram", wxBITMAP_TYPE_RESOURCE);
		bitmaps[10] = new wxBitmap("paste", wxBITMAP_TYPE_RESOURCE);
 		bitmaps[11] = new wxBitmap("new", wxBITMAP_TYPE_RESOURCE);
		bitmaps[12] = new wxBitmap("zoomin", wxBITMAP_TYPE_RESOURCE);
		bitmaps[13] = new wxBitmap("zoomout", wxBITMAP_TYPE_RESOURCE);
		bitmaps[14] = new wxBitmap("internal", wxBITMAP_TYPE_RESOURCE);
    #else
		bitmaps[0] = new wxBitmap( open_xpm );
		bitmaps[1] = new wxBitmap( saveas_xpm );
		bitmaps[2] = new wxBitmap( copy_xpm );
		bitmaps[3] = new wxBitmap( print_preview_xpm );
		bitmaps[4] = new wxBitmap( print_xpm );
		bitmaps[5] = new wxBitmap( previous_xpm );
		bitmaps[6] = new wxBitmap( next_xpm );
		bitmaps[7] = new wxBitmap( slanttree_xpm );
		bitmaps[8] = new wxBitmap( rectangletree_xpm );
		bitmaps[9] = new wxBitmap( phylogram_xpm );
		bitmaps[10] = new wxBitmap( paste_xpm);
		bitmaps[11] = new wxBitmap( new_xpm);
		bitmaps[12] = new wxBitmap( zoomin_xpm);
		bitmaps[13] = new wxBitmap( zoomout_xpm);
		bitmaps[14] = new wxBitmap( internal_xpm);
  	#endif

#ifdef __WXMSW__
	int width = 24;
#else
	int width = 16;
#endif
	int currentX = 5;

#ifdef __WXMAC__
	// The toolbar for the Mac is associated with each document window, not the MDI frame,
	// so we don't display Open and Save as buttons
#else
	toolBar->AddTool( wxID_NEW, *(bitmaps[11]), wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("New document"));
	currentX += width + 5;
	toolBar->AddTool( wxID_OPEN, *(bitmaps[0]), wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Open document"));
	currentX += width + 5;
	toolBar->AddTool(wxID_SAVEAS, *bitmaps[1], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Save document"));
	currentX += width + 5;
#endif
	toolBar->AddTool(wxID_COPY, *bitmaps[2], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Copy tree to clipboard"));
	currentX += width + 5;
	toolBar->AddTool(wxID_PASTE, *bitmaps[10], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Paste tree descriiption from clipboard"));
	currentX += width + 5;
	// --
	toolBar->AddSeparator();
	toolBar->AddTool(wxID_PREVIEW, *bitmaps[3], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Print preview"));
	toolBar->AddTool(wxID_PRINT, *bitmaps[4], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Print"));
	//--
	toolBar->AddSeparator();
	toolBar->AddTool(PREVIOUS_TREE_CMD, *bitmaps[5], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Previous tree"));
	toolBar->AddTool(NEXT_TREE_CMD, *bitmaps[6], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Next tree"));
	//--
	toolBar->AddSeparator();
	toolBar->AddTool(SLANTED_TREE_CMD, *bitmaps[7], wxNullBitmap, TRUE, currentX, -1, (wxObject *) NULL, wxT("Slanted cladogram"));
	toolBar->AddTool(RECTANGULAR_TREE_CMD, *bitmaps[8], wxNullBitmap, TRUE, currentX, -1, (wxObject *) NULL, wxT("Rectangular cladogram"));
	toolBar->AddTool(PHYLOGRAM_CMD, *bitmaps[9], wxNullBitmap, TRUE, currentX, -1, (wxObject *) NULL, wxT("Phylogram"));
	//--
 	toolBar->AddSeparator();
	toolBar->AddTool(INTERNAL_LABELS_CMD, *bitmaps[14], wxNullBitmap, TRUE, currentX, -1, (wxObject *) NULL, wxT("Show internal node labels"));
	//--
	toolBar->AddSeparator();
	toolBar->AddTool(ZOOM_IN_CMD, *bitmaps[12], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Zoom in"));
	toolBar->AddTool(ZOOM_OUT_CMD, *bitmaps[13], wxNullBitmap, FALSE, currentX, -1, (wxObject *) NULL, wxT("Zoom out"));

	toolBar->Realize();

	toolBar->EnableTool (wxID_COPY, false);
	toolBar->EnableTool (wxID_PASTE, false);
	toolBar->EnableTool (PREVIOUS_TREE_CMD, false);
	toolBar->EnableTool (NEXT_TREE_CMD, false);
	toolBar->EnableTool (PREVIOUS_TREE_CMD, false);
	toolBar->EnableTool (SLANTED_TREE_CMD, false);
	toolBar->EnableTool (RECTANGULAR_TREE_CMD, false);
	toolBar->EnableTool (PHYLOGRAM_CMD, false);
	toolBar->EnableTool (ZOOM_IN_CMD, false);
	toolBar->EnableTool (ZOOM_OUT_CMD, false);
	toolBar->EnableTool (INTERNAL_LABELS_CMD, false);

	for (int i = 0; i < NUM_BITMAPS; i++)
		delete bitmaps[i];
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event) )
{
	wxString msg;
	msg = wxT("TreeViewX ");
	msg += wxT("Version ");
	msg += wxT(TV_VERSION);
	msg += wxT("\n\nby Roderic D. M. Page (r.page@bio.gla.ac.uk)");
	msg += wxT("\n\nBuilt using ");
	msg +=  wxVERSION_STRING;
	msg +=  wxT(" for ");
#ifdef __WXMAC__
	msg += wxT("Macintosh");
#else
	msg +=  wxGetOsDescription();
#endif
	(void)wxMessageBox(msg, wxT("About TreeViewX"));
}

// Creates a canvas. Called from view.cpp when a new drawing
// view is created.
MyCanvas *MyFrame::CreateCanvas(wxView *view, wxMDIChildFrame *parent)
{
	int width, height;
	parent->GetClientSize(&width, &height);

	// Non-retained canvas
	MyCanvas *canvas = new MyCanvas(view, parent, wxPoint(0, 0), wxSize(width, height), 0);
	//  canvas->SetCursor(wxCursor(wxCURSOR_PENCIL));

	// Give it scrollbars
	canvas->SetScrollbars(0, 0, 0, 0);

	return canvas;
}

MyFrame *GetMainFrame(void)
{
	return frame;
}

