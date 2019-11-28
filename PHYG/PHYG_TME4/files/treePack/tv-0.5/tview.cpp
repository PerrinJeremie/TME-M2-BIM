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
 
 // $Id: tview.cpp,v 1.40 2007/03/09 13:59:19 rdmp1c Exp $

#ifdef __GNUG__
// #pragma implementation
#endif

#include <strstream>

// TreeLib (must go before any wx stuff to avoid clash between STL and Windows min/max)
#include "TreeLib.h"
#include "Parse.h"
#include "treedrawer.h"
#include "treeorder.h"

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

#include <wx/clipbrd.h>	
#include <wx/fontdlg.h>
#include <wx/metafile.h>
#include <wx/filename.h>

// Use the native Mac OS X font dialog
#ifdef __WXMAC__
	#undef wxFontDialog
	#include "wx/mac/fontdlg.h"
#endif

// Pictures

#ifdef __WXMSW__
#else
#ifdef USE_SVG
	#include <wx/svg/dcsvg.h>
#endif
#endif

// GUI interface
#include "tv.h"
#include "tdoc.h"
#include "tview.h"

// Dialogs
#include "treeorder_dialog.h"

#if wxUSE_PRINTING_ARCHITECTURE
class WXDLLEXPORT wxTVPrintout : public wxDocPrintout
{
public:
	wxTVPrintout (wxView *view = (wxView *) NULL, const wxString& title = wxT("Printout")) : wxDocPrintout (view, title) {};
	bool OnPrintPage(int page);
};
#endif // wxUSE_PRINTING_ARCHITECTURE

//------------------------------------------------------------------------------
// The wxWindows Printing Framework scales the printout to the screen image,
// whereas I want to print the tree on the whole page, regardless of how it looks
// on screen. To do this we extend wxDocPrintout to compute the size of the page
// and to set the flag "printing" in the associated TView window before calling
// the ancestral OnPrintPage. This lets TView:OnDraw know whether to scale the tree 
// for the screen or for the printed page.
bool wxTVPrintout::OnPrintPage(int page)
{
	bool result = false;

	TView *t = (TView *)GetView();
	if (t)
	{
		int w_mm, h_mm;
		GetPageSizeMM (&w_mm, &h_mm);

		// Page size in points
		t->page_width = (int) ((float)(w_mm * 72.0) / 25.4);
		t->page_height = (int) ((float)(h_mm * 72.0) / 25.4);
		t->printing = true;
	}

	result = wxDocPrintout::OnPrintPage (page);

	if (t)
	{
		t->printing = false;
	}

	return result;
}




IMPLEMENT_DYNAMIC_CLASS(TView, wxView)

/*
BEGIN_EVENT_TABLE(DrawingView, wxView)
    EVT_MENU(DOODLE_CUT, DrawingView::OnCut)
END_EVENT_TABLE()   */


	
	


//------------------------------------------------------------------------------
TView::TView ()
{
	canvas = (MyCanvas *) NULL; 
	frame = (wxMDIChildFrame *) NULL; 
	TreeStyleCmd = SLANTED_TREE_CMD;
	Order = DEFAULT_ORDER;
	
	LeafFont = *wxSWISS_FONT;

	printing = false;

	magnification = 1;
	mShowInternalLabels = false;
}

//------------------------------------------------------------------------------
wxPrintout *TView::OnCreatePrintout()
{
	wxDocPrintout *p = new wxTVPrintout (this);
	return p;
}

//------------------------------------------------------------------------------
// What to do when a view is created. Creates actual
// windows for displaying the view.
bool TView::OnCreate(wxDocument *doc, long WXUNUSED(flags) )
{
    frame = wxGetApp().CreateChildFrame(doc, this, TRUE);
    frame->SetTitle(wxT("TView"));

    canvas = GetMainFrame()->CreateCanvas(this, frame);
#ifdef __X__
    // X seems to require a forced resize
    int x, y;
    frame->GetSize(&x, &y);
    frame->SetSize(-1, -1, x, y);
#endif
    frame->Show(TRUE);
    Activate(TRUE);

    return TRUE;
}


//------------------------------------------------------------------------------
void TView::OnSavePicture (wxCommandEvent& WXUNUSED(event))
{
#ifdef __WXMSW__
	// Use metafiles for Windows
	wxString pictureFileName = GetFrame()->GetTitle();
	pictureFileName += wxT(".emf");
	wxFrame *f = GetMainFrame();
   	wxFileDialog dialog((wxWindow *)f, wxT("Save Picture as"), wxT(""), pictureFileName,
        wxT("Enhanced metafile (*.emf)|*.emf"),
        wxSAVE|wxOVERWRITE_PROMPT);
    
    if (dialog.ShowModal() == wxID_OK)
    {
		wxMetafileDC pictureDC (dialog.GetPath(), 600, 650) ;
		OnDraw (&pictureDC);
		pictureDC.Close();
	}
#else
#ifdef USE_SVG
	wxString pictureFileName = frame->GetTitle();
	// Use SVG on other platforms	
	pictureFileName += wxT(".svg");
	#ifdef __WXMAC__
	wxFrame *f = GetMainFrame();
	#else
	wxFrame *f = GetMainFrame();
	#endif 
    wxFileDialog dialog((wxWindow *)f, wxT("Save Picture as"), wxT(""), pictureFileName,
        wxT("SVG vector picture files (*.svg)|*.svg"),
        wxSAVE|wxOVERWRITE_PROMPT);

    if (dialog.ShowModal() == wxID_OK)
    {
		wxSVGFileDC pictureDC (dialog.GetPath(), 600, 650) ;
		OnDraw (&pictureDC);
	}
#endif
#endif
}


//------------------------------------------------------------------------------
// Sneakily gets used for default print/preview
// as well as drawing on the screen.
void TView::OnDraw(wxDC *dc)
{
	if (p.GetNumTrees() > 0)
	{
		dc->SetFont(LeafFont);
		wxPen pen (*wxBLACK_PEN);
		pen.SetCap (wxCAP_PROJECTING);
		dc->SetPen(pen);

		Tree t = p.GetCurrentTree();
		t.Update();
		
		int width, height;
		if (printing)
		{
			width = page_width;
			height = page_height;
#if defined(__WXGTK__) || defined(__WXMOTIF__)
			height -= 72;
#endif 
		}
		else
			if (magnification == 1)
			{
				canvas->GetClientSize (&width, &height);
			}
			else
			{
				int ignore;
				canvas->GetVirtualSize (&ignore, &height);
				canvas->GetClientSize (&width, &ignore);
			}		wxRect r (0, 0, width, height);

#if defined(__WXGTK__) || defined(__WXMOTIF__)
		if (printing)
		{
			r.Offset (0, 72);
		}
#endif 

			
		// Ensure a little space around the tree
		r.Inflate (-20, -20);

   		t.MakeNodeList();
    	int maxlen = 0;
   		for (int i = 0; i < t.GetNumLeaves(); i++)
		{
			wxCoord w, h;
			Node *p;
			p = t[i];
   			wxString s;
#if wxUSE_UNICODE
			char buf[256];
			strcpy (buf, p->GetLabel().c_str());
			wchar_t wbuf[256];
			mbstowcs (wbuf, buf, size_t(wbuf));
			s << wbuf;
#else
   			s << p->GetLabel().c_str();
#endif
   			dc->GetTextExtent (s, &w, &h);
    		if (w > maxlen)
				maxlen =  w;
		}
		r.SetWidth (r.GetWidth() - maxlen);
		
		if (Order != DEFAULT_ORDER)
		{
			TreeOrder *order = NULL;
			switch (Order)
			{
				case ALPHA_ORDER:
					order = new AlphaOrder (&t);
					break;
				case LEFT_ORDER:
					order = new LeftOrder (&t);
					break;
				case RIGHT_ORDER:
					order = new RightOrder (&t);
					break;
			}
			order->Order ();
			delete order;
		}

		
		TreeDrawer *td = NULL;
		switch (TreeStyleCmd)
		{
			case SLANTED_TREE_CMD:
				td = new TreeDrawer (&t);
				break;
			case RECTANGULAR_TREE_CMD:
				td = new RectangleTreeDrawer (&t);
				break;
			case PHYLOGRAM_CMD:
				td = new PhylogramDrawer (&t);
				break;
		}
		td->SetDC (dc);
        td->SetRect (r.GetLeft(), r.GetTop(), r.GetRight(), r.GetBottom());
        td->CalcCoordinates ();
        td->SetPenWidth (2);
        td->SetDrawInternalLabels (mShowInternalLabels);
        td->Draw ();
        delete td;
	}
}

//------------------------------------------------------------------------------
void TView::OnUpdate(wxView *WXUNUSED(sender), wxObject *WXUNUSED(hint))
{
	if (canvas)
		canvas->Refresh();

/* Is the following necessary?
#ifdef __WXMSW__
  if (canvas)
    canvas->Refresh();
#else
  if (canvas)
    {
      wxClientDC dc(canvas);
      dc.Clear();
      OnDraw(& dc);
    }
#endif
*/
}

//------------------------------------------------------------------------------
// Clean up windows used for displaying the view.
bool TView::OnClose(bool deleteWindow)
{
	if (!GetDocument()->Close())
		return FALSE;

	// Clear the canvas in  case we're in single-window mode,
	// and the canvas stays.
#if (wxMINOR_VERSION > 4) // from 2.5 Clear is deprecated
	canvas->ClearBackground();
#else
	canvas->Clear();
#endif	
	
	canvas->view = (wxView *) NULL;
	canvas = (MyCanvas *) NULL;

	wxString s(wxTheApp->GetAppName());
	if (frame)
		frame->SetTitle(s);

	SetFrame((wxFrame*)NULL);

	Activate(FALSE);

	if (deleteWindow)
	{
		delete frame;
		return TRUE;
	}
	return TRUE;
}

//------------------------------------------------------------------------------
bool TView::ReadTrees (const wxString &filename)
{
	bool result = false;
	ifstream f;
	f.open (filename.mb_str(), ios::binary | ios::in);
	if (p.ReadTrees (f))
	{
		// Set window name to that of file
		wxFileName fn (filename);
		wxString title = fn.GetName();
		if (fn.HasExt())
		{
			title += wxT(".");
			title += fn.GetExt();
		}
		frame->SetTitle(title);

		result = true;
	}
	else
	{
		wxMessageBox(wxT("Failed to read trees"), wxT("ReadTrees"),  wxICON_EXCLAMATION | wxOK);
		result = false;
	}
	f.close();
	return result;
}

//------------------------------------------------------------------------------
bool TView::WriteTrees (const wxString &filename)
{
	bool result = false;
	ofstream f (filename.mb_str());
	if (p.WriteTrees (f))
	{
		// Set window name to that of file
		wxFileName fn (filename);
		wxString title = fn.GetName();
		if (fn.HasExt())
		{
			title += wxT(".");
			title += fn.GetExt();
		}
		frame->SetTitle(title);

		result = true;
	}
	else
	{
		wxMessageBox(wxT("Failed to write trees"), wxT("TreeView X"),  wxICON_EXCLAMATION | wxOK);
		result = false;
	}
	f.close();
	return result;
}

//------------------------------------------------------------------------------
void TView::OnActivateView (bool activate, wxView *activeView, wxView *deactiveView)
{
	if (activate)
	{
#ifdef __WXMAC__
		wxToolBar* toolBar = frame->GetToolBar();
#else
		wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
		wxMenuBar* menuBar = frame->GetMenuBar();
				
		if (p.GetNumTrees() > 0)
		{
			toolBar->EnableTool (PREVIOUS_TREE_CMD, (p.GetCurrentTreeNumber() > 0));
			toolBar->EnableTool (NEXT_TREE_CMD, (p.GetCurrentTreeNumber() < (p.GetNumTrees() - 1)));	
			toolBar->EnableTool (SLANTED_TREE_CMD, true);		
			toolBar->EnableTool (RECTANGULAR_TREE_CMD, true);
			toolBar->EnableTool (PHYLOGRAM_CMD, p.CurrentTreeHasEdgeLengths());
			toolBar->ToggleTool (SLANTED_TREE_CMD, FALSE );		
			toolBar->ToggleTool (RECTANGULAR_TREE_CMD, FALSE );		
			toolBar->ToggleTool (PHYLOGRAM_CMD, FALSE );

#ifndef __WXMOTIF__
			menuBar->Enable (COPY_AS_TEXT_CMD, true);
#endif	
			
			menuBar->Enable (SLANTED_TREE_CMD, true);		
			menuBar->Enable (RECTANGULAR_TREE_CMD, true);
			menuBar->Enable (PHYLOGRAM_CMD, p.CurrentTreeHasEdgeLengths());
			menuBar->Check (SLANTED_TREE_CMD, false );		
			menuBar->Check (RECTANGULAR_TREE_CMD, false );		
			menuBar->Check (PHYLOGRAM_CMD, false );				
			menuBar->Enable (CHOOSE_TREE_CMD, (p.GetNumTrees() > 1));

 			menuBar->Enable (INTERNAL_LABELS_CMD, p.CurrentTreeHasInternallabels());
			menuBar->Check (INTERNAL_LABELS_CMD, mShowInternalLabels );
			toolBar->EnableTool (INTERNAL_LABELS_CMD, p.CurrentTreeHasInternallabels());
			toolBar->ToggleTool (INTERNAL_LABELS_CMD, mShowInternalLabels );


			ShowTreeName ();
			TreeStyleMenu (TreeStyleCmd);	
		}
		else
		{
			toolBar->EnableTool (PREVIOUS_TREE_CMD, false);
			toolBar->EnableTool (NEXT_TREE_CMD, false);
			toolBar->EnableTool (PREVIOUS_TREE_CMD, false);
			toolBar->EnableTool (SLANTED_TREE_CMD, false);		
			toolBar->EnableTool (RECTANGULAR_TREE_CMD, false);
			toolBar->EnableTool (PHYLOGRAM_CMD, false);		
			
			menuBar->Enable (COPY_AS_TEXT_CMD, false);
			menuBar->Enable (SLANTED_TREE_CMD, false);		
			menuBar->Enable (RECTANGULAR_TREE_CMD, false);
			menuBar->Enable (PHYLOGRAM_CMD, false);
			menuBar->Enable (INTERNAL_LABELS_CMD, false );

		}		
	}
}

//------------------------------------------------------------------------------
void TView::PreviousTree (wxCommandEvent& WXUNUSED(event))
{
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	int n = p.GetCurrentTreeNumber ();
	if (n >= 0)
	{
		n--;
		p.SetCurrentTreeNumber (n);
		if (n == 0)
			toolBar->EnableTool (PREVIOUS_TREE_CMD, false);
		toolBar->EnableTool (NEXT_TREE_CMD, true);
		TreeStyleMenu(TreeStyleCmd);
		ShowTreeName ();
		OnUpdate (this, NULL); // Ensure we redraw window
	}
}

//------------------------------------------------------------------------------
void TView::NextTree (wxCommandEvent& WXUNUSED(event))
{
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	int n = p.GetCurrentTreeNumber ();
	
	if (n < p.GetNumTrees () - 1)
	{
		n++;
		p.SetCurrentTreeNumber (n);
		if (n == (p.GetNumTrees () - 1))
			toolBar->EnableTool (NEXT_TREE_CMD, false);
		toolBar->EnableTool (PREVIOUS_TREE_CMD, true);
		TreeStyleMenu(TreeStyleCmd);		
		ShowTreeName ();
		OnUpdate (this, NULL); // Ensure we redraw window
	}
}

//------------------------------------------------------------------------------
void TView::ShowTreeName ()
{
	wxString txt = wxT("");
	Tree t = p.GetCurrentTree ();
	if (t.GetName() != "")
	{
#if wxUSE_UNICODE
		char buf[256];
		strcpy (buf, t.GetName().c_str());
		wchar_t wbuf[256];
		mbstowcs (wbuf, buf, size_t(wbuf));
		txt << wbuf;
#else
		txt << t.GetName().c_str();
#endif
		
		txt << wxT(" ");
	}
		
	txt << wxT("(") << (p.GetCurrentTreeNumber() + 1) << wxT("/") << p.GetNumTrees() << wxT(")");
#if defined(__WXMSW__) || defined(__WXMAC__)
	// Status bar for the view window
	frame->SetStatusText (txt, 0);
#else
	GetMainFrame()->SetStatusText (txt, 0);
#endif
}

//------------------------------------------------------------------------------
void TView::TreeStyleMenu (int style)
{
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	wxMenuBar* menuBar = frame->GetMenuBar();

	toolBar->ToggleTool( TreeStyleCmd, FALSE );
	menuBar->Check (TreeStyleCmd, FALSE);

	// Ensure style is consistent with this tree (we may have a mix of trees with and without
	// branch lengths and/or internal labels)
	if (!p.CurrentTreeHasEdgeLengths() && (style == PHYLOGRAM_CMD))
	{
		menuBar->Enable (PHYLOGRAM_CMD, FALSE);
		toolBar->EnableTool ( PHYLOGRAM_CMD, FALSE );
		style = RECTANGULAR_TREE_CMD;
	}
	else if (p.CurrentTreeHasEdgeLengths())
	{
		menuBar->Enable (PHYLOGRAM_CMD, TRUE);
		toolBar->EnableTool ( PHYLOGRAM_CMD, TRUE );	
	}
	menuBar->Enable (INTERNAL_LABELS_CMD, p.CurrentTreeHasInternallabels());
	menuBar->Check (INTERNAL_LABELS_CMD, mShowInternalLabels );
	toolBar->EnableTool (INTERNAL_LABELS_CMD, p.CurrentTreeHasInternallabels());
	toolBar->ToggleTool (INTERNAL_LABELS_CMD, mShowInternalLabels );
	
	toolBar->ToggleTool( style, TRUE );
	menuBar->Check (TreeStyleCmd, FALSE);
	menuBar->Check (style, TRUE);
	TreeStyleCmd = style;
}


//------------------------------------------------------------------------------
void TView::OnSlant (wxCommandEvent& WXUNUSED(event))
{
	TreeStyleMenu (SLANTED_TREE_CMD);
	OnUpdate (this, NULL); // Ensure we redraw window
}

//------------------------------------------------------------------------------
void TView::OnRectangle (wxCommandEvent& WXUNUSED(event))
{
	TreeStyleMenu (RECTANGULAR_TREE_CMD);
	OnUpdate (this, NULL); // Ensure we redraw window
}

//------------------------------------------------------------------------------
void TView::OnPhylogram (wxCommandEvent& WXUNUSED(event))
{
	TreeStyleMenu (PHYLOGRAM_CMD);
	OnUpdate (this, NULL); // Ensure we redraw window
}

//------------------------------------------------------------------------------
void TView::OnOrderTree (wxCommandEvent& WXUNUSED(event))
{
#ifdef __WXMAC__
	TreeOrderDlg dialog((wxWindow *)GetFrame());
#else
	TreeOrderDlg dialog((wxWindow *)GetMainFrame());
#endif
	dialog.SetTreeOrder(Order);
	dialog.Centre (wxBOTH);
	if (dialog.ShowModal() == wxID_OK)
	{
		// Get data
		Order = dialog.GetTreeOrder();
		OnUpdate (this, NULL); // Ensure we redraw window		
	}
}

//------------------------------------------------------------------------------
void TView::OnInternalLabels (wxCommandEvent& WXUNUSED(event))
{
	mShowInternalLabels = !mShowInternalLabels;

#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	toolBar->ToggleTool( INTERNAL_LABELS_CMD, mShowInternalLabels );

	wxMenuBar* menuBar = frame->GetMenuBar();
	menuBar->Check (INTERNAL_LABELS_CMD, mShowInternalLabels);

	OnUpdate (this, NULL); // Ensure we redraw window
}


//------------------------------------------------------------------------------
void TView::OnLeafFont (wxCommandEvent& WXUNUSED(event))
{
#if 0
	TreeOrderDlg *dlg = new TreeOrderDlg(frame);
	dlg->ShowModal();
#else
	wxFontData data;
	data.SetInitialFont (LeafFont);

#ifdef __WXMAC__
	wxFontDialog dialog((wxWindow *)GetFrame(), data);
#else
	wxFontDialog dialog((wxWindow *)GetMainFrame(), &data);
#endif
	if ( dialog.ShowModal() == wxID_OK )
	{
		wxFontData retData = dialog.GetFontData();
		LeafFont = retData.GetChosenFont();
		OnUpdate (this, NULL); // Ensure we redraw window
	}
#endif
}

//------------------------------------------------------------------------------
void TView::OnChooseTree (wxCommandEvent& WXUNUSED(event))
{
	wxString *choices = new wxString[p.GetNumTrees()];
	
	for (int i = 0; i < p.GetNumTrees(); i++)
	{
#if wxUSE_UNICODE
		char buf[256];
		strcpy (buf, p.GetIthTreeName(i).c_str());
		wchar_t wbuf[256];
		mbstowcs (wbuf, buf, size_t(wbuf));

		std::wstring tname = wbuf;

		// If tree has no name create one from it's place in the profile
		if (tname == L"")
		{
			wchar_t buf[32];
		    swprintf (buf, wcslen(buf), L"%d", int (i+1));
			tname += buf;
		}
		choices[i] = tname.c_str();
#else
		std::string tname = p.GetIthTreeName(i);

		// If tree has no name create one from it's place in the profile
		if (tname == "")
		{
			char buf[32];
			sprintf (buf, "%d", (i+1));
			tname += buf;
		}
		choices[i] = tname.c_str();
#endif
	}
		
    wxSingleChoiceDialog dialog(frame, wxT("Tree names"),
        wxT("Choose a tree"), p.GetNumTrees(), choices);

    dialog.SetSelection(p.GetCurrentTreeNumber());
	
/*

From charles@plessy.org 13 February 2007 14:55:04 GMT

An Ubuntu user found a bug in TreeView X, and a patch was made by
William Alexander Grant, an Ubuntu developper (CCed). I tested it and
confirmed that it efficiently fixes the problem.

The bug is described in this page:
https://launchpad.net/ubuntu/+source/treeviewx/+bug/75137
(actually, the segfaulting option name seems to be "choose tree")

You can find the pach here:
http://patches.ubuntu.com/t/treeviewx/treeviewx_0.5.1-2ubuntu1.patch

Or more simply:

diff -urNad treeviewx-0.5.1~/tview.cpp treeviewx-0.5.1/tview.cpp
--- treeviewx-0.5.1~/tview.cpp	2006-12-17 10:33:21.000000000 +1100
+++ treeviewx-0.5.1/tview.cpp	2006-12-17 10:33:53.000000000 +1100
@@ -742,7 +742,7 @@

     dialog.SetSelection(p.GetCurrentTreeNumber());

-    if (dialog.ShowModal() == wxID_OK)
+    if ((dialog.ShowModal() == wxID_OK) && (p.GetNumTrees() != 0))
     {
 		int j = dialog.GetSelection ();				
 #ifdef __WXMAC__


By the way, do not hesitate to advertise the Debian and Ubuntu packages
on your website. They are more up to date than the RPM, and support SVG.
*/
	if ((dialog.ShowModal() == wxID_OK) && (p.GetNumTrees() != 0))    
	{
		int j = dialog.GetSelection ();				
#ifdef __WXMAC__
		wxToolBar* toolBar = frame->GetToolBar();
#else
		wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
		toolBar->EnableTool (PREVIOUS_TREE_CMD, (j != 0));
		toolBar->EnableTool (NEXT_TREE_CMD, (j < p.GetNumTrees()-1));
		p.SetCurrentTreeNumber (j);
		TreeStyleMenu(TreeStyleCmd);
		ShowTreeName ();
		OnUpdate (this, NULL); // Ensure we redraw window
	}
    
    delete [] choices;
}

#ifndef __WXMOTIF__
// ----------------------------------------------------------------------------
void TView::OnUpdateCopyCommand(wxUpdateUIEvent& event)
{
	wxMenuBar* menuBar = frame->GetMenuBar();
	menuBar->Enable (wxID_COPY, (p.GetNumTrees() > 0));	
	menuBar->Enable (COPY_AS_TEXT_CMD, (p.GetNumTrees() > 0));	
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	toolBar->EnableTool (wxID_COPY, (p.GetNumTrees() > 0));

}

// ----------------------------------------------------------------------------
void TView::OnUpdatePasteCommand(wxUpdateUIEvent& event)
{
	wxMenuBar* menuBar = frame->GetMenuBar();
	menuBar->Enable (wxID_PASTE, (p.GetNumTrees() == 0) && wxTheClipboard->IsSupported( wxDF_TEXT ));
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	toolBar->EnableTool (wxID_PASTE, (p.GetNumTrees() == 0) && wxTheClipboard->IsSupported( wxDF_TEXT ));
}

// ----------------------------------------------------------------------------
void TView::OnCopy (wxCommandEvent& WXUNUSED(event))
{
#if defined(__WXMAC__) || defined(__WXMSW__)

    // Get current tree
    Tree t = p.GetCurrentTree();
    t.Update();

    // Write tree to clipboard as a metafile
    wxMetafileDC dc;
    if (dc.Ok())
    {
        OnDraw (&dc);
        wxMetafile *mf = dc.Close();
        if (mf)
        {
            bool success = mf->SetClipboard((int)(dc.MaxX() + 10), (int)(dc.MaxY() + 10));
            delete mf;
        }
    }
#endif
}

// ----------------------------------------------------------------------------
void TView::OnCopyAsText (wxCommandEvent& WXUNUSED(event))
{
	// Get current tree
	Tree t = p.GetCurrentTree();
	t.Update();
	
	// Write tree to a string stream
	ostrstream f;	
	NewickTreeWriter tw (&t);
	tw.SetStream (&f);
	tw.Write();	
	f << '\0';

	// Store on clipboard
	if (wxTheClipboard->Open())
	{
		// Data objects are held by the clipboard, 
		// so do not delete them in the app.
		string str = f.str();
		wxTheClipboard->SetData( new wxTextDataObject( wxString(str.c_str(), *wxConvCurrent)) );
		wxTheClipboard->Close();
	}
}



// ----------------------------------------------------------------------------
void TView::OnPaste (wxCommandEvent& WXUNUSED(event))
{
	// Paste tree description from clipboard into this window -- duh! need to make this a call to the main frame too

	if (wxTheClipboard->Open())
	{
		if (wxTheClipboard->IsSupported( wxDF_TEXT ))
		{
			wxTextDataObject data;
			wxTheClipboard->GetData( data );
			
			// Temporary stream to hold clipbard text
			istrstream f (data.GetText().mb_str());
			if (p.ReadTrees (f))
			{
				Activate (TRUE);				// ensure menus are set correctly
				OnUpdate (this, NULL);			// ensure tree is drawn
				GetDocument()->Modify (TRUE);	// tell document that content has changed
			}
			
		}  
		wxTheClipboard->Close();
	}
}
#endif

//------------------------------------------------------------------------------
void TView::ZoomToFit (wxCommandEvent& WXUNUSED(event))
{
	magnification = 1;
	Resize();
}


//------------------------------------------------------------------------------
void TView::ZoomIn (wxCommandEvent& WXUNUSED(event))
{
	magnification++;
	Resize();
}

//------------------------------------------------------------------------------
void TView::ZoomOut (wxCommandEvent& WXUNUSED(event))
{
	magnification--;
	Resize();
}

//------------------------------------------------------------------------------
void TView::Resize ()
{
	canvas->SetMagnification (magnification);
	canvas->Resize ();
	OnUpdate (this, NULL); // Ensure we redraw window
}

// ----------------------------------------------------------------------------
void TView::OnUpdateZoomToFitCommand (wxUpdateUIEvent& event)
{
	wxMenuBar* menuBar = frame->GetMenuBar();
	menuBar->Enable (ZOOM_TO_FIT_CMD, (magnification != 1));	
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	//toolBar->EnableTool (ZOOM_TO_FIT_CMD, (magnification != 1));
}


// ----------------------------------------------------------------------------
void TView::OnUpdateZoomInCommand (wxUpdateUIEvent& event)
{
	wxMenuBar* menuBar = frame->GetMenuBar();
	menuBar->Enable (ZOOM_IN_CMD, (magnification < MAX_MAGNIFICATION));	
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	toolBar->EnableTool (ZOOM_IN_CMD, (magnification < MAX_MAGNIFICATION));
}

// ----------------------------------------------------------------------------
void TView::OnUpdateZoomOutCommand (wxUpdateUIEvent& event)
{
	wxMenuBar* menuBar = frame->GetMenuBar();
	menuBar->Enable (ZOOM_OUT_CMD, (magnification > 1));	
#ifdef __WXMAC__
	wxToolBar* toolBar = frame->GetToolBar();
#else
	wxToolBar* toolBar = GetMainFrame()->GetToolBar();
#endif
	toolBar->EnableTool (ZOOM_OUT_CMD, (magnification > 1));
}

		
//------------------------------------------------------------------------------			
// Map commands to methods in TView

BEGIN_EVENT_TABLE(TView, wxView)
    EVT_MENU(NEXT_TREE_CMD, TView::NextTree)
    EVT_MENU(PREVIOUS_TREE_CMD, TView::PreviousTree)
    EVT_MENU(SLANTED_TREE_CMD, TView::OnSlant)
    EVT_MENU(RECTANGULAR_TREE_CMD, TView::OnRectangle)
    EVT_MENU(PHYLOGRAM_CMD, TView::OnPhylogram)
    EVT_MENU(INTERNAL_LABELS_CMD, TView::OnInternalLabels)
    EVT_MENU(LEAF_FONT_CMD, TView::OnLeafFont)
    EVT_MENU(CHOOSE_TREE_CMD, TView::OnChooseTree)
    EVT_MENU(ORDER_TREE_CMD, TView::OnOrderTree)

#if defined  __WXMSW__ || defined USE_SVG
    EVT_MENU(SAVEAS_PICTURE_CMD, TView::OnSavePicture)
#endif

#ifndef __WXMOTIF__
    EVT_MENU(wxID_COPY, TView::OnCopy)
    EVT_MENU(COPY_AS_TEXT_CMD, TView::OnCopyAsText)
    EVT_MENU(wxID_PASTE, TView::OnPaste)
#endif

    EVT_MENU(ZOOM_TO_FIT_CMD, TView::ZoomToFit)
    EVT_MENU(ZOOM_IN_CMD, TView::ZoomIn)
    EVT_MENU(ZOOM_OUT_CMD, TView::ZoomOut)

        
#ifndef __WXMOTIF__
    EVT_UPDATE_UI(wxID_COPY, TView::OnUpdateCopyCommand)
    EVT_UPDATE_UI(COPY_AS_TEXT_CMD, TView::OnUpdateCopyCommand)
    EVT_UPDATE_UI(wxID_PASTE, TView::OnUpdatePasteCommand)
#endif

    EVT_UPDATE_UI(ZOOM_TO_FIT_CMD, TView::OnUpdateZoomToFitCommand)
    EVT_UPDATE_UI(ZOOM_IN_CMD, TView::OnUpdateZoomInCommand)
    EVT_UPDATE_UI(ZOOM_OUT_CMD, TView::OnUpdateZoomOutCommand)

	
END_EVENT_TABLE()

/*
 * Window implementations
 */

BEGIN_EVENT_TABLE(MyCanvas, wxScrolledWindow)
    EVT_MOUSE_EVENTS(MyCanvas::OnMouseEvent)
	EVT_SIZE(MyCanvas::OnSize)
END_EVENT_TABLE()

// Define a constructor for my canvas
MyCanvas::MyCanvas(wxView *v, wxMDIChildFrame *frame, const wxPoint& pos, const wxSize& size, long style):
 wxScrolledWindow(frame, -1, pos, size, style)
{
	SetBackgroundColour(wxColour(wxT("WHITE")));
	view = v;
	magnification = 1;
}

//------------------------------------------------------------------------------
// Define the repainting behaviour
void MyCanvas::OnDraw(wxDC& dc)
{
	if (view)
		view->OnDraw(& dc);
}

// ----------------------------------------------------------------------------
void MyCanvas::OnSize(wxSizeEvent& event)
{
	Resize ();
    event.Skip();
}


// ----------------------------------------------------------------------------
void MyCanvas::Resize()
{
	int width, height;
	GetClientSize (&width, &height);
	int range = height * magnification; 
	int steps = range / 50;
	SetScrollbars(0, steps, 0, 50);
}


//------------------------------------------------------------------------------
// This implements a tiny doodling program. Drag the mouse using
// the left button.
void MyCanvas::OnMouseEvent(wxMouseEvent& event)
{
}



