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
 
 // $Id: tv.h,v 1.4 2004/08/05 15:55:26 rdmp1c Exp $
 
 #ifdef __GNUG__
// #pragma interface "docview.h"
#endif

#ifndef __DOCVIEWSAMPLEH__
#define __DOCVIEWSAMPLEH__

#include "wx/docmdi.h"

class wxDocManager;

// Define a new application
class MyApp: public wxApp
{
public:
	MyApp(void);
	bool OnInit(void);
	int OnExit(void);

	wxMDIChildFrame *CreateChildFrame(wxDocument *doc, wxView *view, bool isCanvas);

#ifdef __WXMAC__
    virtual void MacOpenFile(const wxString &fileName);
#endif

protected:
	wxDocManager* m_docManager;
};

DECLARE_APP(MyApp)

// Define a new frame
class MyCanvas;
class MyFrame: public wxDocMDIParentFrame
{
	DECLARE_CLASS(MyFrame)
 public:
//  wxMenu *editMenu;
  
	MyFrame(wxDocManager *manager, wxFrame *frame, const wxString& title, const wxPoint& pos, const wxSize& size,
		long type);
	void InitToolBar(wxToolBar* toolBar);
	void OnAbout(wxCommandEvent& event);
	MyCanvas *CreateCanvas(wxView *view, wxMDIChildFrame *parent);

DECLARE_EVENT_TABLE()
};

extern MyFrame *GetMainFrame(void);

#define DOCVIEW_CUT     1
#define DOCVIEW_ABOUT   2

extern bool singleWindowMode;

#endif
