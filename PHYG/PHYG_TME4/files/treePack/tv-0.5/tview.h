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
 
 // $Id: tview.h,v 1.21 2005/08/31 08:55:01 rdmp1c Exp $
 
#ifdef __GNUG__
// #pragma interface
#endif

#include "tproject.h"

#ifndef __VIEWSAMPLEH__
#define __VIEWSAMPLEH__

#include "wx/docview.h"

#define PREVIOUS_TREE_CMD		100
#define NEXT_TREE_CMD			101
#define SLANTED_TREE_CMD		102
#define RECTANGULAR_TREE_CMD		103
#define PHYLOGRAM_CMD			104
#define LEAF_FONT_CMD			105
#define CHOOSE_TREE_CMD			106
#define ORDER_TREE_CMD			107
#define SAVEAS_PICTURE_CMD		108
#define COPY_AS_TEXT_CMD		109

#define ZOOM_IN_CMD				110
#define ZOOM_OUT_CMD				111
#define ZOOM_TO_FIT_CMD			112

#define INTERNAL_LABELS_CMD	113

const int DEFAULT_ORDER = 0;
const int ALPHA_ORDER 	= 1;
const int LEFT_ORDER 	= 2;
const int RIGHT_ORDER 	= 3;

#define MAX_MAGNIFICATION 5

class MyCanvas: public wxScrolledWindow
{
public:
	wxView *view;

	MyCanvas(wxView *v, wxMDIChildFrame *frame, const wxPoint& pos, const wxSize& size, long style);
	virtual void OnDraw(wxDC& dc);
	virtual void OnSize(wxSizeEvent& event);
	void OnMouseEvent(wxMouseEvent& event);
	void Resize ();
	void SetMagnification (int m) { magnification = m; };
	
protected:
	int magnification;

private:
	DECLARE_EVENT_TABLE()
};


class TView: public wxView
{
public:
	wxMDIChildFrame *frame;
	MyCanvas *canvas;

	TView();
	~TView() {}

	void OnActivateView (bool activate, wxView *activeView, wxView *deactiveView);
	bool OnCreate(wxDocument *doc, long flags);
	void OnDraw(wxDC *dc);
	void OnUpdate(wxView *sender, wxObject *hint = (wxObject *) NULL);
	bool OnClose(bool deleteWindow = TRUE);
	virtual wxPrintout *OnCreatePrintout();

	void OnSlant (wxCommandEvent& WXUNUSED(event) ); 
	void OnRectangle (wxCommandEvent& WXUNUSED(event) ); 
	void OnPhylogram (wxCommandEvent& WXUNUSED(event) ); 
	
	void OnInternalLabels (wxCommandEvent& WXUNUSED(event) );

	void OnLeafFont (wxCommandEvent& WXUNUSED(event) );

	void OnChooseTree (wxCommandEvent& WXUNUSED(event) ); 

	void OnOrderTree (wxCommandEvent& WXUNUSED(event) ); 

	void OnSavePicture (wxCommandEvent& WXUNUSED(event) ); 
	
	void NextTree (wxCommandEvent& WXUNUSED(event) );    
	void PreviousTree (wxCommandEvent& WXUNUSED(event) );
	void ShowTreeName ();
	
	void TreeStyleMenu (int style);
	
	// Clipboard
#ifndef __WXMOTIF__
	void OnUpdateCopyCommand(wxUpdateUIEvent& event);
	void OnUpdatePasteCommand(wxUpdateUIEvent& event);
	void OnCopy (wxCommandEvent& event);
	void OnCopyAsText (wxCommandEvent& event);
	void OnPaste(wxCommandEvent& event);
#endif

	// Zoom
	void ZoomToFit (wxCommandEvent& WXUNUSED(event) );    
	void ZoomIn (wxCommandEvent& WXUNUSED(event) );    
	void ZoomOut (wxCommandEvent& WXUNUSED(event) );
	void Resize ();
	void OnUpdateZoomToFitCommand (wxUpdateUIEvent& event);
	void OnUpdateZoomInCommand (wxUpdateUIEvent& event);
	void OnUpdateZoomOutCommand (wxUpdateUIEvent& event);

	TProject <Tree>	p;
	unsigned int	TreeFlags; 
	int		TreeStyleCmd;
	wxFont		LeafFont;
	bool		mShowInternalLabels;

	int 		page_width;
	int 		page_height;
	bool 		printing;
	
	int 		magnification;

	bool ReadTrees (const wxString &filename);
	bool WriteTrees (const wxString &filename);
	
protected:
	int Order;

private:
	DECLARE_DYNAMIC_CLASS(TView)
	DECLARE_EVENT_TABLE()
};


#endif
