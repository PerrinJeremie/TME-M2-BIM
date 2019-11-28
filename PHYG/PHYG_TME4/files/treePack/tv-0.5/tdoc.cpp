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
 
 // $Id: tdoc.cpp,v 1.6 2001/07/25 10:37:06 rdmp1c Exp $
 
#ifdef __GNUG__
// #pragma implementation
#endif

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


IMPLEMENT_DYNAMIC_CLASS(TDocument, wxDocument)

TDocument::TDocument(void)
{
}

TDocument::~TDocument(void)
{
}




// Since text windows have their own method for saving to/loading from files,
// we override OnSave/OpenDocument instead of Save/LoadObject
bool TDocument::OnSaveDocument(const wxString& filename)
{
	TView *view = (TView *)GetFirstView();
	if (!view->WriteTrees (filename))
        	return FALSE;
	Modify(FALSE);
	return TRUE;
}

bool TDocument::OnOpenDocument(const wxString& filename)
{
	TView *view = (TView *)GetFirstView();
	if (!view->ReadTrees (filename))
		return FALSE;

	SetFilename(filename, TRUE);

	Modify(FALSE);
	UpdateAllViews();
	view->Activate(TRUE); // Ensures TView::OnActivate is called  
	return TRUE;
}


