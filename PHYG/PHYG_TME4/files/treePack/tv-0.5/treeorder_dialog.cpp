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
 
// $Id: treeorder_dialog.cpp,v 1.2 2005/02/21 19:41:24 rdmp1c Exp $

#include "treeorder_dialog.h"

TreeOrderDlg::TreeOrderDlg (wxWindow *parent): wxDialog (parent, -1, wxT("Tree Order"),wxDefaultPosition, wxSize(185,185))
{
	// Top sizer
	wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );
	
	// Radiobuttons	
	wxString *choices = new wxString[4];
	choices[0] = wxT("Default");
	choices[1] = wxT("Alphabetical");
	choices[2] = wxT("\"Ladderise\" left");
	choices[3] = wxT("\"Ladderise\" right");

	radio_buttons = new wxRadioBox (this, -1, wxT("Order"), 
		wxDefaultPosition, wxDefaultSize, 4, choices, 1, wxRA_SPECIFY_COLS);
	topsizer->Add (radio_buttons,
		0,
		wxALIGN_CENTER);

	// Bottom panel has OK and Cancel buttons
    wxBoxSizer *okcancel_sizer = new wxBoxSizer( wxHORIZONTAL );
	wxButton *ok_button = new wxButton( this, wxID_OK, wxT("OK") );
	ok_button->SetDefault();
	okcancel_sizer->Add (ok_button,
		0,
		wxALL,
		10 );

	okcancel_sizer->Add(
		new wxButton( this, wxID_CANCEL, wxT("Cancel") ),
		0,
		wxALL,
		10 );
	 
	topsizer->Add(
		okcancel_sizer,
		0,
		wxALIGN_CENTER );
	
	SetAutoLayout( TRUE );
	SetSizer( topsizer );
	
	topsizer->Fit( this );
	topsizer->SetSizeHints( this );
}

