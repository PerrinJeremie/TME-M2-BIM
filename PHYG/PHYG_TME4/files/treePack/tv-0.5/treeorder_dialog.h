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
 
// $Id: treeorder_dialog.h,v 1.2 2003/10/22 16:46:26 rdmp1c Exp $

#ifndef TREEORDER_DIALOGH
#define TREEORDER_DIALOGH

#include <wx/wx.h>

class TreeOrderDlg : public wxDialog
{
public:
	TreeOrderDlg (wxWindow *parent);
	int GetTreeOrder () { return radio_buttons->GetSelection(); };
	void SetTreeOrder (int n) { radio_buttons->SetSelection(n); };
private:
	wxRadioBox *radio_buttons;
};

#endif
