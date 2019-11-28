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
 
 // $Id: tproject.h,v 1.9 2005/08/31 08:55:25 rdmp1c Exp $
 
#ifndef TPROJECTH
#define TPROJECTH

#include "profile.h"

/**
 * @class TProject
 * Encapsulates a profile of trees and operations on those trees. This is basically
 * a wrapper around the Profile class, and is designed for use in GUI programs.
 *
 * @see Profile
 *
 */

template <class T> class TProject
{
public:
	/**
	 * Constructor
	 */	
	TProject () { CurrentTree = 0; };
	/**
	 * Destructor
	 */	
	virtual ~TProject () {};
	/**
	 * True if the current tree has edge lengths
	 */	
	virtual bool CurrentTreeHasEdgeLengths () { return GetCurrentTree().GetHasEdgeLengths(); };
	/**
	 * True if the current tree has internal node labels
	 */	
	virtual bool CurrentTreeHasInternallabels () { return GetCurrentTree().GetHasInternalLabels(); };
	/**
	 * Get the current tree object
	 * @return The tree
	 */	

	virtual T GetCurrentTree () { return p.GetIthTree (CurrentTree); };
	/**
	 * Get the index of the current tree
	 * @return The index
	 */	
	virtual int GetCurrentTreeNumber () const { return CurrentTree; };
	/**
	 * @return The number of trees in the profile
	 */	
	virtual int GetNumTrees() { return p.GetNumTrees (); };
	/**
	 * @brief The name of the ith tree in the profile
	 *
	 * @param i the index of the tree in the range 0 - (n-1)
	 * @return The tree
	 */
	virtual std::string GetIthTreeName (int i) { return p.GetIthTreeName (i); };
	/**
	 * Read one or more trees from an input stream
	 * @param f an input stream
	 * @return True if trees read in successfully
	 */	
	virtual bool ReadTrees (std::istream &f) { return p.ReadTrees (f); };
	/**
	 * Set the index of the current tree
	 * @param i the index of the tree in the range 0 - (n-1)
	 */	
	virtual void SetCurrentTreeNumber (const int i) { CurrentTree = i; };
	/**
	 * Write trees to output stream
	 * @param f an output stream
	 * @return True if trees written successfully
	 */	
	virtual bool WriteTrees (std::ostream &f, const int format = 0, const char *endOfLine = "\n") 
	{ return p.WriteTrees (f, format, endOfLine); };

protected:
	/**
	 * The index of the current tree
	 */
	int CurrentTree;
	/**
	 * The profile of trees
	 */	
	Profile <T> p;
};



#endif
