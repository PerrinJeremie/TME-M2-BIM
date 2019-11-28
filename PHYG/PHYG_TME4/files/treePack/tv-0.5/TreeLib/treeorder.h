/*
 * TreeLib
 * A library for manipulating phylogenetic trees.
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

 // $Id: treeorder.h,v 1.4 2005/02/23 14:19:00 rdmp1c Exp $
 
/**
 * @file treeorder.h
 *
 * Classes to reorder nodes in a tree
 *
 */

#ifndef TREEORDERH
#define TREEORDERH

#ifdef __BORLANDC__
	// Undefine __MINMAX_DEFINED so that min and max are correctly defined
	#ifdef __MINMAX_DEFINED
		#undef __MINMAX_DEFINED
	#endif
    // Ignore "Cannot create precompiled header: code in header" message
    // generated when compiling string.cc
    #pragma warn -pch
#endif

#include <map>

#include "TreeLib.h"
#include "nodeiterator.h"


/**
 * @class TreeOrder
 * Base class for ordering tree. By overiding the MustSwap fucntion, any
 * ordering can be described.
 */
class TreeOrder
{
public:
	TreeOrder () {};
	/**
	 * Constructor takes a tree as a parameter.
     * @param tree tree to be ordered.
     */
	TreeOrder (Tree *tree) { t = tree; };
        virtual ~TreeOrder() {};
 	/**
	 * Order the tree. Call this method to reorder nodes. The tree is traversed
     * in post order, and for each internal node SortDescendants is called. How the
     * nodes are sorted is determined by MustSwap, so normally Order will not need
     * to be overridden in descendant classes.
	 */
	virtual void Order ();
protected:
 	Tree *t;
 	/**
	 * Test whether nodes p and q need to be swapped. This is an abstract
     * function that must be overriden in descendant classes. This function
     * determines how the tree is ordered.
     * @param p a node in the tree
     * @param q another node in the tree that is a sibling of q
	 * @return True if p and q must be swapped.
	 */
	virtual bool MustSwap (NodePtr p, NodePtr q) = 0;
 	/**
	 * Sort descendants of node according to criterion defined in MustSwap.
     * @param node the node whose descendants are being sorted.
	 */
	virtual void SortDescendants (NodePtr node);
};



/**
 * @class LeftOrder
 * Extends Treeorder to order trees by placing "heavier" nodes to the left.
 */
class LeftOrder : public  TreeOrder
{
public:
	LeftOrder (Tree *tree) : TreeOrder (tree) {};
        virtual ~LeftOrder() {};
protected:
	virtual bool MustSwap (NodePtr p, NodePtr q);

};

/**
 * @class RightOrder
 * Extends Treeorder to order trees by placing "heavier" nodes to the right.
 */
class RightOrder : public  TreeOrder
{
public:
	RightOrder (Tree *tree) : TreeOrder (tree) {};
        virtual ~RightOrder() {};
protected:
	virtual bool MustSwap (NodePtr p, NodePtr q);

};


/**
 * @class AlphaOrder
 * Extends Treeorder to order trees by sorting nodes alphabetically.
 */
class AlphaOrder : public  TreeOrder
{
public:
	AlphaOrder (Tree *tree) : TreeOrder (tree) {};
        virtual ~AlphaOrder() {};
	/**
	 * Order the tree. Extends ancestral method by storing "lowest" label
     * found in any descendant of a given node using a map.
	 */
    virtual void Order ();
protected:
	std::map<Node *, std::string, std::less<Node *> > labels;
	virtual bool MustSwap (NodePtr p, NodePtr q);

};

/**
 * @class ThorneOrder
 * Extends Treeorder to order in form required by Jeff Thorne's programs
 * estabranches and multidivtime. That is, the outgroup is the rightmost
 * node of the tree.
 */
class ThorneOrder : public  TreeOrder
{
public:
	ThorneOrder (Tree *tree) : TreeOrder (tree) { outgroup = t->GetRoot()->GetChild(); };
	virtual ~ThorneOrder() {};
	virtual void SetOutgroup (NodePtr o) { outgroup = o; };
	/**
	 * Order the tree. Extends ancestral method by storing "lowest" label
     * found in any descendant of a given node using a map.
	 */
    virtual void Order ();
protected:
	virtual bool MustSwap (NodePtr p, NodePtr q);
	NodePtr outgroup;

};


#if __BORLANDC__
	// Redefine __MINMAX_DEFINED so Windows header files compile
	#ifndef __MINMAX_DEFINED
    		#define __MINMAX_DEFINED
	#endif
#endif


#endif
