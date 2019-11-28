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
 
// $Id: treewriter.h,v 1.5 2005/02/23 14:36:08 rdmp1c Exp $

#ifndef TREEWRITER_H
#define TREEWRITER_H

#include "TreeLib.h"
#include "nodeiterator.h"

#include <iostream>

/**
 * @class NewickTreeWriter
 * Base class for writing Newick format tree descriptions.
 *
 */
class NewickTreeWriter
{
public:
	NewickTreeWriter (Tree *tree) { t = tree; writeEdgeLengths = true; };
        virtual ~NewickTreeWriter () {};
	
	/**
	 * Toggle on or off the writing of edge lengths in the tree description.
	 */	
	virtual void SetWriteEdgeLengths (bool on) { writeEdgeLengths = on; };
	/**
	 * Set output stream to which the tree description is written
	 * @param s pointer to the stream
	 */	
	virtual void SetStream (std::ostream *s) { f = s; };
	/**
	 * Write the tree description.
	 */	    
    virtual void Write ();
protected:
	std::ostream *f;
	Node *cur;
	std::stack < Node *, std::vector<Node *> > stk;
	Tree * t;
	
	/**
	 * Write the symbol signalling the end of the tree description. By default this is ';'
	 */	
	virtual void WriteEndOfTree ();
	/**
	 * Write '('
	 */	
    virtual void WriteLeftParenthesis ();
	/**
	 * Write ')'
	 */	
    virtual void WriteRightParenthesis ();
	/**
	 * Write the sibling symbol. By default this is ','
	 */	
    virtual void WriteSiblingSymbol ();
	/**
	 * Write leaf node.
	 */	
    virtual void WriteLeaf ();
	/**
	 * Write internal node.
	 */	
    virtual void WriteInternal ();
    
    bool writeEdgeLengths;

};

/**
 * @class XMLTreeWriter
 * Abstract base class for writing XML tree descriptions.
 *
 */
class XMLTreeWriter
{
public:
	XMLTreeWriter (Tree *tree) { t = tree; count = 0; };
        virtual ~XMLTreeWriter () {};
	
	/**
	 * Set output stream to which the tree description is written
	 * @param s pointer to the stream
	 */	
	virtual void SetStream (std::ostream *s) { f = s; };
	/**
	 * Write the tree description.
	 */	    
        virtual void Write ();
        
        virtual void OpenNode() = 0;
        virtual void NodeInfo () = 0;
        virtual void CloseNode() = 0;
        virtual void OpenTree() = 0;
        virtual void CloseTree() = 0;

protected:
	std::ostream *f;
	Node *cur;
	std::stack < Node *, std::vector<Node *> > stk;
	Tree * t;
        int count;

};

/**
 * @class TreeBolicTreeWriter
 * Writing TreeBolic XML tree description.
 *
 */
class TreeBolicTreeWriter : public XMLTreeWriter
{
public:
	TreeBolicTreeWriter (Tree *tree) : XMLTreeWriter (tree) {};
        virtual ~TreeBolicTreeWriter () {};
        virtual void OpenNode();
        virtual void NodeInfo ();
        virtual void CloseNode();
        virtual void OpenTree();
        virtual void CloseTree();
};


/**
 * @class ThorneTreeWriter
 * Writing tree for Jeff Thorne's multidivtime tree description.
 *
 */
class ThorneTreeWriter : public NewickTreeWriter
{
public:
	ThorneTreeWriter (Tree *tree) : NewickTreeWriter (tree) {};
	virtual ~ThorneTreeWriter () {};
	/**
	 * Write leaf node.
	 */	
    virtual void WriteLeaf ();
	/**
	 * Write internal node.
	 */	
    virtual void WriteInternal ();
protected:
	std::string formatLabel (std::string s);
};



#endif


