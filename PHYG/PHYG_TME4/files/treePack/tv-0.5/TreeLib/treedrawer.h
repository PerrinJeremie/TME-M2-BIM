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

// $Id: treedrawer.h,v 1.10 2007/06/07 15:48:02 rdmp1c Exp $

 /**
 * @file treedrawer.h
 *
 * Classes to draw trees to a graphics device
 *
 */

#ifndef TREEDRAWERH
#define TREEDRAWERH

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
#include "gport.h"

#if USE_VC2
	#define USE_PORT 0
	#include "VPort.h"
#elif USE_WXWINDOWS
	#define USE_PORT 0
   	#include "wx/wx.h"
	#ifdef wxSTRING_MAXLEN
	#else
		#define wxSTRING_MAXLEN 255
	#endif

#else
	#define USE_PORT 1
#endif

/**
 * @class point
 * An x, y point in real number coodinates
 */
class point
{
public:
    double x;
    double y;
    point () { x = y = 0.0; };
};

/**
 * @class TreeDrawer
 * Base class for tree drawer object
 */
class TreeDrawer
{
public:
    TreeDrawer (Tree *tree);
    virtual ~TreeDrawer () {};
    /**
     * Calculate coordinates for internal nodes. For a slanted cladogram
     * internal node are spaced along the x-axis by their "weight", and along the
     * y-axis they are placed in the middle of their descendant leaves
     * @param p the internal node being visited
     */
    virtual void CalcInternal (Node *p);
    /**
     * Calculate coordinates for leaf nodes. For a slanted cladogram
     * all leaves have the same coordinate on the x-aixs, and are spaced
     * evenly along the y-axis
     * @param p the leaf being visited
     */
    virtual void CalcLeaf (Node *p);
 	/**
	 * Calculate coordinates for tree drawing. Tree is traversed in post order.
	 */
    virtual void CalcCoordinates ();
 	/**
	 * Draw tree by traversing it in post-order and calling DrawLeaf and
     * DrawInternal.
     */
    virtual void Draw();
 	/**
	 * Draw a line connecting the leaf with its immediate ancestor.
     */
    virtual void DrawLeaf (Node *p);
 	/**
	 * Draw leaf label (if showLeafLabels is true)
     */
	virtual void DrawLeafLabel (Node *p);
 	/**
	 * Draw a line connecting internal node with its immediate ancestor.
     */
    virtual void DrawInternal (Node *p);
 	/**
	 * Draw intenral label (if it exists and if showInternalLabels is true)
     */
	virtual void DrawInternalLabel (Node *p);
    /**
     * Draw a line to the graphics device connecting the two points pt1 and pt2.
     * @param pt1 start of the line
     * @param pt2 end of the line
     */
    virtual void DrawLine (point pt1, point pt2);
    /**
     * Draw a line representing the root of the tree.
     */
    virtual void DrawRoot ();
    /**
     * Draw a text to the graphics device.
     * @param pt posuition for text
     * @param s text to draw
     */
    virtual void DrawText (point pt, std::string s);
	
	virtual point GetCoordinates (NodePtr p) { return node_coordinates[p]; }

    virtual void SetRect (int l, int t, int right, int bottom)
    {
    	left = (double) l;
        top = (double) t;
        width = (double)(right - l);
        height = (double)(bottom - t);
    }
    
#if USE_WXWINDOWS
	virtual void SetDC (wxDC *wxdc) { dc = wxdc; };
#endif

#if USE_PORT
    virtual void SetFont (GBaseFont f) { font = f; };
#endif
    /**
     * Toggle drawing of tree as rooted. For a cladogram this toggles the display
     * of a short branch below the root node.
     */
    virtual void SetDrawRooted (bool on = true) { rooted = on; };
    /**
     * Toggle drawing of internal labels.
     */
    virtual void SetDrawInternalLabels (bool on = true) { showInternalLabels = on; };
    /**
     * Toggle drawing of leaf labels.
     */
    virtual void SetDrawLeafLabels (bool on = true) { showLeafLabels = on; };
    
    virtual void SetPenWidth (int width);

protected:
    Tree * t;
    std::map<Node *, point, std::less<Node *> > node_coordinates;
    double left, top, width, height;
    double leafGap, lastY, nodeGap;
    int leafCount;
    bool rooted;
    bool showInternalLabels;
    bool showLeafLabels;

#if USE_PORT
    GBaseFont font;
#endif

#if USE_WXWINDOWS
    wxDC	*dc;
#endif

};

/**
 * @class RectangleTreeDrawer
 * Draw rectangular trees
 */
class RectangleTreeDrawer : public TreeDrawer
{
public:
    RectangleTreeDrawer (Tree *tree) : TreeDrawer (tree) {};
    virtual ~RectangleTreeDrawer () {};

    /**
     * Calculate coordinates for internal nodes. For a rectangular cladogram we
     * find the longest path from the root to a leaf, and space the internal nodes
     * evenly. This gives a nicer look than using weights.
     * @param p the internal node being visited
     */
    virtual void CalcInternal (Node *p);
    virtual void CalcCoordinates ();
    virtual void DrawLeaf (Node *p);
    virtual void DrawInternal (Node *p);
protected:
    int maxDepth;

};

/**
 * @class PhylogramDrawer
 * Draw phylograms
 */
class PhylogramDrawer : public RectangleTreeDrawer
{
public:
    PhylogramDrawer (Tree *tree) : RectangleTreeDrawer (tree) { mUltrametric = false; };
    virtual ~PhylogramDrawer () {};
    /**
     * Calculate coordinates for internal nodes.
     * @param p the internal node being visited
     */
    virtual void CalcInternal (Node *p);
    virtual void CalcLeaf (Node *p);
    virtual void CalcCoordinates ();
    /**
     * Draw tree by traversing it in post-order and calling DrawLeaf and
     * DrawInternal.
     */
    virtual void Draw();
//    virtual void DrawLeaf (Node *p);
//    virtual void DrawInternal (Node *p);
    virtual void DrawRoot ();
    virtual void DrawScaleBar ();
protected:
	double mMaxPathLength;
	int scalebar_space;
    bool mUltrametric;
};

/**
 * @class CircleTreeDrawer
 * Draw circle trees
 */
class CircleTreeDrawer  : public RectangleTreeDrawer
{
public:
    CircleTreeDrawer (Tree *tree) : RectangleTreeDrawer (tree) {};
    virtual ~CircleTreeDrawer () {};

    virtual void CalcInternal (Node *p);
    virtual void CalcLeaf (Node *p);
    virtual void CalcCoordinates ();

//    virtual void Draw();

    virtual void DrawLeaf (Node *p);
//	virtual void DrawLeafLabel (Node *p);
    virtual void DrawInternal (Node *p);
//	virtual void DrawInternalLabel (Node *p);
//   virtual void DrawText (point pt, std::string s);

protected:
    double leaf_angle;
    double leaf_radius;
    point origin;
    std::map<Node *, double, std::less<Node *> > node_angle;
    std::map<Node *, double, std::less<Node *> > node_radius;
    std::map<Node *, point, std::less<Node *> > node_backarc;
};

/**
 * @class KMLTreeDrawer
 * Draw KML trees for Google Earth
 */
 
class KMLTreeDrawer  : public RectangleTreeDrawer
{
public:
    KMLTreeDrawer (Tree *tree) : RectangleTreeDrawer (tree) { altitudeFactor = 20000.0; treeStyle = 0; };
    virtual ~KMLTreeDrawer () {};

    virtual void CalcInternal (Node *p);
    virtual void CalcLeaf (Node *p);
    virtual void CalcCoordinates ();

    virtual void Draw();

    virtual void DrawLeaf (Node *p);
	virtual void DrawLeafLabel (Node *p);
    virtual void DrawInternal (Node *p);
	virtual void DrawInternalLabel (Node *p);
    virtual void DrawRoot ();

//   virtual void DrawText (point pt, std::string s);

	/**
	 * Set output stream to which the tree description is written
	 * @param s pointer to the stream
	 */	
	virtual void SetStream (std::ostream *s) { f = s; };

	virtual void SetStyle (int s) { treeStyle=s; };

protected:
    std::map<Node *, double, std::less<Node *> > nodeLatitude;
    std::map<Node *, double, std::less<Node *> > nodeLongitude;
    std::map<Node *, double, std::less<Node *> > nodeAltitude;
	
	std::ostream *f;
	Node *cur;
	std::stack < Node *, std::vector<Node *> > stk;
	int count;
	double altitudeFactor;
   	double mMaxPathLength;
	
	int treeStyle;
};



#if __BORLANDC__
	// Redefine __MINMAX_DEFINED so Windows header files compile
	#ifndef __MINMAX_DEFINED
    		#define __MINMAX_DEFINED
	#endif
#endif

#endif
