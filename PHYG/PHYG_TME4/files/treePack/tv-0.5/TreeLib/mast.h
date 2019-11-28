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

 // $Id: mast.h,v 1.1 2002/05/16 16:07:32 rdmp1c Exp $
 
/**
 * @file mast.h
 *
 * Compute Maximum Agreement Subtree (MAST) between a pair of trees
 *
 */
#ifndef MASTH
#define MASTH

#ifdef __BORLANDC__
	// Undefine __MINMAX_DEFINED so that min and max are correctly defined
	#ifdef __MINMAX_DEFINED
		#undef __MINMAX_DEFINED
	#endif
    // Ignore "Cannot create precompiled header: code in header" message
    // generated when compiling string.cc
    #pragma warn -pch
#endif

#include "nodeiterator.h"
#include "ntree.h"

/**
 * @fn int MAST (NTree &T1, NTree &T2)
 * @brief Computes the size of the maximum agreement subtree between two rooted trees
 *
 * MAST uses the dynamic programming algorithm of Steel and Warnow (1993) to compute
 * the size of the largest agreement subtree. The algorithm has time complexity
 * @f${\rm O}(n^{4.5} \log n + V)@f$ for two trees with @em n leaves, where @em V
 * is the maximum number of nodes in the trees. This implementation would be much improved
 * by using an efficient algorithm for maximum weight bipartite matching (rather than
 * simply enumerating all @f$(k - 1)!@f$ possible matchings), and by using the
 * sparsification techniques of Farach and Thorup (1997).
 *
 * @param T1 one tree
 * @param T2 another tree
 */
int MAST (NTree &T1, NTree &T2);

#if __BORLANDC__
	// Redefine __MINMAX_DEFINED so Windows header files compile
	#ifndef __MINMAX_DEFINED
    		#define __MINMAX_DEFINED
	#endif
#endif

#endif
