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

 // $Id: mast.cpp,v 1.5 2005/05/24 22:55:16 rdmp1c Exp $
 
/**
 * @file mast.cpp
 *
 * Compute Maximum Agreement Subtree (MAST) between a pair of trees
 *
 */

// Vital that this is included before <iomanip.h> or Borland C++ 5.5 complains that
// min and max are alreday defined
#include "mast.h" 
 
#include <iomanip.h>
 
#ifdef __GNUC__
	#include <algo.h>
#endif

#define USE_MATCHING 1
#if USE_MATCHING
	#include "mwbmatching.h"
#endif

typedef std::vector <NNodePtr> NodeVector;

int MAST (NTree &T1, NTree &T2)
{
	int result = 0;

	// 1. create lists of the nodes in T1 and T2 in postorder
    int count = 0;
    NodeVector pot1;
    NodeIterator <Node> n1 (T1.GetRoot());
    Node *q = n1.begin();
    while (q)
    {
    	q->SetIndex (count);
		pot1.push_back ((NNode *)q);
		count++;
		q = n1.next();
    }

    count = 0;
    NodeVector pot2;
    NodeIterator <Node> n2 (T2.GetRoot());
    q = n2.begin();
    while (q)
    {
    	q->SetIndex (count);
		pot2.push_back ((NNode *)q);
		count++;
        q = n2.next();
    }

	// Matrix to hold solutions
    int **m;

    m = new int *[T1.GetNumNodes()];
    for (int i = 0; i < T1.GetNumNodes(); i++)
    	m[i] = new int [T2.GetNumNodes()];

    for (int i = 0; i < T1.GetNumNodes(); i++)
    	for (int j = 0; j <T2.GetNumNodes(); j++)
        	m[i][j] = 0;


	// 2. Visit all pairs of nodes in T1 and T2
     for (int i = 0; i < T1.GetNumNodes(); i++)
     {
    	for (int j = 0; j <T2.GetNumNodes(); j++)
        {
			if (pot1[i]->IsLeaf() || pot2[j]->IsLeaf())
            {
            	// Both are leaves, so MAST[i,j] is 1 if labels are identical
            	if (pot1[i]->IsLeaf() && pot2[j]->IsLeaf())
				{
                	if ( pot1[i]->GetLabel() == pot2[j]->GetLabel())
                    	m[i][j] = 1;
				}
				else
				{
            		// Only one is a leaf, so MAST[i,j] is 1 if leaf is element in cluster
					IntegerSet common;
					std::set_intersection (pot1[i]->Cluster.begin(), pot1[i]->Cluster.end(),
						pot2[j]->Cluster.begin(), pot2[j]->Cluster.end(),
						std::inserter (common, common.end()));
					int w = common.size();
                	m[i][j] = w;
				}
        	}
            else
            {
            	// Both are internals so MAST[i,j] is MAX (diag, match)
                std::vector <NodePtr> pchildren, qchildren;

                // diag
                int diag = 0;

                // Get MAST of base of subtree in t1 and immediate children of
                // base of subtree in t2, and at the same time store list of
                // immediate children
                NodePtr p = pot2[j]->GetChild();
                while (p)
                {
                	qchildren.push_back(p);
                	diag = std::max (diag, m[i][p->GetIndex()]);
                    p = p->GetSibling();
                }
                // get MAST of base of subtree in t2 and immediate children of
                // base of subtree in t1, and at the same time store list of
                // immediate children
                NodePtr q = pot1[i]->GetChild();
                while (q)
                {
                	pchildren.push_back(q);
                	diag = std::max (diag, m[q->GetIndex()][j]);
                    q = q->GetSibling();
                }

                // maximum weighted bipartite matching
				
#if USE_MATCHING
				int match = 0;
				graph g;
				g.make_directed();
				
				// Nodes for p and q children
				map <int, node, less <int> > p_node;
				map <int, node, less <int> > q_node;
				for (int k = 0; k < pchildren.size(); k++)
                {
					node n = g.new_node();
					p_node[k] = n;
                }
				for (int k = 0; k < qchildren.size(); k++)
                {
					node n = g.new_node();
					q_node[k] = n;
                }
				// Edges
				edge_map<int> weights(g, 0);

				
				for (int k = 0; k < pchildren.size(); k++)
				{
					for (int r = 0; r < qchildren.size(); r++)
					{
						int v = pchildren[k]->GetIndex();
                        int w = qchildren[r]->GetIndex();
						
						// It seems that if the partition "from" is much larger than "to,
						// the matching algorithm can blow up with a memory access error
						// in fheap.c. Reversing the graph seems to help.
						
						edge e;
						if (pchildren.size() < qchildren.size())
							e = g.new_edge (p_node[k], q_node[r]);
						else
							e = g.new_edge (q_node[r], p_node[k]);
						
						weights[e] = m[v][w];
					}
				}
				
//				cout << "g" << endl;
//				cout << g;
//				g.save();
//				cout << "Start matching...";
				
				if (g.number_of_nodes() == 0)
				{
					match = 0;
				}
				
				else
				{
				
				mwbmatching mwbm;
				mwbm.set_vars (weights);
				
				if (mwbm.check(g) != algorithm::GTL_OK)
				{
					cout << "Maximum weight bipartite matching algorithm check failed" << endl;
					exit(1);
				}
				else
				{
					if (mwbm.run(g) != algorithm::GTL_OK)
					{
						cout << "Error running maximum weight bipartite matching algorithm" << endl;
						exit(1);
					}
					else
					{
						match = mwbm.get_mwbm();
					}
				}
				}
//				cout << "matching done (" << match << ")" << endl;
#else
                // For now (sigh) brute force generation of all matchings. Eventually
                // will need to do this more efficiently
                int n = std::max (pchildren.size(), qchildren.size());

                // Store a vector of indices of children of subtree in t2.
                // We will permute this to generate all matchings. If t2
                // has fewer children than subtree in t1, vector will contain
                // one or more "null" (-1) values.
                std::vector <int> perm;
                for (int k = 0; k < n; k++)
                {
                	if (k < qchildren.size())
	                	perm.push_back (k);
                    else
                    	perm.push_back (-1);
                }

                // Generate all matchings

                // First matching
                int match = 0;
               	for (int k = 0; k < n; k++)
                {
                    if ((k < pchildren.size()) && (perm[k] != -1))
                    {
                    	int v = pchildren[k]->GetIndex();
                        int w = qchildren[perm[k]]->GetIndex();
                        match += m[v][w];
                    }
                }

                // Remaining matchings
                while  (next_permutation (perm.begin(), perm.end()) )
                {
                	int this_match = 0;
                    for (int k = 0; k < n; k++)
                    {
                        if ((k < pchildren.size()) && (perm[k] != -1))
                        {
                            int v = pchildren[k]->GetIndex();
                            int w = qchildren[perm[k]]->GetIndex();
                            this_match += m[v][w];
                        }
                    }
                    match = std::max (match, this_match);
                }
#endif
                m[i][j] = std::max (diag, match);
            }
		}
     }

     result = m[T1.GetNumNodes() - 1][T2.GetNumNodes() - 1];
	 
	 // Show matrix
/*	 for (int i = 0; i < T1.GetNumNodes(); i++)
	 {
	 	cout << setw(3) << i << "|";
		for (int j = 0; j < T2.GetNumNodes(); j++)
			cout << setw(3) << m[i][j];
		cout << endl;
	}
*/

     // clean up
    for (int i = 0; i < T1.GetNumNodes(); i++)
       	delete [] m[i];

    delete [] m;

	return result;
}

