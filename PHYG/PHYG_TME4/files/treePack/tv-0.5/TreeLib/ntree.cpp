#include "ntree.h"

#ifdef __GNUC__
#include <algo.h>
#endif

#include <iostream>

enum SetRelations {IDENTITY = 0, DISJOINT, SUBSET, SUPERSET, OVERLAPPING};


SetRelations Relationship (IntegerSet &one, IntegerSet &two);


SetRelations Relationship (IntegerSet &one, IntegerSet &two)
{
	SetRelations result;

/*        cout << " {";
        std::copy (one.begin(), one.end(),
                std::ostream_iterator<int>(cout, " "));
        cout << "}";

       cout << " {";
        std::copy (two.begin(), two.end(),
                std::ostream_iterator<int>(cout, " "));
        cout << "}";
*/
    if (one == two)
    	result = IDENTITY;
    else if (std::includes (two.begin(), two.end(), one.begin(), one.end()))
    	result = SUBSET; // one is a subset of two
    else if (std::includes (one.begin(), one.end(), two.begin(), two.end()))
    	result = SUPERSET; // two is a subset of one
    else
    {
        IntegerSet common;
        std::set_intersection (one.begin(), one.end(),
            two.begin(), two.end(),
            std::inserter (common, common.end()));
        if (common.size() == 0)
        	result = DISJOINT;
        else
        	result = OVERLAPPING;
    }

//    cout << "result=" << result << endl;

    return result;
}


void NTree::BuildClustersTraverse (NNodePtr p)
{
	if (p)
    {
    	p->Cluster.erase (p->Cluster.begin(), p->Cluster.end());
        BuildClustersTraverse ((NNodePtr)(p->GetChild()));
        BuildClustersTraverse ((NNodePtr)(p->GetSibling()));
        if (p->IsLeaf())
        {
        	switch (use)
            {
            	case useLeafNumber:
        			p->Cluster.insert (p->GetLeafNumber());
                    break;

                case useLabelNumber:
        			p->Cluster.insert (p->GetLabelNumber());
                    break;
            }
        }
        if (p !=Root)
        {
        	NNodePtr anc = (NNodePtr)(p->GetAnc());
			std::set<int, std::less<int> > temp_set;
			std::set_union(anc->Cluster.begin(), anc->Cluster.end(),
				p->Cluster.begin(), p->Cluster.end(),
				std::inserter(temp_set, temp_set.begin()));
			anc->Cluster.swap(temp_set);
        }
    }
}

void NTree::BuildLeafLabelsTraverse (NNodePtr p)
{
	if (p)
    {
        BuildLeafLabelsTraverse ((NNodePtr)(p->GetChild()));
        BuildLeafLabelsTraverse ((NNodePtr)(p->GetSibling()));
        if (p->IsLeaf())
        {
			leaf_labels[p->GetLabel()] = p;
        }
    }
}

// debugging
//------------------------------------------------------------------------------
void NTree::ShowClustersTraverse (NNodePtr p)
{
	if (p)
	{
        ShowClustersTraverse ((NNodePtr)(p->GetChild()));
        // Show clusters
		std::cout << "{ ";
		std::copy (p->Cluster.begin(), p->Cluster.end(),
			std::ostream_iterator<int>(std::cout, " "));
		std::cout << "}" << std::endl;
        ShowClustersTraverse ((NNodePtr)(p->GetSibling()));
	}
}

//------------------------------------------------------------------------------
// Create a star tree with n leaves
void NTree::StarTree (int n)
{
	Leaves = n;
	Internals = 1;

	Root = NewNode();
	Root->SetWeight (n);
	Root->SetDegree (n);

	CurNode = NewNode();
	CurNode->SetLeaf(true);
	CurNode->SetLeafNumber(1);
	CurNode->SetLabelNumber(1);

	Root->SetChild (CurNode);
	CurNode->SetAnc (Root);

	// Remaining leaves
	for (int i = 1; i < n; i++)
	{
		NodePtr q = NewNode ();
		q->SetLeaf(true);
		q->SetLeafNumber(i+1);
		q->SetLabelNumber(i+1);
		q->SetAnc (Root);
		CurNode->SetSibling (q);;
		CurNode = q;
	}

    MakeNodeList();
    Update();
    BuildLeafClusters ();

}

//------------------------------------------------------------------------------
// Add cluster to a tree
int NTree::AddCluster (IntegerSet &b, char *label)
{
	Error = 0;
	int n = b.size();
	if (n > 1 && n < ((NNodePtr)Root)->Cluster.size())
	{
		// non trivial set
		NNodePtr p = (NNodePtr)Root;
		NNodePtr q = (NNodePtr)(p->GetChild());
        bool done = false;
		SetRelations r;

		// locate insertion point
		while (q && !done)
		{
			r = Relationship (b, q->Cluster);

			switch (r)
			{
				case IDENTITY:
					done = true;
					Error = 1; // we already have this cluster
					break;

				case SUBSET:
					p = q;
					q = (NNodePtr)(q->GetChild());
					break;

				case DISJOINT:
					p = q;
					q = (NNodePtr)(q->GetSibling());
					break;

				case SUPERSET:
					done = true;
					break;

				case OVERLAPPING:
					done = true;
					Error = 2; // can't be a part of a tree
					break;

				default:
					break;
			}
		}

		if (Error == 0)
		{
			// Create node
			NNodePtr nunode = (NNodePtr)NewNode ();
			Internals++;
			nunode->SetWeight(n);
			nunode->Cluster = b;

			if (label)
				nunode->SetLabel (label);


			// Insert node in tree
			switch (r)
			{
				case SUBSET:
                	p->SetChild (nunode);
                    nunode->SetAnc (p);
					break;

				case DISJOINT:
					p->SetSibling (nunode);
                    nunode->SetAnc (p->GetAnc());
					break;

				case SUPERSET:
					if (q == p->GetChild())
					{
                    	p->SetChild (nunode);
                        nunode->SetChild (q);
                        nunode->SetAnc (p);
                        q->SetAnc (nunode);
					}
					else
					{
	                   	p->SetSibling (nunode);
                        nunode->SetChild (q);
                        nunode->SetAnc (p->GetAnc());
                        q->SetAnc (nunode);
					}

                    NNodePtr tmp = q;
                    NNodePtr s = (NNodePtr)(q->GetSibling());
                    NNodePtr t = (NNodePtr)(q->GetAnc());

                    while (s)
                    {
                    	r = Relationship (s->Cluster, b);
                        if ((r == IDENTITY) || (r == SUBSET))
                        {
                        	s->SetAnc(nunode);
                            tmp = s;
                            s = (NNodePtr)(s->GetSibling());
                        }
                        else
                        {
                        	t->SetSibling (s);
                            tmp->SetSibling (s->GetSibling());
                            t = s;
                            t->SetSibling (NULL);
                            s = (NNodePtr)(tmp->GetSibling());
                       }
                    }
                    break;

//					checkSiblings (q, nunode->cluster);
				}
//			nunode->GetDegreeOf ();
//			nunode->anc->degree = nunode->degree - nunode->degree - 1;
		}
	}
	return Error;
}


