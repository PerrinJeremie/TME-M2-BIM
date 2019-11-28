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

// $Id: treeorder.cpp,v 1.2 2005/02/23 14:19:00 rdmp1c Exp $

#include "treeorder.h"

#include <iostream>

//------------------------------------------------------------------------------
void TreeOrder::Order ()
{
    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (!q->IsLeaf ())
			SortDescendants (q);

        q = n.next();
    }
}

//------------------------------------------------------------------------------
void TreeOrder::SortDescendants (NodePtr node)
{
	NodePtr head = node->GetChild ();
	NodePtr tail = head;
	while (tail->GetSibling () != NULL)
	{
		NodePtr p = tail->GetSibling ();
		if (MustSwap (head, p))
		{
			tail->SetSibling (p->GetSibling ());
			p->SetSibling (head);
			head = p;
			p->GetAnc()->SetChild (p);
		}
		else
		{
			NodePtr q = head;
			NodePtr r = q->GetSibling ();
			while (MustSwap (p, r))
			{
				q = r;
				r = q->GetSibling ();
			}
			if (p == r)
				tail = p;
			else
			{
				tail->SetSibling (p->GetSibling ());
				p->SetSibling (r);
				q->SetSibling (p);
			}
		}
	}
}

//------------------------------------------------------------------------------
bool LeftOrder::MustSwap (NodePtr p, NodePtr q)
{
	return (p->GetWeight() < q->GetWeight());
}


//------------------------------------------------------------------------------
bool RightOrder::MustSwap (NodePtr p, NodePtr q)
{
	return (p->GetWeight() > q->GetWeight());
}

//------------------------------------------------------------------------------
bool AlphaOrder::MustSwap (NodePtr p, NodePtr q)
{
	return (labels[p] > labels[q]);
}

//------------------------------------------------------------------------------
void AlphaOrder::Order ()
{
    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        	labels[q] = q->GetLabel();
        q = n.next();
    }
	q = n.begin();
    while (q)
    while (q)
    {
    	if (!q->IsLeaf ())
        {
			SortDescendants (q);
            labels[q] = labels[q->GetChild()];
        }

        q = n.next();
    }
}

//------------------------------------------------------------------------------
void ThorneOrder::Order ()
{
	if (outgroup->GetSibling() != NULL)
	{
		NodePtr q = t->GetRoot()->GetChild()->GetRightMostSibling();
		
		//std::cout << "OG" << outgroup->GetLabel() << std::endl;

		if (outgroup == t->GetRoot()->GetChild())
		{
			t->GetRoot()->SetChild(outgroup->GetSibling());
		}
		else
		{
			NodePtr p = outgroup->LeftSiblingOf();
			p->SetSibling (outgroup->GetSibling());
		}
		q->SetSibling (outgroup);
		outgroup->SetSibling (NULL);
	}
	
}

//------------------------------------------------------------------------------
bool ThorneOrder::MustSwap (NodePtr p, NodePtr q)
{
	return (p == outgroup);
}

