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
 
 // $Id: treewriter.cpp,v 1.4 2005/02/23 14:32:44 rdmp1c Exp $

#include "treewriter.h"

//------------------------------------------------------------------------------
void NewickTreeWriter::WriteLeftParenthesis ()
{
	*f << '(';
}

//------------------------------------------------------------------------------
void NewickTreeWriter::WriteRightParenthesis ()
{
	*f << ')';
}

//------------------------------------------------------------------------------
void NewickTreeWriter::WriteSiblingSymbol ()
{
	*f << ',';
}

//------------------------------------------------------------------------------
void NewickTreeWriter::WriteLeaf ()
{
	*f << NEXUSString (cur->GetLabel());
	if (t->GetHasEdgeLengths () && writeEdgeLengths)
		*f << ":" << cur->GetEdgeLength();
}


//------------------------------------------------------------------------------
void NewickTreeWriter::WriteInternal ()
{
	if (cur->GetLabel() != "")
		*f << NEXUSString (cur->GetLabel());
	if (t->GetHasEdgeLengths () && writeEdgeLengths)
		*f << ":" << cur->GetEdgeLength();
}

//------------------------------------------------------------------------------
void NewickTreeWriter::WriteEndOfTree ()
{
	*f << ';';
}


//------------------------------------------------------------------------------
void NewickTreeWriter::Write ()
{
	cur = t->GetRoot();

    while (cur)
    {
        if (cur->GetChild())
        {
            WriteLeftParenthesis ();
            stk.push (cur);
            cur = cur->GetChild();
        }
        else
        {
            WriteLeaf ();
            while (!stk.empty() && (cur->GetSibling() == NULL))
            {
                WriteRightParenthesis ();
                cur = stk.top();
                WriteInternal ();
                stk.pop();
            }
            if (stk.empty())
                cur = NULL;
            else
            {
                WriteSiblingSymbol ();
                cur = cur->GetSibling();
            }
        }
    }
    WriteEndOfTree ();
}


//------------------------------------------------------------------------------
void XMLTreeWriter::Write ()
{
    OpenTree();

	cur = t->GetRoot();
        count = 0;

    while (cur)
    {
        if (cur->GetChild())
        {
             count++;
             OpenNode();
             NodeInfo();
                
            stk.push (cur);
            cur = cur->GetChild();
        }
        else
        {
            count++;
            OpenNode();
            NodeInfo();
            CloseNode();
            while (!stk.empty() && (cur->GetSibling() == NULL))
            {
                CloseNode();
                cur = stk.top();
                stk.pop();
            }
            if (stk.empty())
                cur = NULL;
            else
            {
                cur = cur->GetSibling();
            }
        }
    }
    CloseTree();

}

void TreeBolicTreeWriter::OpenTree ()
{
    *f << "<tree>" << std::endl;
    *f << "<nodes>" << std::endl;
}

void TreeBolicTreeWriter::CloseTree ()
{
    *f << "</tree>" << std::endl;
    *f << "</nodes>" << std::endl;
}

void TreeBolicTreeWriter::OpenNode ()
{
    *f << std::string (stk.size() * 3, ' ') << "<node id = \"" << count << "\">" << std::endl;
}

void TreeBolicTreeWriter::CloseNode ()
{
    *f << std::string (stk.size() * 3, ' ') << "</node>" << std::endl;
}

void TreeBolicTreeWriter::NodeInfo ()
{
    // label
    if (cur->GetLabel() != "")
        *f << std::string (stk.size() * 3 + 3, ' ') << "<label>" << cur->GetLabel() << "</label>" << std::endl;
}

//------------------------------------------------------------------------------
void ThorneTreeWriter::WriteLeaf ()
{
	*f << formatLabel (cur->GetLabel());
}


//------------------------------------------------------------------------------
void ThorneTreeWriter::WriteInternal ()
{
}

std::string ThorneTreeWriter::formatLabel(std::string s)
{	
	std::string  sub = s.substr(0,30);
	std::string thorne_label = "";
	int k = 0;
	while (k < sub.length())
	{
		switch (sub[k])
		{
			case '\'': thorne_label += '_'; break;
			case ' ': thorne_label += '_'; break;
			default: thorne_label += sub[k]; break;
		}
		k++;
	}
	return thorne_label;
}


