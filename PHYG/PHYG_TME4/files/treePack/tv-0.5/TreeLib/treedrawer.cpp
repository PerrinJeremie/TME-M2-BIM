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

// $Id: treedrawer.cpp,v 1.20 2007/06/07 15:48:02 rdmp1c Exp $

#include <math.h>
#include <stdio.h>

#include "treedrawer.h"

#ifndef M_PI
	#define M_PI		3.14159265358979323846	// define pi
#endif

#define DEGREES_TO_RADIANS(d) ((d / 180.0) * M_PI)
#define RADIANS_TO_DEGREES(r) ((180.0 * r) / M_PI)

//------------------------------------------------------------------------------
TreeDrawer::TreeDrawer (Tree *tree)
{
	t = tree;
    rooted = true;
    showInternalLabels = true;
    showLeafLabels = true;
    left = 0.0;
    top = 0.0;
    width = 400.0;
    height = 400.0;
    leafCount = 0;
    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        {
        	point p;
        	node_coordinates[q] = p;
        }
        q = n.next();
    }
}

//------------------------------------------------------------------------------
void TreeDrawer::CalcLeaf (Node *p)
{
	node_coordinates[p].y = top + (double)leafCount * leafGap;
    lastY = node_coordinates[p].y;
    leafCount++;

    // Cladogram
	node_coordinates[p].x = left + width;
}

//------------------------------------------------------------------------------
void TreeDrawer::CalcInternal (Node *p)
{
	// Cladogram
	node_coordinates[p].x = left + (double)nodeGap * (double)(t->GetNumLeaves() - p->GetWeight());

    // Slant
    node_coordinates[p].y = lastY - (((double)(p->GetWeight() - 1) * leafGap) / 2.0);
}


//------------------------------------------------------------------------------
void TreeDrawer::CalcCoordinates ()
{
	double l = t->GetNumLeaves();
    leafGap = height / (l - 1.0);
	if (rooted)
    	nodeGap = width / l;
    else
    	nodeGap = width / (l - 1.0);
	leafCount = 0;

    if (rooted)
    {
    	// Allow for edge below root
    	left += nodeGap;
        width -= nodeGap;
    }

    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        {
        	CalcLeaf (q);
        }
        else
        {
			CalcInternal (q);
        }

        q = n.next();
    }

}

//------------------------------------------------------------------------------
void TreeDrawer::DrawLeaf (Node *p)
{
	// cladogram (slant)
    NodePtr anc = p->GetAnc();
    if (anc)
    {
    	DrawLine (node_coordinates[p], node_coordinates[anc]);
    }
    DrawLeafLabel (p);
}

//------------------------------------------------------------------------------
void TreeDrawer::DrawInternal (Node *p)
{
	// cladogram (slant)
    NodePtr anc = p->GetAnc();
    if (anc)
    {
    	DrawLine (node_coordinates[p], node_coordinates[anc]);
    }
    DrawInternalLabel (p);
}


//------------------------------------------------------------------------------
void TreeDrawer::DrawLeafLabel (Node *p)
{
	if (showLeafLabels && (p->GetLabel() != ""))
    {
    	point pt = node_coordinates[p];

#if USE_PORT
		pt.x += font.GetSize()/2;
       	pt.y += font.GetSize()/3;
#endif
		DrawText (pt, p->GetLabel());
	}
}

//------------------------------------------------------------------------------
void TreeDrawer::DrawInternalLabel (Node *p)
{
	if (showInternalLabels && (p->GetLabel() != ""))
		DrawText (node_coordinates[p], p->GetLabel());
}



//------------------------------------------------------------------------------
void TreeDrawer::Draw ()
{
    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        {
        	DrawLeaf (q);
        }
        else
        {
			DrawInternal (q);
        }

        q = n.next();
    }
    if (rooted)
    {
    	DrawRoot ();
    }
}

//------------------------------------------------------------------------------
void TreeDrawer::DrawRoot ()
{
    	// Draw edge below root
    	point pt2 = node_coordinates[t->GetRoot()];
        point pt1 = pt2;
        pt1.x -= nodeGap;
        DrawLine (pt1, pt2);
}

//------------------------------------------------------------------------------
void TreeDrawer::DrawLine (point pt1, point pt2)
{
#if USE_WXWINDOWS
	dc->DrawLine ((int)pt1.x, (int)pt1.y, (int)pt2.x, (int)pt2.y);
#endif

#if USE_VC2
	Port.DrawLine ((int)pt1.x, (int)pt1.y, (int)pt2.x, (int)pt2.y);
#endif

#if USE_PORT
	Port->DrawLine ((int)pt1.x, (int)pt1.y, (int)pt2.x, (int)pt2.y);
#endif

}

//------------------------------------------------------------------------------
void TreeDrawer::DrawText (point pt, std::string s)
{
#if wxUSE_UNICODE
	std::wstring formatedString = L"";
#else
	std::string formatedString = "";
#endif
	
	int i = 0;
	while (i < s.size())
	{
		if (s[i] == '_')
			formatedString += ' ';
		else
			formatedString += s[i];
		i++;
	}

#if USE_VC2
	pt.x += Port.GetMaxCharWidth()/2;
	pt.y += Port.GetFontHeight()/2;
	Port.DrawText (pt.x, pt.y, (char *)formatedString.c_str());
#endif		

#if USE_WXWINDOWS
	// We add a new level of scope to avoid a "Declaration of 's' shadows a parameter"
	// error in gcc, which is probably a gcc bug
	{
		wxCoord w, h, descent;
		wxString s (formatedString.c_str(), wxSTRING_MAXLEN);
		pt.x += dc->GetCharWidth();
		pt.y -= dc->GetCharHeight()/2;
		dc->DrawText (s, (int)pt.x, (int)pt.y);
	}
#endif		


#if USE_PORT
	Port->DrawText (pt.x, pt.y, (char *)formatedString.c_str());
#endif
}


//------------------------------------------------------------------------------
void TreeDrawer::SetPenWidth (int width)
{
#if USE_VC2
	::PenSize (width, width);
#endif

#if USE_PORT
#endif
}


//------------------------------------------------------------------------------
void RectangleTreeDrawer::CalcCoordinates ()
{
	t->MakeNodeList();
    maxDepth = 0;
    // Clear internal node depths
    for (int i = t->GetNumLeaves(); i < t->GetNumNodes(); i++)
    {
    	(*t)[i]->SetDepth(0);
    }
    for (int i = 0; i < t->GetNumLeaves();  i++)
    {
    	NodePtr p = (*t)[i]->GetAnc();
        int count = 1;
        while (p)
        {
        	if (count > p->GetDepth())
            {
            	p->SetDepth(count);
                if (count > maxDepth)
                	maxDepth = count;
            }
            count++;
            p = p->GetAnc();
        }
    }

 	double l = t->GetNumLeaves();
    leafGap = height / (l - 1.0);
	l = maxDepth + 1.0;
	if (rooted)
    	nodeGap = width / l;
    else
    	nodeGap = width / (l - 1.0);
	leafCount = 0;

    if (rooted)
    {
    	// Allow for edge below root
    	left += nodeGap;
        width -= nodeGap;
    }

    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        {
        	CalcLeaf (q);
        }
        else
        {
			CalcInternal (q);
        }

        q = n.next();
    }

}


//------------------------------------------------------------------------------
void RectangleTreeDrawer::CalcInternal (Node *p)
{
	// Cladogram
	node_coordinates[p].x = left +  (double)nodeGap * (double)(maxDepth - p->GetDepth());

    // Rectangular style
    node_coordinates[p].y =
	    node_coordinates[p->GetChild()].y + (node_coordinates[p->GetChild()->GetRightMostSibling()].y - node_coordinates[p->GetChild()].y)/2.0;
}

//------------------------------------------------------------------------------
void RectangleTreeDrawer::DrawLeaf (Node *p)
{
    NodePtr anc = p->GetAnc();
    if (anc)
    {
    	point pt;
        pt.x = node_coordinates[anc].x;
        pt.y = node_coordinates[p].y;
    	DrawLine (node_coordinates[p], pt);
    }
    DrawLeafLabel (p);
}

//------------------------------------------------------------------------------
void RectangleTreeDrawer::DrawInternal (Node *p)
{
	// cladogram (slant)
    NodePtr anc = p->GetAnc();
    if (anc)
    {
   		point pt;
        pt.x = node_coordinates[anc].x;
        pt.y = node_coordinates[p].y;
    	DrawLine (node_coordinates[p], pt);
    }

    point pt1, pt2;
    pt1.x = node_coordinates[p].x;
    pt1.y = node_coordinates[p->GetChild()].y;
    pt2.x = node_coordinates[p].x;
    pt2.y = node_coordinates[p->GetChild()->GetRightMostSibling()].y;
    DrawLine (pt1, pt2);

    DrawInternalLabel (p);
}



//------------------------------------------------------------------------------
void PhylogramDrawer::CalcCoordinates ()
{
	// 1. Get path lengths
	mMaxPathLength = 0.0;
	t->GetRoot()->SetPathLength (t->GetRoot()->GetEdgeLength()); // modify for rooted?


    // duh! this needs to be preorder!!!!!!!!!

    PreorderIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
       	double d = q->GetEdgeLength();
        if (d < 0.00001)
        	d = 0.0;
        if (q != t->GetRoot())
	    	q->SetPathLength (q->GetAnc()->GetPathLength() + d);
		if (q->GetPathLength() > mMaxPathLength)
			mMaxPathLength = q->GetPathLength();
        q = n.next();
    }

    // Is the tree ultrametric? (should really be a method of the tree class...)
    mUltrametric = true;
    NodeIterator <Node> u (t->GetRoot());
    q = u.begin();
    while (q && mUltrametric)
    {
        if (q->IsLeaf())
        {
            double d = q->GetPathLength() - mMaxPathLength;
#if (__GNUC__ == 3 && __GNUC_MINOR__ > 3) || __GNUC__ > 3	    
		// fabs is no longer in the std namespace in gcc 3.4 or greater (sigh)
            mUltrametric = (fabs(d) <= 0.0001);
#else
            mUltrametric = (std::fabs(d) <= 0.0001);
#endif	    
//            cout << mMaxPathLength << ":" << q->GetPathLength() << " " << d << endl;
        }
        q = u.next();
    }
            

    // Allow for scale bar

#if USE_VC2
	scalebar_space = Port.GetFontHeight() * 2;
#endif		

#if USE_WXWINDOWS
	scalebar_space = dc->GetCharHeight() * 2;
#endif		


#if USE_PORT
	scalebar_space = font.GetSize() * 2;
#endif

	height -= scalebar_space;

 	double l = t->GetNumLeaves();
    leafGap = height / (l - 1.0);
	leafCount = 0;

    NodeIterator <Node> po (t->GetRoot());
    q = po.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        {
        	CalcLeaf (q);
        }
        else
        {
			CalcInternal (q);
        }

        q = po.next();
    }

}


//------------------------------------------------------------------------------
void PhylogramDrawer::CalcInternal (Node *p)
{
    node_coordinates[p].x = left +  (p->GetPathLength() / mMaxPathLength) * width;

    // Rectangular style
    node_coordinates[p].y =
    node_coordinates[p->GetChild()].y + (node_coordinates[p->GetChild()->GetRightMostSibling()].y - node_coordinates[p->GetChild()].y)/2.0;
}

//------------------------------------------------------------------------------
void PhylogramDrawer::CalcLeaf (Node *p)
{
    node_coordinates[p].x = left +  (p->GetPathLength() / mMaxPathLength) * width;
    node_coordinates[p].y = top + (double)leafCount * leafGap;
    lastY = node_coordinates[p].y;
    leafCount++;

}

//------------------------------------------------------------------------------
void PhylogramDrawer::Draw ()
{
    TreeDrawer::Draw ();
    DrawScaleBar ();
}

//------------------------------------------------------------------------------
void PhylogramDrawer::DrawScaleBar ()
{
    point pt1, pt2;
    // Draw scale bar using "log" scale
    float m = log10 (mMaxPathLength);
    int i = (int) m;
    if (!mUltrametric)
        i -= 1;
    float bar = pow (10.0, i);
    float scalebar = (bar / mMaxPathLength) * (float)width;
    //int scalebar = (int)((bar / mMaxPathLength) * width);


    // Pen to draw scale line
#if USE_WXWINDOWS
    wxPen scale_pen (*wxBLACK, 1, wxSOLID);
    dc->SetPen (scale_pen);
#endif

    
    if (mUltrametric)
    {
      // Draw scale bar that is the length of the tree
      pt1.x = left;
      pt1.y = top + height + scalebar_space;
      pt2.x = left + width;
      pt2.y = pt1.y;

      // Draw
      DrawLine (pt1, pt2);


      // Draw ticks and labels
     int num_ticks = (int)(mMaxPathLength/bar);
     int tick_height = scalebar_space/3;
     for (int k = 0; k <= num_ticks; k++)
     {
         pt1.x = left + width - scalebar * (float)k;
         pt1.y = top + height + scalebar_space;
         pt2.x = pt1.x;
         pt2.y = pt1.y - tick_height;
         DrawLine (pt1, pt2);

         // Label
#if wxUSE_UNICODE
        wchar_t buf[16];
        if (i >= 0)
        {
            swprintf (buf, wcslen(buf), L"%d", int (bar * k));
        }
        else
        {
            int j = abs (i);
            swprintf (buf, wcslen(buf), L"%.*f", j, bar * (float)k);
        }
#else
        char buf[16];
        if (i >= 0)
        {
            sprintf (buf, "%d", int (bar * k));
        }
        else
        {
            int j = abs (i);
            sprintf (buf, "%.*f", j, bar * (float)k);
        }
#endif

#if USE_WXWINDOWS
        wxCoord w, h;
        wxString s (buf, wxSTRING_MAXLEN);
        dc->GetTextExtent (s, &w, &h);
        int x = (int)pt2.x;
        int y = (int)pt2.y;
        x -= w/2;
        y -= h;
        dc->DrawText (s, x, y);
#else
        std::string s = buf;
        DrawText (pt2, buf);

#endif
     }
    }
    else
    {
        // Simple bar to indicate scale
        pt1.x = left;
        pt1.y = top + height + scalebar_space/2;
        pt2 = pt1;
        pt2.x += scalebar;

        // Draw
        DrawLine (pt1, pt2);

        // Text to the right a la PAUP*
        char buf[64];
        if (i >= 0)
        {
            sprintf (buf, "%d", int (bar));
        }
        else
        {
            int j = abs (i);
            sprintf (buf, "%.*f", j, bar);
        }
        std::string s = buf;
        DrawText (pt2, buf);
    }
}


//------------------------------------------------------------------------------
void PhylogramDrawer::DrawRoot ()
{
    // Draw edge below root
    point pt2 = node_coordinates[t->GetRoot()];
    point pt1 = pt2;
    pt1.x -= (t->GetRoot()->GetPathLength() / mMaxPathLength) * width;
    DrawLine (pt1, pt2);
}




//------------------------------------------------------------------------------
void CircleTreeDrawer::CalcLeaf (Node *p)
{
    node_angle[p] = leaf_angle * (double)leafCount;
    node_radius[p] = leaf_radius;
    leafCount++;

    node_coordinates[p].x = node_radius[p] * cos (node_angle[p]);
    node_coordinates[p].y = node_radius[p] * sin (node_angle[p]);
}

//------------------------------------------------------------------------------
void CircleTreeDrawer::CalcInternal (Node *p)
{
    double left_angle = node_angle[p->GetChild()];
    double right_angle = node_angle[p->GetChild()->GetRightMostSibling()];
    node_angle[p] = left_angle + (right_angle - left_angle)/2.0;
    node_radius[p] = nodeGap * (double)(maxDepth - p->GetDepth());

    node_coordinates[p].x = node_radius[p] * cos (node_angle[p]);
    node_coordinates[p].y = node_radius[p] * sin (node_angle[p]);

    NodePtr q = p->GetChild();
    while (q)
    {
    	point pt;
        node_backarc[q] = pt;
    	node_backarc[q].x = node_radius[p] * cos (node_angle[q]);
        node_backarc[q].y = node_radius[p] * sin (node_angle[q]);
        q = q->GetSibling();
    }
}


//------------------------------------------------------------------------------
void CircleTreeDrawer::CalcCoordinates ()
{
    t->MakeNodeList();
    maxDepth = 0;
    // Clear internal node depths
    for (int i = t->GetNumLeaves(); i < t->GetNumNodes(); i++)
    {
    	(*t)[i]->SetDepth(0);
    }
    for (int i = 0; i < t->GetNumLeaves();  i++)
    {
    	NodePtr p = (*t)[i]->GetAnc();
        int count = 1;
        while (p)
        {
        	if (count > p->GetDepth())
            {
            	p->SetDepth(count);
                if (count > maxDepth)
                	maxDepth = count;
            }
            count++;
            p = p->GetAnc();
        }
    }

    leaf_angle = 2 * M_PI / t->GetNumLeaves();
    left = top = 0.0;
    width = height = 400.0;
    leaf_radius = width / 2.0;
    leafCount = 0;
    nodeGap = leaf_radius / double(maxDepth);
    origin.x = 0.0;
    origin.y = 0.0;

    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        {
        	CalcLeaf (q);
        }
        else
        {
			CalcInternal (q);
        }

        q = n.next();
    }

    // Translate
    origin.x = left + (width/2.0);
    origin.y = top + (height/2.0);
    q = n.begin();
    while (q)
    {
    	node_coordinates[q].x += origin.x;
        node_coordinates[q].y += origin.y;
    	node_backarc[q].x += origin.x;
        node_backarc[q].y += origin.y;
        q = n.next();
    }


}

//------------------------------------------------------------------------------
void CircleTreeDrawer::DrawLeaf (Node *p)
{
	DrawLine (node_coordinates[p], node_backarc[p]);
    DrawLeafLabel (p);
}

//------------------------------------------------------------------------------
void CircleTreeDrawer::DrawInternal (Node *p)
{
    NodePtr anc = p->GetAnc();
    if (anc)
    {
    	DrawLine (node_coordinates[p], node_backarc[p]);
        double left_angle = node_angle[p->GetChild()];
    	double right_angle = node_angle[p->GetChild()->GetRightMostSibling()];
#if USE_PORT
        GPoint pt;
        pt.SetX((int)origin.x);
        pt.SetY((int)origin.y);
        Port->DrawArc (pt, node_radius[p], RADIANS_TO_DEGREES(left_angle), RADIANS_TO_DEGREES(right_angle));
#endif
    }
    DrawInternalLabel (p);
}



//------------------------------------------------------------------------------
void KMLTreeDrawer::CalcLeaf (Node *p)
{
	nodeLatitude[p] = p->GetLatitude();
	nodeLongitude[p] = p->GetLongitude();
	
	switch (treeStyle)
	{
		case 0:
		case 1:
			nodeAltitude[p] = p->GetWeight() * altitudeFactor;
			break;
			
		case 2:
			nodeAltitude[p] = (mMaxPathLength - p->GetPathLength()) * 500 * altitudeFactor;
			break;
	}
}

//------------------------------------------------------------------------------
void KMLTreeDrawer::CalcInternal (Node *p)
{
	if (treeStyle == 0)
	{
		// angle
		nodeLatitude[p] =  nodeLatitude[p->GetChild()] + (nodeLatitude[p->GetChild()->GetRightMostSibling()] -  nodeLatitude[p->GetChild()])/2.0;
		nodeLongitude[p] =  nodeLongitude[p->GetChild()] + (nodeLongitude[p->GetChild()->GetRightMostSibling()] -  nodeLongitude[p->GetChild()])/2.0;
	}
	else
	{
		// rectangular, so we need great circles
			double lat1 = DEGREES_TO_RADIANS(nodeLatitude[p->GetChild()]);
			double lon1 = DEGREES_TO_RADIANS(nodeLongitude[p->GetChild()]);
			double lat2 = DEGREES_TO_RADIANS(nodeLatitude[p->GetChild()->GetRightMostSibling()]);
			double lon2 = DEGREES_TO_RADIANS(nodeLongitude[p->GetChild()->GetRightMostSibling()]);
			
			// d=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2))
			double d = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));
		
			//cout << d << endl;
			
			if (d == 0.0)
			{
				// descendants are on same spot
				nodeLatitude[p] =  nodeLatitude[p->GetChild()] ;
				nodeLongitude[p] =  nodeLongitude[p->GetChild()];
			}
			else
			{ 
				double frac = 0.5;

				double A =sin((1-frac)*d)/sin(d);
				double B =sin(frac*d)/sin(d);
				double x = A*cos(lat1)*cos(lon1) +  B*cos(lat2)*cos(lon2);
				double y = A*cos(lat1)*sin(lon1) +  B*cos(lat2)*sin(lon2);
				double z = A*sin(lat1)           +  B*sin(lat2);
				double lat=atan2(z,sqrt(x*x+y*y));
				double lon=atan2(y,x);
				
				nodeLongitude[p] = RADIANS_TO_DEGREES(lon);
				nodeLatitude[p] = RADIANS_TO_DEGREES(lat);
			}
	
	
	}

	switch (treeStyle)
	{
		case 0:
			nodeAltitude[p] = p->GetWeight() * altitudeFactor;
			break;
		case 1:
			nodeAltitude[p] = p->GetWeight() * altitudeFactor;
			break;
			
		case 2:
			nodeAltitude[p] = (mMaxPathLength - p->GetPathLength()) * 500 * altitudeFactor;
			break;
	}
}


//------------------------------------------------------------------------------
void KMLTreeDrawer::CalcCoordinates ()
{
    t->MakeNodeList();
    maxDepth = 0;
    // Clear internal node depths
    for (int i = t->GetNumLeaves(); i < t->GetNumNodes(); i++)
    {
    	(*t)[i]->SetDepth(0);
    }
    for (int i = 0; i < t->GetNumLeaves();  i++)
    {
    	NodePtr p = (*t)[i]->GetAnc();
        int count = 1;
        while (p)
        {
        	if (count > p->GetDepth())
            {
            	p->SetDepth(count);
                if (count > maxDepth)
                	maxDepth = count;
            }
            count++;
            p = p->GetAnc();
        }
    }
	
	if (treeStyle == 2)
	{
		// 1. Get path lengths
		mMaxPathLength = 0.0;
		t->GetRoot()->SetPathLength (t->GetRoot()->GetEdgeLength()); // modify for rooted?


		// duh! this needs to be preorder!!!!!!!!!

		PreorderIterator <Node> pin (t->GetRoot());
		Node *q = pin.begin();
		while (q)
		{
			double d = q->GetEdgeLength();
			if (d < 0.00001)
				d = 0.0;
			if (q != t->GetRoot())
				q->SetPathLength (q->GetAnc()->GetPathLength() + d);
			if (q->GetPathLength() > mMaxPathLength)
				mMaxPathLength = q->GetPathLength();
			q = pin.next();
		}	
	}
    leafCount = 0;

    NodeIterator <Node> n (t->GetRoot());
    Node *q = n.begin();
    while (q)
    {
    	if (q->IsLeaf ())
        {
        	CalcLeaf (q);
        }
        else
        {
			CalcInternal (q);
        }

        q = n.next();
    }

}

//------------------------------------------------------------------------------
void KMLTreeDrawer::DrawLeaf (Node *p)
{
	DrawLeafLabel (p);
	
    NodePtr anc = p->GetAnc();
    if (anc)
    {
		*f << "<Placemark>" << endl;
		
		*f << "<styleUrl>#treeLine</styleUrl>" << endl;	
				

		
		*f << "<LineString>" << endl;
		*f << "<altitudeMode>absolute</altitudeMode>" << endl;
		*f << "<coordinates>" << endl;
		
		*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
		*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[anc] << endl;
		
		if (treeStyle == 0)
		{
			*f << nodeLongitude[anc] << "," << nodeLatitude[anc] << "," << nodeAltitude[anc] << endl;
		}
		else
		{
			*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[anc] << endl;
		}
		
		*f << "</coordinates>" << endl;
		*f << "</LineString>" << endl;
		*f << "</Placemark>" << endl;
		
	}
}

//------------------------------------------------------------------------------
void KMLTreeDrawer::DrawInternal (Node *p)
{
	DrawInternalLabel (p);
	
    NodePtr anc = p->GetAnc();	
	
	
		*f << "<Placemark>" << endl;
		
		*f << "<styleUrl>#treeLine</styleUrl>" << endl;	
		
		*f << "<LineString>" << endl;
		*f << "<altitudeMode>absolute</altitudeMode>" << endl;
		*f << "<coordinates>" << endl;
		


		if (treeStyle == 0)
		{
			*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
			if (anc)
			{
				*f << nodeLongitude[anc] << "," << nodeLatitude[anc] << "," << nodeAltitude[anc] << endl;
			}
		}
		else
		{
			// Code to handle curvature...
		
			double lat1 = DEGREES_TO_RADIANS(nodeLatitude[p->GetChild()]);
			double lon1 = DEGREES_TO_RADIANS(nodeLongitude[p->GetChild()]);
			double lat2 = DEGREES_TO_RADIANS(nodeLatitude[p->GetChild()->GetRightMostSibling()]);
			double lon2 = DEGREES_TO_RADIANS(nodeLongitude[p->GetChild()->GetRightMostSibling()]);
			
			// d=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2))
			double d = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));
		
			cout << d << endl;
			
			if (d == 0.0)
			{
				// descendants are on same spot, draw stem only
				if (anc)
				{
					*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[anc] << endl;
					*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
				}
			}
			else
			{ 
				// Draw arc and stem
				for (double frac =0.0; frac <= 1.0; frac += 0.1)
				{
					// A=sin((1-f)*d)/sin(d)
					double A =sin((1-frac)*d)/sin(d);
					double B =sin(frac*d)/sin(d);
					double x = A*cos(lat1)*cos(lon1) +  B*cos(lat2)*cos(lon2);
					double y = A*cos(lat1)*sin(lon1) +  B*cos(lat2)*sin(lon2);
					double z = A*sin(lat1)           +  B*sin(lat2);
					double lat=atan2(z,sqrt(x*x+y*y));
					double lon=atan2(y,x);
					
					*f << RADIANS_TO_DEGREES(lon) << "," << RADIANS_TO_DEGREES(lat) << "," << nodeAltitude[p] << endl;
					
					if (frac == 0.5)
					{
						if (anc)
						{
							*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[anc] << endl;
							*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
						}						
					}
				}
			}
			
		
			// Original straight code
			/*
			
		
			// left
			*f << nodeLongitude[p->GetChild()] << "," << nodeLatitude[p->GetChild()] << "," << nodeAltitude[p] << endl;
			// centre
			*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
			// stem
			if (anc)
			{
				*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[anc] << endl;
				*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
			}
			// right
			*f << nodeLongitude[p->GetChild()->GetRightMostSibling()] << "," << nodeLatitude[p->GetChild()->GetRightMostSibling()] << "," << nodeAltitude[p] << endl;
			*/
		}
		
		*f << "</coordinates>" << endl;
		*f << "</LineString>" << endl;
		
		
		*f << "</Placemark>" << endl;
}

//------------------------------------------------------------------------------
void KMLTreeDrawer::DrawRoot ()
{
}

//------------------------------------------------------------------------------
void KMLTreeDrawer::Draw ()
{
	cur = t->GetRoot();
	count = 0;

    while (cur)
    {
        if (cur->GetChild())
        {
             count++;
					
			*f << "<Folder>" << endl;
			
			if (cur == t->GetRoot())
			{
				*f << "<name>Tree";
				if (t->GetName() != "")
				{
					*f << ": " << t->GetName();
				}
				*f << "</name>" << endl;
			}
			else
			{
				if (cur->GetLabel() != "")
				{
					*f << "<name>";
					*f << cur->GetLabel();
					*f << "</name>" << endl;
				}
			}
				
			
			if (cur->IsLeaf()) DrawLeaf(cur);
			else DrawInternal(cur);
			
            stk.push (cur);
            cur = cur->GetChild();
        }
        else
        {
            count++;
			
			*f << "<Folder>" << endl;
			if (cur->IsLeaf()) DrawLeaf(cur);
			else DrawInternal(cur);
			*f << "</Folder>" << endl;
			
            while (!stk.empty() && (cur->GetSibling() == NULL))
            {
			*f << "</Folder>" << endl;
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
	
	//
	// Taxon labels
	*f << "<Folder>" << endl << "<name>Taxon Labels</name>" << endl;
    for (int i = 0; i < t->GetNumLeaves();  i++)
    {
		NodePtr p = (*t)[i];
		*f << "<Placemark>" << endl;
		*f << "<name>" << p->GetLabel() << "</name>" << endl;
		*f << "<styleUrl>#whiteBall</styleUrl>" << endl;
		*f << "<Point>" << endl;
		*f << "<altitudeMode>absolute</altitudeMode>" << endl;
		*f << "<extrude>1</extrude>" << endl;
		*f << "<coordinates>" << endl;
		*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
		*f << "</coordinates>" << endl;
		*f << "</Point>" << endl;
		*f << "</Placemark>" << endl;
    }
	*f << "</Folder>" << endl;


}


//------------------------------------------------------------------------------
void KMLTreeDrawer::DrawLeafLabel (Node *p)
{
	if (showLeafLabels && (p->GetLabel() != ""))
    {
    	*f << "<name>" << p->GetLabel() << "</name>" << endl;
	}
}

//------------------------------------------------------------------------------
void KMLTreeDrawer::DrawInternalLabel (Node *p)
{
	if (showInternalLabels && (p->GetLabel() != ""))
    {
		*f << "<Placemark>" << endl;
		*f << "<name>" << p->GetLabel() << "</name>" << endl;
		*f << "<styleUrl>#whiteBall</styleUrl>" << endl;
		*f << "<Point>" << endl;
		*f << "<altitudeMode>absolute</altitudeMode>" << endl;
		*f << "<coordinates>" << endl;
		*f << nodeLongitude[p] << "," << nodeLatitude[p] << "," << nodeAltitude[p] << endl;
		*f << "</coordinates>" << endl;
		*f << "</Point>" << endl;
		*f << "</Placemark>" << endl;
	}
}



