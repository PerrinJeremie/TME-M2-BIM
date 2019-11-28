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
 
 // $Id: profile.h,v 1.31 2007/06/05 13:26:17 rdmp1c Exp $
 
/**
 * @file profile.h
 *
 * Storage of trees
 *
 */
#ifndef PROFILE_H
#define PROFILE_H

#ifdef __BORLANDC__
	// Undefine __MINMAX_DEFINED so that min and max are correctly defined
	#ifdef __MINMAX_DEFINED
		#undef __MINMAX_DEFINED
	#endif
    // Ignore "Cannot create precompiled header: code in header" message
    // generated when compiling string.cc
    #pragma warn -pch
#endif

#include "TreeLib.h"

#include "treereader.h"
#include "treewriter.h"

// NCL includes
#include "nexusdefs.h"
#include "xnexus.h"
#include "nexustoken.h"
#include "nexus.h"
#include "taxablock.h"
#include "assumptionsblock.h"
#include "treesblock.h"
#include "discretedatum.h"
#include "discretematrix.h"
#include "charactersblock.h"
#include "datablock.h"

#if USE_XML
	#include "xml.h"
#endif

#if USE_VC2
	#include "VMsg.h"
#endif

#if USE_WXWINDOWS
   	#include "wx/wx.h"
#endif

#if (__BORLANDC__ < 0x0550)
	#include <time.h>
#else
	#include <ctime>
#endif    



/**
 *@typedef std::map <std::string, int, std::less<std::string> > LabelMap;
 */
typedef std::map <std::string, int, std::less<std::string> > LabelMap;


/**
 * @class Profile
 * Encapsulates reading and storing a set of trees.
 *
 */
template <class T> class Profile
{
public:
	/**
	 * A map between leaf labels in the profile and a unique integer index
	 */	
	LabelMap Labels;
	/**
	 * For each leaf label the number of trees in the profile that have the corresponding leaf
	 */	
	LabelMap LabelFreq;
	
	/**
	 * Constructor
	 */
	Profile () {};
	/**
	 * Destructor
	 */
	virtual ~Profile () {};

	/**
	 * @brief The ith tree in the profile
	 *
	 * @param i the index of the tree in the range 0 - (n-1)
	 * @return The tree
	 */
	 virtual T GetIthTree (int i) { return Trees[i]; };
	/**
	 * @brief The name of the ith tree in the profile
	 *
	 * @param i the index of the tree in the range 0 - (n-1)
	 * @return The tree
	 */
	 virtual std::string GetIthTreeName (int i) { return Trees[i].GetName(); };
	/**
	 * @return The number of labels in the profile
	 */
	virtual int GetNumLabels () { return Labels.size(); };
	/**
	 * @return The number of trees in the profile
	 */
	virtual int GetNumTrees () { return Trees.size(); };
	/**
	 * @brief The index of a leaf label
	 * @param s A leaf label
	 * @return The index of the leaf label
	 */
	virtual int GetIndexOfLabel (std::string s);
	/**
	 * @brief The unique index of a leaf label
	 * @param i A leaf index
	 * @return The ith leaf label
	 * @sa Profile::GetIndexOfLabel
	 */
	virtual std::string GetLabelFromIndex (int i) { return LabelIndex[i]; };

	/**
	 * @brief Assign a unique integer index to each leaf label in the profile
	 */
	virtual void MakeLabelList ();
	/**
	 * @brief Count the number of trees each leaf label occurs in
	 */
	virtual void MakeLabelFreqList ();
	/**
	 * @brief Read a NEXUS file and store the trees in Profile::trees. 
	 
	 * If a TAXA block is
	 * present then the leaf labels are stored in Labels in the same order as in
	 * the TAXA block, otherwise Profile::MakeLabelList is called to assign a unique integer
	 * index to each label.
	 *
	 * @param f input stream in NEXUS format
	 * @return true if sucessful
	 */
	virtual bool ReadNEXUS (std::istream &f);
	/**
	 * @brief Read a PHYLIP tree file and store the trees in Profile::trees. 
	 *
	 * @param f input stream in PHYLIP format
	 * @return true if sucessful
	 */
 	virtual bool ReadPHYLIP (std::istream &f);
	/**
	 * @brief Read a set of trees from an input stream
	 * At present PHYLIP, NEXUS, and a subset of PhyloXML
	 * formats are supported.
	 * @param f input stream 
	 * @return true if successful
	 */
	virtual bool ReadTrees (std::istream &f);
#if USE_XML
	/**
	 * @brief Read a XML file and store the trees in Profile::trees. 
	 * XML format is based on "XML for phylogenetic and population genetic data"
	 * (see http://evolve.zoo.ox.ac.uk/PhyloXML/).
	 *
	 * @param f input stream in PHYLIP format
	 * @return true if sucessful
	 */
 	virtual bool ReadXML (std::istream &f);
#endif
	/**
	 * @brief Output leaf labels.
	 *
	 * @param f output stream
	 */
	virtual void ShowLabelList (std::ostream &f);

	/**
	 * @brief Output trees as dendrograms.
	 *
	 * @param f output stream
	 */
    	virtual void ShowTrees (std::ostream &f);
	/**
	 * @brief Write a set of trees to an output stream
	 * @param f output stream 
	 * @param format file format to use (at present nexus only)
	 * @return true if successful
	 */
	virtual bool WriteTrees (std::ostream &f, const int format = 0, const char *endOfLine = "\n");

protected:
	/**
	 * The trees
	 *
	 */
	std::vector <T> Trees;
	/**
	 * The leaf labels stored as a vector so they can be accessed by an index value
	 *
	 */
	std::vector <std::string> LabelIndex;
	
#if USE_XML
	virtual bool xmlTraverse (XMLElementPtr p);
	virtual bool xmlCleanup (XMLElementPtr p);
	virtual void postProcessXMLTree (int tree_num);
#endif
	
	int curTreeNumber;
	IntegerNodeMap node_map;
};

//------------------------------------------------------------------------------
template <class T> int Profile<T>::GetIndexOfLabel (std::string s)
{
	return Labels[s];
}

//------------------------------------------------------------------------------
template <class T> void Profile<T>::MakeLabelList ()
{
	for (int i = 0; i < Trees.size(); i++)
	{
    	T t = Trees[i];
		t.MakeNodeList();
		for (int j = 0; j < t.GetNumLeaves(); j++)
		{
			std::string s = t[j]->GetLabel();
			
			if (Labels.find (s) == Labels.end ())
			{
				int index = Labels.size();
//				cout << "Labels.size() = " << Labels.size();
				Labels[s] = index;
				LabelIndex.push_back (s);
//				cout << "Labels[s] =" << Labels[s] << endl;
			}
		} 
	}
/*	cout << endl << "Number of labels = " << Labels.size() << endl;
	LabelMap::iterator it = Labels.begin();
	LabelMap::iterator end = Labels.end();
	while (it != end)
	{
		cout << (*it).first << " - " << (*it).second << endl;
		it++;
	}*/
}


//------------------------------------------------------------------------------
template <class T> void Profile<T>::ShowLabelList (std::ostream &f)
{
	f << std::endl << "Number of labels = " << Labels.size() << std::endl;
	LabelMap::iterator it = Labels.begin();
	LabelMap::iterator end = Labels.end();
	while (it != end)
	{
		f << (*it).first << " - " << (*it).second << std::endl;
		it++;
	}
	f << std::endl;
}

//------------------------------------------------------------------------------
template <class T> void Profile<T>::MakeLabelFreqList ()
{
//	cout << "MakeLabelListFreq" << endl;
	for (int i = 0; i < Trees.size(); i++)
	{
		T t = Trees[i];
		t.MakeNodeList();
		for (int j = 0; j < t.GetNumLeaves(); j++)
		{
			std::string s = t[j]->GetLabel();
			LabelMap::iterator there = LabelFreq.find (s);
			
			if (there == Labels.end ())
			{
				LabelFreq[s] = 1;
			}
			else
			{
				LabelFreq[s] += 1;
			}
		} 
	}
/*	cout << "Frequency of labels" << endl;
	LabelMap::iterator it = LabelFreq.begin();
	LabelMap::iterator end = LabelFreq.end();
	while (it != end)
	{
		cout << (*it).first << " - " << (*it).second << endl;
		it++;
	}*/
}

/**
 * @class MyNexus
 * Extends Nexus class to output progress to cout
 *
 */
class MyNexus : public Nexus
{
public:
	MyNexus () : Nexus() { isOK = true; };
	
#if (USE_VC2 || USE_WXWINDOWS)

	#if USE_VC2
		virtual void EnteringBlock( nxsstring blockName ) { }
		virtual void ExitingBlock( nxsstring blockName ) { };
		virtual void SkippingBlock( nxsstring blockName ) {  };
		virtual void SkippingDisabledBlock( nxsstring blockName ) { };
		virtual void ExecuteStarting() { };
		virtual void ExecuteStopping() { };
		virtual void OutputComment( nxsstring comment ) {};
	 	virtual void NexusError( nxsstring& msg, streampos pos, long line, long col )
		{
			char buf[256];
			sprintf (buf,"%s at line %d, column %d", msg.c_str(), line, col);
			Message (MSG_ERROR, "Error reading NEXUS file", buf);
			isOK = false;
		};
	#endif

	
	#if USE_WXWINDOWS
	#if 1
	   	 virtual void EnteringBlock( nxsstring blockName ) { }
	   	 virtual void ExitingBlock( nxsstring blockName ) { };
	   	 virtual void SkippingBlock( nxsstring blockName ) {  };
	   	 virtual void SkippingDisabledBlock( nxsstring blockName ) { };
		virtual void ExecuteStarting() { };
		virtual void ExecuteStopping() { };
	#else // debugging NEXUS reader
	   	 virtual void EnteringBlock( nxsstring blockName ) { wxLogMessage ("Entering %s block", blockName.c_str()); }
	   	 virtual void ExitingBlock( nxsstring blockName ) { wxLogMessage ("Exiting %s block", blockName.c_str()); };
	   	 virtual void SkippingBlock( nxsstring blockName ) { wxLogWarning ("Skipping %s block", blockName.c_str());  };
	   	 virtual void SkippingDisabledBlock( nxsstring blockName ) {  wxLogWarning ("Skipping disabled %s block", blockName.c_str());  };
		virtual void ExecuteStarting() { wxLogMessage ("Starting to execute NEXUS file"); };
		virtual void ExecuteStopping() { wxLogMessage ("Finished executing NEXUS file"); };
	#endif
		virtual void OutputComment( nxsstring comment ) { std::cout << comment << std::endl;};
	 	virtual void NexusError( nxsstring& msg, std::streampos pos, long line, long col )
		{
			wxLogError (wxT("%s at line %d, column %d"), msg.c_str(), line, col);
			isOK = false;
		};
	#endif

#else
	virtual void EnteringBlock( nxsstring blockName ) { cout << "   Entering " << blockName << " block..."; };
	virtual void ExitingBlock( nxsstring blockName ) { cout << "done" << endl; };
	virtual void SkippingBlock( nxsstring blockName ) { cout << "   (Skipping " << blockName << " block)" << endl; };
	virtual void SkippingDisabledBlock( nxsstring blockName ) { cout << "   (Skipping disabled " << blockName << " block)" << endl; };
	virtual void ExecuteStarting() { cout << "Starting to execute NEXUS file" << endl; };
	virtual void ExecuteStopping() { cout << "Finished executing NEXUS file" << endl; };
	virtual void OutputComment( nxsstring comment ) { cout << comment << endl;};
	virtual void NexusError( nxsstring& msg, streampos pos, long line, long col )
	{
   		cerr << "Error: " << msg << " line " << line << ", col " << col << endl;
		isOK = false;
	};
#endif
	bool GetIsOK () { return isOK; };

protected:
	bool isOK;
};

//------------------------------------------------------------------------------
template <class T> bool Profile<T>::ReadNEXUS (std::istream &f)
{
	bool result = false;

	TaxaBlock* taxa;
	DataBlock *data;
	CharactersBlock *characters;
	AssumptionsBlock *assumptions;
	TreesBlock *trees;

	taxa = new TaxaBlock();
	assumptions = new AssumptionsBlock (*taxa);
	data = new DataBlock (*taxa, *assumptions);
	characters = new CharactersBlock (*taxa, *assumptions);
	trees = new TreesBlock (*taxa);

	MyNexus nexus;
	nexus.Add( taxa );
	nexus.Add( data );
	nexus.Add( characters );
	nexus.Add( trees );
    


    // Set to binary to handle Mac and Unix files
#ifdef __MWERKS__
#elif __BORLANDC__
	f.setf (ios::binary);
#elif __GNUC__
	#if __GNUC__ < 3
		f.setf (ios::binary);
	#endif
#endif

	NexusToken token (f);
    

	try 
	{
    		nexus.Execute (token);
	}
	catch (XNexus x)
	{
		std::cout << x.msg << " (line " << x.line << ", column " << x.col << ")" << std::endl;
	}    	

	if (nexus.GetIsOK() && (trees->GetNumTrees() > 0))
	{

		// Display information about the trees
#if (USE_WXWINDOWS || USE_VC2)
#else
		trees->Report (std::cout);
		std::cout << endl;
#endif
		// Store the trees themselves
		int i = 0;
		int error = 0;
		while ((i < trees->GetNumTrees()) && (error == 0))
		{ 
			T t;
			std::string tstr;
//			if (trees->HasTranslationTable())
//				tstr = trees->GetTranslatedTreeDescription (i);
//			else
				tstr = trees->GetTreeDescription (i);
			tstr += ";";
			error = t.Parse (tstr.c_str());
			if (error == 0)
			{
				t.SetName (trees->GetTreeName (i));
				t.SetRooted (trees->IsRootedTree (i));
				t.SetWeight (trees->GetTreeWeight (i));
				
				if (trees->HasTranslationTable())
				{
					t.MakeNodeList();
					for (int k = 0; k < t.GetNumLeaves (); k++)
					{
						std::string skey = t[k]->GetLabel();
						if (skey != "")
						{
							nxsstring svalue = trees->GetTranslatedLabel(skey);
							if (svalue != "")
								t[k]->SetLabel (svalue);
						}
					}
				}
				
				
				
				Trees.push_back (t);
			}
			else
			{
#if USE_WXWINDOWS
				wxLogError ( wxT("Error in description of tree %d: %s"), (i+1), t.GetErrorMsg().c_str());
#elif USE_VC2
				char buf[256];
				sprintf (buf, "Reading tree %d: %s", (i+1), t.GetErrorMsg().c_str());
				Message (MSG_ERROR, "Error in tree description", buf);
#else
				cerr << "Error in tree description " << (i + 1) << t.GetErrorMsg() << endl;
#endif
				return false;
			}           
			 i++;
		}
        
		// Assign each label a unique index
		if (taxa->GetNumTaxonLabels() == 0)
		{
			// No taxa block in NEXUS file
			MakeLabelList ();
		}
		else
		{
			// Store the labels in the same order encountered in the
			// NEXUS file
			for (int i = 0; i < taxa->GetNumTaxonLabels (); i++)
			{
				Labels[taxa->GetTaxonLabel (i)] = i;
				LabelIndex.push_back (taxa->GetTaxonLabel (i));
			}
		}
		result = true;
	}
	return result;
}


//------------------------------------------------------------------------------
template <class T> bool Profile<T>::ReadTrees (std::istream &f)
{
	bool result = false;

	char ch = (char)f.peek ();
	if (ch == '#')
		result = ReadNEXUS (f);
	else if (strchr ("([", ch))
		result = ReadPHYLIP (f);
#if USE_XML
	else if (ch == '<')
		result = ReadXML (f);
#endif
	return result;
}

//------------------------------------------------------------------------------
template <class T> bool Profile<T>::ReadPHYLIP (std::istream &f)
{
	Tokeniser p (f);
	PHYLIPReader tr (p);
	bool ok = true;
	while (ok)
	{
		T t;

		try
		{
			ok = tr.Read (&t);
		}
		catch (XTokeniser x)
		{
#if USE_WXWINDOWS 
			wxLogError (wxT("%s at line %d, column %d"), x.msg.c_str(), x.line, x.col);           
#elif USE_VC2
			char buf[256];
			sprintf (buf, "%s at line %d, column %d", x.msg.c_str(), x.line, x.col);
			Message (MSG_ERROR, "Error reading tree file", buf);
#else
			std::cerr << x.msg << " (line " << x.line << ", column " << x.col << ")" << std::endl;
#endif
		 	return false;
		}

		if (ok)
			Trees.push_back (t);
	}
	
	bool result = (Trees.size() > 0);
	
	if (result)
	{
		// Build a list of labels in the profile, such that each label 
		// is assigned a unique index    	
		MakeLabelList ();
	}
	return result;
}

//------------------------------------------------------------------------------
template <class T> void Profile<T>::ShowTrees (std::ostream &f)
{
	for (int i = 0; i < Trees.size(); i++)
	{
		T t = Trees[i];
		t.Update ();
		t.Draw (f);
	}
}

//------------------------------------------------------------------------------
template <class T> bool Profile<T>::WriteTrees (std::ostream &f, const int format, const char *endOfLine)
{
	bool result = true;
	
	// Simple nexus tree file
	f << "#nexus" << endOfLine;
	f << endOfLine;
	f << "begin trees;";
		
	// Date the file
	f << " [Treefile written ";
	time_t timer = time(NULL);
	struct tm* tblock = localtime(&timer);
	char time_buf[64];
	strncpy (time_buf, asctime(tblock), sizeof (time_buf));
	char *q = strrchr (time_buf, '\n');
	if (q)
		*q = '\0';
	f << time_buf << "]" << endOfLine;

	for (int i = 0; i < Trees.size(); i++)
	{
		T t = Trees[i];
		f << "\ttree ";
		if (t.GetName() != "")
			f << NEXUSString (t.GetName());
		else
			f << "tree_" << (i+1);
		f << " = ";
		if (t.IsRooted())
			f << "[&R] ";
		else
			f << "[&U] ";
			
		// Tree
		NewickTreeWriter tw (&t);
		tw.SetStream (&f);
		tw.Write();	
		f << endOfLine;		
			
//		f << t << endOfLine;
	}
	f << "end;" << endOfLine;
	
	return result;
}


#if USE_XML
//------------------------------------------------------------------------------
template <class T> void Profile<T>::postProcessXMLTree (int tree_num)
{
	IntegerNodeMap::iterator it = node_map.begin();
	IntegerNodeMap::iterator end = node_map.end();
	
	if (it != end)
	{
	    // Dump tree
	    int nodeCount = 0;
	    int leafCount = 0;
		NodePtr root = NULL;
		while (it != end)
		{
			nodeCount++;
			int index = it->first;
			NodePtr p = it->second;
			
			if (p->GetAnc() == NULL)
				Trees[tree_num].SetRoot (p);

			if (p->GetChild() == NULL)
			{
				p->SetLeaf (true);
				p->SetLeafNumber (++leafCount);
			}
			else
			{
				if (p->GetLabel() != "")
					Trees[tree_num].SetInternalLabels (true);
			}
			
			it++;
		}
		Trees[tree_num].SetNumLeaves (leafCount);
		Trees[tree_num].SetNumInternals (nodeCount - leafCount);
	}
}



//------------------------------------------------------------------------------
template <class T> bool Profile<T>::xmlTraverse (XMLElementPtr p)
{
	if (p)
    {
        string attrname, attrvalue;
        if (p->name == "tree")
        {
          	attrname = "rooted";
            std::string rooted = p->attributes[attrname];

            // Finish any previous tree. If we are processing a new <tree> tag,
            // then now is our chance to finish the previous tre
            if (curTreeNumber > 0)
            {
            	postProcessXMLTree (curTreeNumber - 1);
            }
            
            // Clean up node list
            node_map.erase(node_map.begin(), node_map.end());
            
            // Create new tree
            T t;
            t.SetRooted (true);
			Trees.push_back(t); 
	        Trees[curTreeNumber].SetRooted (rooted == "true");
			
			curTreeNumber++;
	    }
        
        if (p->name == "node")
        {
          	attrname = "name";
            std::string node_name = p->attributes[attrname];
          	attrname = "id";
            std::string node_id = p->attributes[attrname];

            // Create node, set its attributes, use a map to access it by its id
            int i = atoi (node_id.c_str());
            Node *p = Trees[curTreeNumber-1].NewNode();
            node_map[i] =  p;
            if (node_name != "")
            	p->SetLabel (node_name);
            p->SetIndex (i);
	    }
        if (p->name == "edge")
        {
          	attrname = "source";
            attrvalue = p->attributes[attrname];
            int source = atoi (attrvalue.c_str());
            
          	attrname = "target";
            attrvalue = p->attributes[attrname];
            int target = atoi (attrvalue.c_str());

            // join ancestor-descendant pair, set any attributes of this edge
			Node *s = node_map[target];
			Node *q = node_map[source];
			s->SetAnc (q);
			if (q->GetChild())
			{
				Node *r = q->GetChild();
				r = r->GetRightMostSibling();
				r->SetSibling (s);
			}
			else 
				q->SetChild (s);   
			
			
          	attrname = "length";
            attrvalue = p->attributes[attrname];
            if (attrvalue != "")
            {
	            float length = atof (attrvalue.c_str());
	            s->SetEdgeLength (length);
	            Trees[curTreeNumber-1].SetEdgeLengths (true);
			}
	    }
	    
        if (p->name == "newick")
        {
            // parse Newick description 
    		Trees[curTreeNumber - 1].Parse (p->data.c_str());
	    }

		xmlTraverse (p->child);
        xmlTraverse (p->sib);
    }
}


//------------------------------------------------------------------------------
template <class T> bool Profile<T>::xmlCleanup (XMLElementPtr p)
{
	if (p)
    {
		xmlCleanup (p->child);
        xmlCleanup (p->sib);
        delete p;
    }
}


//------------------------------------------------------------------------------
template <class T> bool Profile<T>::ReadXML (istream &f)
{
	bool result = true;
  
    XMLElementPtr Root;
	myParser parser;

	char buf[BUFSIZ];
  
	while (f.good())
	{
		f.getline (buf, sizeof (buf));
		if (!parser.XML_Parse(buf, strlen (buf), false)) 
		{
#if USE_WXWINDOWS 
			wxLogError ("%s at line %d", XML_ErrorString(parser.XML_GetErrorCode()), parser.XML_GetCurrentLineNumber());           
#elif USE_VC2
			char buf[256];
			sprintf (buf, "%s at line %d", XML_ErrorString(parser.XML_GetErrorCode()), parser.XML_GetCurrentLineNumber());
			Message (MSG_ERROR, "Error reading tree file", buf);
#else
			cerr << XML_ErrorString(parser.XML_GetErrorCode()) << " at line " << parser.XML_GetCurrentLineNumber() << endl;
#endif
			return false;
		}
	}  
	
	// Parsed XML file OK, so now extract trees	and clean up
	curTreeNumber = 0;
	xmlTraverse (parser.Root);	
	postProcessXMLTree (curTreeNumber - 1);	
	xmlCleanup (parser.Root);

	return result;
}
#endif // USE_XML

#if __BORLANDC__
	// Redefine __MINMAX_DEFINED so Windows header files compile
	#ifndef __MINMAX_DEFINED
    		#define __MINMAX_DEFINED
	#endif
#endif

#endif

