#include "gport.h"


// System specific defines

#if GPORT_MAC
	// From MacTech 4(6) "Comments about PICTs"
	#define picGrpBeg	140
	#define picGrpEnd	141
#endif


// The global port
GBasePort *Port = NULL;


// A sensible default font
GBaseFont::GBaseFont ()
{
    description = "Times-Roman";
    name = "Times-Roman";
    size = 10;
    bold = false;
    italic = false;
}

GPostscriptPort::GPostscriptPort ()
{
    PenWidth = 1;
    DocumentFonts = "";
    Device = devPostscript;
    DisplayRect.SetRect (0, 0, 595-144, 842-144);
}

void GPostscriptPort::DrawArc (const GPoint &pt, const int radius,
	const double startAngleDegrees, const double endAngleDegrees)
{
	portStream << "newpath" << endl;
	portStream << pt.GetX() << " " << -pt.GetY()
    	<< " " << radius
        << " " << (360.0 -startAngleDegrees)
        << " " << (360.0 - endAngleDegrees)
        << " arcn" 		<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
}


void GPostscriptPort::DrawLine (const int x1, const int y1, const int x2, const int y2)
{
	portStream  << x2 << " " << -y2 << " " << x1 << " " << -y1 << " " << PenWidth << " DrawLine" << endl;
}

void GPostscriptPort::DrawCircle (const GPoint &pt, const int radius)
{
	portStream << "newpath" << endl;
	portStream << pt.GetX() << " " << -pt.GetY() << " " << radius << " 0 360 arc" 		<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
}


void GPostscriptPort::DrawRect (const GRect &r)
{
	portStream << r.GetLeft() << " " << -r.GetTop() << " moveto" 	<< endl;
	portStream << r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << -r.GetHeight() << " rlineto" 			<< endl;
	portStream << -r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << r.GetHeight() << " rlineto"				<< endl;
	portStream << "closepath" 									<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
}

void GPostscriptPort::DrawFilledRect (const GRect &r, double graylevel)
{
	portStream << r.GetLeft() << " " << -r.GetTop() << " moveto" 	<< endl;
	portStream << r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << -r.GetHeight() << " rlineto" 			<< endl;
	portStream << -r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << r.GetHeight() << " rlineto"				<< endl;
	portStream << "closepath" 									<< endl;
	portStream << "gsave"                                         << endl;
	portStream << graylevel << " setgray fill"                    << endl;
	portStream << "grestore"                                      << endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
}

void GPostscriptPort::DrawText (const int x, const int y, const char *text)
{
	portStream  << "(" << text << ") " << x << " " << -y << " DrawText" << endl;
}


void GPostscriptPort::GetPrintingRect (GRect &r)
{
    // A4, with 1" margin
    r.SetRect (0, 0, 595-144, 842-144);
}


void GPostscriptPort::SetCurrentFont (GBaseFont &font)
{
    std::string face = font.GetName();
    if (font.IsBold() || font.IsItalic())
    {
		face += "-";
        if (font.IsBold())
        	face += "Bold";
        if (font.IsItalic())
        	face += "Italic";
    }
/*
	// Duh -- need to do this earlier, perhaps scan the list of
    // fonts already created and output those...
	// Store this font in the list of fonts we need for our document
    int found = DocumentFonts.find_first_of (face, 0);
    if ((found < 0) || (found > DocumentFonts.length()))
    {
    	if (DocumentFonts.length() > 0)
        	DocumentFonts += ", ";
		DocumentFonts += face;
    }
*/
    portStream << endl;
    portStream << "/" << face << " findfont" << endl;
    portStream << font.GetSize () << " scalefont" << endl;
    portStream << "setfont" << endl;
    portStream << endl;
}


// Mac
// Win
// Postscript

void GPostscriptPort::SetPenWidth (int w)
{
    PenWidth = w;
    portStream << w << " setlinewidth" 						<< endl;
    portStream << endl;
}

void GPostscriptPort::AddPSCommands(char* c)
{
	portStream << c;
}

void GPostscriptPort::StartPicture (char *pictFileName)
{
    portStream.open (pictFileName);
        

    // Postscript header
    portStream << "%!PS-Adobe-2.0" 							<< endl;
    portStream << "%%Creator: Roderic D. M. Page" 			<< endl;
    portStream << "%%DocumentFonts: Times-Roman" 		 		<< endl;
    portStream << "%%Title:" <<  pictFileName 				<< endl;
    portStream << "%%BoundingBox: 0 0 595 842" 				<< endl; // A4
    portStream << "%%Pages: 1" 								<< endl;
    portStream << "%%EndComments" 							<< endl;
    portStream << endl;

    // Move origin to top left corner
    portStream << "0 842 translate" << endl;
    portStream << "72 -72 translate" << endl; // one inch margin

    // Some definitions for drawing lines, etc.

    // Drawline draws text with encaps that project...
    portStream << "% Encapsulate drawing a line" 				<< endl;
    portStream << "%    arguments x1 y1 x2 xy2 width" 		<< endl;
    portStream << "/DrawLine {" 								<< endl;
    portStream << "   gsave" 									<< endl;
    portStream << "   setlinewidth" 							<< endl;
    // We may not always want to set this as it works best with rectangular trees...
//    portStream << "   2 setlinecap"							<< endl;
    portStream << "   0 setgray" 								<< endl;
    portStream << "   moveto" 								<< endl;
    portStream << "   lineto" 								<< endl;
    portStream << "   stroke" 								<< endl;
    portStream << "   grestore" 								<< endl;
    portStream << "   } bind def" 							<< endl;
    portStream << endl;

    portStream << "% Encapsulate drawing text" 				<< endl;
    portStream << "%    arguments x y text" 					<< endl;
    portStream << "/DrawText {" 								<< endl;
    portStream << "  gsave 1 setlinewidth 0 setgray" 			<< endl;
    portStream << "  moveto" 									<< endl;
    portStream << "  show grestore" 							<< endl;
	portStream << "} bind def" 								<< endl;
    portStream << endl;

}

void GPostscriptPort::EndPicture ()
{
    portStream << "showpage" 									<< endl;
    portStream << "%%Trailer" 								<< endl;
    portStream << "%%end" 									<< endl;
    portStream << "%%EOF" 									<< endl;
    portStream.close ();
}



#if GPORT_MAC
// Macintosh
void GMacPort::BeginGroup ()
{
//	::PicComment (picGrpBeg, 0, NULL);
}

void GMacPort::EndGroup ()
{
//	::PicComment (picGrpEnd, 0, NULL);
}
#endif

SVGPort::SVGPort ()
{ 
    fontString = "font-family:Times;font-size:12"; 
    DisplayRect.SetRect (0, 0, 400, 400);
}

void SVGPort::DrawLine (const int x1, const int y1, const int x2, const int y2)
{
    portStream << "<path style=\"stroke:black;";
    portStream << ";stroke-width:" << PenWidth;
    portStream << ";stroke-linecap:square";
    portStream << "\" "; 
    
    
    portStream << "d=\"M";
    portStream << x1 << " " << y1 << " " << x2 << " " << y2 << "\"/>";


//    portStream << "<g style=\"fill:none;stroke:black\"><path d=\"M";
//    portStream << x1 << " " << y1 << " " << x2 << " " << y2 << "\"/>";
//    portStream << "</g>" << endl;
}

void SVGPort::DrawCircle (const GPoint &pt, const int radius)
{
/*	portStream << "newpath" << endl;
	portStream << pt.GetX() << " " << -pt.GetY() << " " << radius << " 0 360 arc" 		<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
*/
}


void SVGPort::DrawRect (const GRect &r)
{
/*	portStream << r.GetLeft() << " " << -r.GetTop() << " moveto" 	<< endl;
	portStream << r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << -r.GetHeight() << " rlineto" 			<< endl;
	portStream << -r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << r.GetHeight() << " rlineto"				<< endl;
	portStream << "closepath" 									<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
*/
}

void SVGPort::DrawText (const int x, const int y, const char *text)
{
    portStream << "<text x=\"" << x << "\" y=\"" << y << "\" style=\"" << fontString << "\" >"
            << text << "</text>" << endl;
}


void SVGPort::StartPicture (char *pictFileName)
{
    portStream.open (pictFileName);
    
    portStream << "<?xml version=\"1.0\" ?>" << endl;

	portStream << "<svg xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
			<< "xmlns=\"http://www.w3.org/2000/svg\" "
			<< " width=\"" << DisplayRect.GetWidth() << "px\" "
        	<< "height=\"" << DisplayRect.GetHeight() << "px\" >" << endl;
        
// <title>test</title>
// <desc>test</desc>
}

void SVGPort::EndPicture ()
{
    portStream << "</svg>" << endl;
    portStream.close ();
}


void SVGPort::GetPrintingRect (GRect &r)
{
    r = DisplayRect;
}

void SVGPort::SetCurrentFont (GBaseFont &font)
{
    fontString = "";
    fontString += "font-family:";
    fontString += font.GetName();
	
    char buf[32];
    sprintf (buf, ";font-size:%dpx", font.GetSize());
    fontString += buf;
    
    if (font.IsItalic())
    {
        fontString += ";font-style:italic";
    }
    if (font.IsBold())
    {
        fontString += ";font-weight:bold";
    }

}

MVGPort::MVGPort ()
{ 
	// for now
	fontString = "font '/users/rpage/Sites/CLARITY_.TTF'\n";
	fontString += "font-size 12\n";
	fontString += "text-antialias 0\n"; 
    DisplayRect.SetRect (0, 0, 400, 400);
}

void MVGPort::DrawLine (const int x1, const int y1, const int x2, const int y2)
{
	portStream << "line " <<  x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
}

void MVGPort::DrawCircle (const GPoint &pt, const int radius)
{
/*	portStream << "newpath" << endl;
	portStream << pt.GetX() << " " << -pt.GetY() << " " << radius << " 0 360 arc" 		<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
*/
}


void MVGPort::DrawRect (const GRect &r)
{
/*	portStream << r.GetLeft() << " " << -r.GetTop() << " moveto" 	<< endl;
	portStream << r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << -r.GetHeight() << " rlineto" 			<< endl;
	portStream << -r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << r.GetHeight() << " rlineto"				<< endl;
	portStream << "closepath" 									<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
*/
}

void MVGPort::DrawText (const int x, const int y, const char *text)
{
	portStream << "text " << x << " " << y << "\"" << text << "\"" << std::endl;
}


void MVGPort::StartPicture (char *pictFileName)
{
    portStream.open (pictFileName);
	
	portStream << "viewbox 0 0 "  << DisplayRect.GetWidth() << " " << DisplayRect.GetHeight() <<  std::endl;
	portStream << fontString;
}

void MVGPort::EndPicture ()
{
     portStream.close ();
}


void MVGPort::GetPrintingRect (GRect &r)
{
    r = DisplayRect;
}

void MVGPort::SetCurrentFont (GBaseFont &font)
{
	// for now
	fontString = "font '/users/rpage/Sites/CLARITY_.TTF'\n";
	fontString += "font-size 12\n";
	fontString += "text-antialias 0\n";

/*    fontString = "";
    fontString += "font-family:";
    fontString += font.GetName();
	
    char buf[32];
    sprintf (buf, ";font-size:%dpx", font.GetSize());
    fontString += buf;
    
    if (font.IsItalic())
    {
        fontString += ";font-style:italic";
    }
    if (font.IsBold())
    {
        fontString += ";font-weight:bold";
    }
*/
}

ImageMapPort::ImageMapPort ()
{ 
//    fontString = "font-family:Times;font-size:12"; 
    DisplayRect.SetRect (0, 0, 400, 400);
}

void ImageMapPort::DrawLine (const int x1, const int y1, const int x2, const int y2)
{
//	portStream << "line " <<  x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
}

void ImageMapPort::DrawCircle (const GPoint &pt, const int radius)
{
/*	portStream << "newpath" << endl;
	portStream << pt.GetX() << " " << -pt.GetY() << " " << radius << " 0 360 arc" 		<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
*/
}


void ImageMapPort::DrawRect (const GRect &r)
{
/*	portStream << r.GetLeft() << " " << -r.GetTop() << " moveto" 	<< endl;
	portStream << r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << -r.GetHeight() << " rlineto" 			<< endl;
	portStream << -r.GetWidth() << " 0 rlineto" 					<< endl;
	portStream << "0 " << r.GetHeight() << " rlineto"				<< endl;
	portStream << "closepath" 									<< endl;
	portStream << "stroke" 										<< endl;
	portStream << endl;
*/
}

void ImageMapPort::DrawText (const int x, const int y, const char *text)
{
	// x and y specify the bottom left corner of the text string
	
	portStream
		<< "<area shape=\"rect\""
		<< " coords=\"" << x << "," << (y-6) << "," << (x + strlen(text) * 5) << "," << y << "\""
		<< " href=\"http://www.google.com/search?q=" << text << "\" target=\"_new\""
		<< " title = \"" << text << "\""
		<< " onMouseOver=\"writeText('" << text << "')\""
		<< " onMouseOut=\"writeText('')\""
		<< "/>" << std::endl;
}


void ImageMapPort::StartPicture (char *pictFileName)
{
    portStream.open (pictFileName);
	portStream << "<map name=\"leaves\">" << std::endl;
}

void ImageMapPort::EndPicture ()
{
	portStream << "</map>"  << std::endl;
     portStream.close ();
}


void ImageMapPort::GetPrintingRect (GRect &r)
{
    r = DisplayRect;
}

void ImageMapPort::SetCurrentFont (GBaseFont &font)
{

/*    fontString = "";
    fontString += "font-family:";
    fontString += font.GetName();
	
    char buf[32];
    sprintf (buf, ";font-size:%dpx", font.GetSize());
    fontString += buf;
    
    if (font.IsItalic())
    {
        fontString += ";font-style:italic";
    }
    if (font.IsBold())
    {
        fontString += ";font-weight:bold";
    }
*/
}


