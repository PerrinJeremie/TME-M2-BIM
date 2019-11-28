#ifndef __XNEXUS_H
#define __XNEXUS_H

//
// XNexus exception class
//
class XNexus
{
public:
	nxsstring msg;
	std::streampos pos;
	long line;
	long col;

	XNexus( nxsstring s, std::streampos fp = 0, long fl = 0L, long fc = 0L );
};

#endif
