#include "mli.h"
#include "engine.h"

#include <iostream>
using namespace std;

Engine * ml;


mli::mli(unsigned int init, unsigned int maxnb):m_nb(init),m_maxNb(maxnb)
{
	ml = engOpen(NULL);
	if(ml==NULL)
	{
		cerr<<"Can't open matlab"<<endl;
			exit(-1);
	}

	string a = "cd 'C:\\Documents and Settings\\jasminl\\My Documents\\Brian_project\\videos\\xmas'";
	engEvalString(ml,a.c_str());
	engEvalString(ml,"addpath('C:\\Documents and Settings\\jasminl\\My Documents\\Brian_project\\code\\track')");
}

mli::~mli(void)
{
	if (ml != NULL)
	{
		engClose(ml);
	}
}

void mli::loadMovie(std::string filename)
{
	string a;
	
	//Get file info
	a = "is = aviinfo('";
	a += filename;
	a += "')";
	engEvalString(ml,a.c_str());
	
	//Load all avi file
	a = "nbFrames = is.NumFrames - 1;";
	engEvalString(ml,a.c_str());
	a = "mov = aviread('";
	a += filename;
	a += "',1:nbFrames)";
	if (0!= engEvalString(ml,a.c_str()))
	{
		cout<<"Couldn't load movie\n";
		exit(-1);
	}
	
}
void mli::getFrame(string nb, unsigned char*& frame, unsigned int& width, unsigned int& height)
{
	//Delete frame if previously allocated
	delete frame;

	//Get frame
	string a = "x = mov(";
	a += nb;
	a += ").cdata;";
	engEvalString(ml,a.c_str());

	//Get size
	a = "[height,width,depth] = size(x);";
	engEvalString(ml,a.c_str());
	mxArray* w = engGetVariable(ml,"width");
	mxArray* h = engGetVariable(ml,"height");
	double wb,hb;
	memcpy((char*)&wb,mxGetPr(w),sizeof(double));
	memcpy((char*)&hb,mxGetPr(h),sizeof(double));
	width = static_cast<unsigned int>(wb);
	height = static_cast<unsigned int>(hb);

	//Get data (this should be 3D array of unsigned char)
	a = "y = interleaveImage(x);";
	engEvalString(ml,a.c_str());	//Reshape
	mxArray* dat = engGetVariable(ml,"y");

	frame = new unsigned char[width * height * 3];
	memcpy((unsigned char*)frame ,(unsigned char*) mxGetData(dat),sizeof(unsigned char) * width * height * 3);

	//Clear up
	mxDestroyArray(w);
	mxDestroyArray(h);
	mxDestroyArray(dat);


}

void mli::getTarget(double& x, double& y, double& hx, double& hy, unsigned char*& roi, unsigned char*& bw)
{
	engEvalString(ml,"[center, bw, bb] = getTarget(x);");
	engEvalString(ml,"xc = center(1); yc = center(2);");

	engEvalString(ml,"hx = bb.BoundingBox(3);");           
	engEvalString(ml,"hy = bb.BoundingBox(4);");

	//Get bounding ellipse size
	mxArray* w = engGetVariable(ml,"hx");
	mxArray* h = engGetVariable(ml,"hy");
	memcpy((char*)&hx,mxGetPr(w),sizeof(double));
	memcpy((char*)&hy,mxGetPr(h),sizeof(double));

	//Get center
	mxArray* xc = engGetVariable(ml,"xc");
	mxArray* yc = engGetVariable(ml,"yc");
	memcpy((char*)&x,mxGetPr(xc),sizeof(double));
	memcpy((char*)&y,mxGetPr(yc),sizeof(double));

	//Get ROI (for feature calculations)
	roi = new unsigned char[hx*hy * 3];		//Assum RGB image
	bw  = new unsigned char[hx*hy];			//Boolean values
	engEvalString(ml,"[outImage outBW] = getroi(x,bw,bb);");	
	mxArray* mroi = engGetVariable(ml,"outImage");
	mxArray* mbw  = engGetVariable(ml,"outBW");
	memcpy((unsigned char*)&roi[0],mxGetPr(mroi),sizeof(unsigned char) * 3 * hx * hy);
	memcpy((unsigned char*)&bw[0],mxGetPr(mbw),sizeof(unsigned char) * hx * hy);

	//Clean up
	mxDestroyArray(w);
	mxDestroyArray(h);
	mxDestroyArray(mroi);
	mxDestroyArray(mbw);
}

void mli::showTarget(double x, double y, double hx, double hy, double scale)
{
	engEvalString(ml,"figure(2)");
	engEvalString(ml,"image(x)");

	double majA = hy* scale;
	double minA = hx* scale;
	mxArray* majAxis = mxCreateDoubleMatrix(1, 1, mxREAL); 
	memcpy((char *) mxGetPr(majAxis), (char *)&majA, sizeof(double));
	mxArray* minAxis = mxCreateDoubleMatrix(1, 1, mxREAL); 
	memcpy((char *) mxGetPr(minAxis), (char *)&minA, sizeof(double));

	engPutVariable(ml,"majAxis",majAxis);
	engPutVariable(ml,"minAxis",minAxis);

    double xLoc = x - minA/2;
    double yLoc = y - majA/2;

	mxArray* pxLoc =  mxCreateDoubleMatrix(1, 1, mxREAL); 
	memcpy((char *) mxGetPr(pxLoc), (char *)&xLoc, sizeof(double));
	
	mxArray* pyLoc =  mxCreateDoubleMatrix(1, 1, mxREAL); 
	memcpy((char *) mxGetPr(pyLoc), (char *)&yLoc, sizeof(double));
	engPutVariable(ml,"xLoc",pxLoc);
	engPutVariable(ml,"yLoc",pyLoc);

    engEvalString(ml,"r = rectangle('position',[xLoc yLoc minAxis majAxis]);");
    engEvalString(ml,"set(r,'EdgeColor',[0 1 0],'LineWidth',2);");
}