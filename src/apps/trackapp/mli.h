#pragma once


//This defines the matlab interface

#include "engine.h"
extern Engine * ml;		//Single global object

#include <iostream>
#include "string.h"

class mli
{

protected:
	unsigned int m_nb;		/*<< Current frame number */
	unsigned int m_maxNb;	/*<< Max number of frames to track */

public:
	mli(unsigned int init, unsigned int maxnb);
 
	/**
		Evaluates a command
	*/
	void eval(std::string command)
	{
		engEvalString(ml,command.c_str());
	}

	unsigned int nbIter() const
	{
		return m_maxNb;
	}

	/**
		Returns one frame from the movie
	*/
	void getFrame(std::string nb, unsigned char*& frame, unsigned int& width, unsigned int& height);

	/**
		Loads movie
	*/
	void loadMovie(std::string filename);

	/**
		Get subsequent frames
	*/
	void nextFrame(unsigned char*& frame)
	{
		unsigned int b1,b2;	//Useless buffers
		char b3[10]={0};
		_itoa_s(m_nb,b3,10,10);
		getFrame(std::string(b3),frame,b1,b2);
		m_nb++;	//Go to next frame
	}

	/**
		Obtain target ellipse, and region of interest (for feature calculations)
	*/
	void getTarget(double& x, double& y, double& hx, double& hy, unsigned char*& roi, unsigned char*& bw);

	/**
		Show target candidate
	*/
	void showTarget(double x, double y, double hx, double hy, double scale);

public:
	~mli(void);
};
