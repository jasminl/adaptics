#include "visio/Buf.h"

namespace visio
{
Buffer::Buffer(unsigned int vsize, unsigned int hsize)
{
	_buffer0 = new double[vsize*hsize];
	_buffer1 = new double[vsize*hsize];
	_buffer2 = new double[vsize*hsize];
	_buffer3 = new double[vsize*hsize];
	_buffer4 = new double[vsize*hsize];
	_buffer5 = new double[vsize*hsize];
	_buffer6 = new double[vsize*hsize];
	_buffer7 = new double[vsize*hsize];
	_buffer8 = new double[vsize*hsize];
	_buffer9 = new double[vsize*hsize];
	_buffer10= new double[vsize*hsize];
	_buffer11= new double[vsize*hsize];
	_buffer12= new double[vsize*hsize];
	_buffer13= new double[vsize*hsize];
	_buffer14= new double[vsize*hsize];
	_buffer15= new double[vsize*hsize];
	_buffer16= new double[vsize*hsize];
	_buffer17= new double[vsize*hsize];
	_buffer18= new double[vsize*hsize];
	_buffer19= new double[vsize*hsize];
	_buffer20= new double[vsize*hsize];
	_buffer21= new double[vsize*hsize];
	_buffer22= new double[vsize*hsize];
	_buffer23= new double[vsize*hsize];
	_buffer24= new double[vsize*hsize];

#if NUM_THREADS == 4 || NUM_THREADS == 16

	_buffer25= new double[vsize*hsize];
	_buffer26= new double[vsize*hsize];
	_buffer27= new double[vsize*hsize];
	_buffer28= new double[vsize*hsize];
	_buffer29= new double[vsize*hsize];
	_buffer30= new double[vsize*hsize];
	_buffer31= new double[vsize*hsize];

#endif

#if NUM_THREADS == 16

	_buffer32= new double[vsize*hsize];
	_buffer33= new double[vsize*hsize];
	_buffer34= new double[vsize*hsize];
	_buffer35= new double[vsize*hsize];
	_buffer36= new double[vsize*hsize];
	_buffer37= new double[vsize*hsize];
	_buffer38= new double[vsize*hsize];
	_buffer39= new double[vsize*hsize];
	_buffer40= new double[vsize*hsize];
	_buffer41= new double[vsize*hsize];
	_buffer42= new double[vsize*hsize];
	_buffer43= new double[vsize*hsize];
	_buffer44= new double[vsize*hsize];
	_buffer45= new double[vsize*hsize];
	_buffer46= new double[vsize*hsize];
	_buffer47= new double[vsize*hsize];
	_buffer48= new double[vsize*hsize];
	_buffer49= new double[vsize*hsize];
	_buffer50= new double[vsize*hsize];
	_buffer51= new double[vsize*hsize];
	_buffer52= new double[vsize*hsize];
	_buffer53= new double[vsize*hsize];
	_buffer54= new double[vsize*hsize];
	_buffer55= new double[vsize*hsize];
	_buffer56= new double[vsize*hsize];
	_buffer57= new double[vsize*hsize];
	_buffer58= new double[vsize*hsize];
	_buffer59= new double[vsize*hsize];
	_buffer60= new double[vsize*hsize];
	_buffer61= new double[vsize*hsize];
	_buffer62= new double[vsize*hsize];
	_buffer63= new double[vsize*hsize];
	_buffer64= new double[vsize*hsize];
	_buffer65= new double[vsize*hsize];
	_buffer66= new double[vsize*hsize];
	_buffer67= new double[vsize*hsize];
	_buffer68= new double[vsize*hsize];
	_buffer69= new double[vsize*hsize];
	_buffer70= new double[vsize*hsize];
	_buffer71= new double[vsize*hsize];
#endif
}

Buffer::~Buffer(void)
{
	delete[] _buffer0 ;
	delete[] _buffer1 ;
	delete[] _buffer2 ;
	delete[] _buffer3 ;
	delete[] _buffer4 ;
	delete[] _buffer5 ;
	delete[] _buffer6 ;
	delete[] _buffer7 ;
	delete[] _buffer8 ;
	delete[] _buffer9 ;
	delete[] _buffer10;
	delete[] _buffer11;
	delete[] _buffer12;
	delete[] _buffer13;
	delete[] _buffer14;
	delete[] _buffer15;
	delete[] _buffer16;
	delete[] _buffer17;
	delete[] _buffer18;
	delete[] _buffer19;
	delete[] _buffer20;
	delete[] _buffer21;
	delete[] _buffer22;
	delete[] _buffer23;
	delete[] _buffer24;

#if NUM_THREADS == 4 || NUM_THREADS == 16

	delete[] _buffer25;
	delete[] _buffer26;
	delete[] _buffer27;
	delete[] _buffer28;
	delete[] _buffer29;
	delete[] _buffer30;
	delete[] _buffer31;

#endif

#if NUM_THREADS == 16

	delete[] _buffer32;
	delete[] _buffer33;
	delete[] _buffer34;
	delete[] _buffer35;
	delete[] _buffer36;
	delete[] _buffer37;
	delete[] _buffer38;
	delete[] _buffer39;
	delete[] _buffer40;
	delete[] _buffer41;
	delete[] _buffer42;
	delete[] _buffer43;
	delete[] _buffer44;
	delete[] _buffer45;
	delete[] _buffer46;
	delete[] _buffer47;
	delete[] _buffer48;
	delete[] _buffer49;
	delete[] _buffer50;
	delete[] _buffer51;
	delete[] _buffer52;
	delete[] _buffer53;
	delete[] _buffer54;
	delete[] _buffer55;
	delete[] _buffer56;
	delete[] _buffer57;
	delete[] _buffer58;
	delete[] _buffer59;
	delete[] _buffer60;
	delete[] _buffer61;
	delete[] _buffer62;
	delete[] _buffer63;
	delete[] _buffer64;
	delete[] _buffer65;
	delete[] _buffer66;
	delete[] _buffer67;
	delete[] _buffer68;
	delete[] _buffer69;
	delete[] _buffer70;
	delete[] _buffer71;

#endif
}
}
