#pragma once

#include <stdexcept>
#include "SparseFrame.h"

namespace visio
{
/**
 *	Input Layer
 */
template <typename T>
class Input
{
public:
	enum{USE_MOTION, NO_USE_MOTION};	//!<Whether or not to use motion
	enum{USE_V1, NO_USE_V1};			//!<Whether or not to use V1 input
	enum{USE_V2, NO_USE_V2}; 			//!<Whether or not to use V2 input
	enum{USE_ATTEN, NO_USE_ATTEN}; 		//!<Whether or not to use attention

public:
	/**
	 *	Constructor
	 *	@param useMotion Determines whether to input motion signals or not.
	 *	@param useV2     Determines whether to input V2 boundaries or not (default false).
	 *	@param useV1     Determines whether to input V1 boundaries or not(default false).
	 */
	Input(int use_motion, int use_V2, int use_V1, int use_attention);

	/**
	 *	Destructor
	 */
	virtual
	~Input<T>();

	/**
	 *	Shows the current frame on stdout
	 *	\param field The name of the field to display {"ONMOT","OFFMOT","V2D1","V2D2","V2D3"}.
	 */
	void show_frame(const char* field) const;

	/**
	 *	Indicates the number of frames loaded.
	 *	\return The number of frames.
	 */
	unsigned int nb_frames() const
	{
		return _nb_frames;
	}

	/**
	 *	Gives the width of the stimulus frames.
	 *	\return The width of the frame.
	 */
	unsigned int width() const
	{
		return _frame_width;
	}

	/**
	 *	Gives the height of the stimulus frames
	 *	\return The height of the frame.
	 */
	unsigned int height() const
	{
		return _frame_height;
	}

	/**
	 *	Returns the number of cycles per second for current stimulus
	 *	\return Number of cycles per second.
	 */
	double cycle_per_second() const
	{
		return _cycle_per_second;
	}

	/**
	 *	Imports an xml file of frames stored sparsely, must be called once at beginning.
	 *	\param inputFile Name of file to load.
	 */
	void import(const std::string& inputFile);

	/**
	 *	To access ON motion input.
	 *	\return A vector containing decoded motion input in the ON channel.
	 */
	const std::vector<T>& on_motion() const
	{
		return _on_motion;
	}

	/**
	 *	To access OFF motion input.
	 *	\return A vector containing decoded motion input in the OFF channel.
	 */
	const std::vector<T>& off_motion() const
	{
		return _off_motion;
	}

	/**
	 *	To access V2 input at each depth
	 *	\param depth The depth plane {1,2,3} for which to access form input.
	 *	\return Decoded form input at the specified depth.
	 */
	std::vector<T> V2(unsigned int depth);

	/**
	 *	Exports current V2 depth plane into array
	 *	\param depth The depth plane to process {1,2,3}.
	 *	\param array The array in which to store decoded values.
	 *	\return Pointer to the decoded array.
	 */
	double* V2(unsigned int depth, double*& array);

	/**
	 *	Exports current attentional top down input
	 *	\param direction Motion direction for which to get attentional input.
	 *	\param array Vector in which to store the attentional input parameters.
	 */
	void attention(unsigned int direction, std::vector<double>& array);

	/**
	 *	\return The number of V2 planes.
	 */
	unsigned int nb_V2_planes() const
	{
		return _nb_V2_planes;
	}

	/**
	 *	To go to next stored frame.
	 */
	void next_frame()
	{
		++_cur_frame_nb;	//Go to next frame
		load();				//Load that frame
	}

	/**
	 *	Loads currently pointed-to frame.
	 */
	void load();

	/**
	 *	Indicates whether there is attentional top down in that time frame and direction.
	 *	\param d Motion direction to probe.
	 *	\return Boolean value indicating whether the motion direction is being attended.
	 */
	bool attend_direction(unsigned int d);

private:

	/**
	 *	Allocate buffer arrays to store current frame.
	 *	\param height Frame height
	 *	\param width  Frame width
	 */
	void create_buffer(unsigned int height, unsigned int width);

	/**
	 *	Decodes an input frame stored sparsely.
	 *	\param Input (encoded) frame.
	 *	\param Output (decoded) frame.
	 */
	void decode(std::vector<T>*& source, std::vector<T>& end);

	//todo: why is the latter needed?
	/**
	 *	Decodes a frame sparsely stored as double values.
	 *	\param Input (encoded) frame.
	 *	\param Output (decoded) frame.
	 */
	void decode_double(const std::vector<double>* source, std::vector<double>& end);

private:

	/** @name Motion input
		Input arrays for ON and OFF channels.
	*/
	//@{
	std::vector<T> _on_motion;
	std::vector<T> _off_motion;
	//@}

	/** @name V2 input
		Input arrays for V2 boundaries at different depths/scales
	*/
	std::vector<double> _V2D1;     //!< V2 boundaries, depth 1 (largest scale)
	std::vector<double> _V2D2;	   //!< V2 boundaries, depth 2 (smallest scale in 2-scales simulations, otherwise middle scale)
	std::vector<double> _V2D3;	   //!< V2 boundaries, depth 3 (nearest, although this depth may not be used)

	/** @name Top down attentional parameters (scale 1)
		Vectors storing the state of attention top down in each motion direction.
	*/
	//@{
	std::vector<double> _attd0s1;
	std::vector<double> _attd1s1;
	std::vector<double> _attd2s1;
	std::vector<double> _attd3s1;
	std::vector<double> _attd4s1;
	std::vector<double> _attd5s1;
	std::vector<double> _attd6s1;
	std::vector<double> _attd7s1;
	//@}

	/** @name Motion input to V1
		Input arrays to Level 1 of form system.
	*/
	//@{
	std::vector<T> _on_V1;
	std::vector<T> _off_V1;
	//@}

	std::vector<SparseFrame<T>>* _frame;	//!< Sparsely encoded array of frames

	SparseFrame<T>* _current_frame;			//!< The currently loaded frame

	/** @name Import flags
     * Boolean flags that indicate which input to load.
	 */
	//@{
	int _import_motion;				//!<Indicates whether to load motion or not
	int _import_V2;					//!<Indicates whether to load V2 or not
	int _import_V1;					//!<Indicates whether to load V1 or not
	int _import_attention;			//!<Indicates whether to load attention or not
	//@}

	unsigned int _nb_frames;		//!<The total number of frames
	unsigned int _frame_width;		//!<Input frame width in pixels
	unsigned int _frame_height;		//!<Input frame height in pixels
	unsigned int _nb_V2_planes;		//!<The total number of depth planes
	double _cycle_per_second;		//!<The number of cycles of the stimulus per second

	unsigned int _cur_frame_nb;		//!< The current Frame number
};

template<typename T>
Input<T>::Input(int use_motion, int use_V2, int use_V1, int use_attention)
: _import_motion(use_motion), _import_V2(use_V2), _import_V1(use_V1), _import_attention(use_attention)
  {};

template <typename T>
Input<T>::~Input()
{
	delete _frame;
}

template <typename T>
void Input<T>::show_frame(const char* field) const
{
	const std::vector<T>* p;
	if (!strcmp("ONMOT",field))
		p = &_on_motion;
	else if (!strcmp("OFFMOT",field))
		p = &_off_motion;
	else if(!strcmp("V2D1",field))
		p = &_V2D1;
	else if(!strcmp("V2D2",field))
		p = &_V2D2;
	else if(!strcmp("V2D3",field))
		p = &_V2D3;
	for(int i = 0; i < _frame_height; i++)
	{
		for(int j = 0; j < _frame_width; j++)
			std::cout<<p->at(i*_frame_height + j)<<" ";
		std::cout<<endl;
	}
}

template<typename T>
void Input<T>::create_buffer(unsigned int height, unsigned int width)
{
	if(_import_motion == USE_MOTION)
	{
		_on_motion = std::vector<T>(height * width);
		fill(_on_motion.begin(), _on_motion.end(), 0);

		_off_motion = std::vector<T>(height * width);
		fill(_off_motion.begin(), _off_motion.end(), 0);
	}

	if(_import_V2 == USE_V2)
	{
		switch(_nb_V2_planes)
		{
			case 1:
				_V2D1 = std::vector<double>(height * width);
				fill(_V2D1.begin(), _V2D1.end(), 0);
				break;
			case 2:
			{
				_V2D1 = std::vector<double>(height * width);
				_V2D2 = std::vector<double>(height * width);
				fill(_V2D1.begin(), _V2D1.end(), 0);
				fill(_V2D2.begin(), _V2D2.end(), 0);
				break;
			}
			case 3:
			{
				_V2D1 = std::vector<double>(height * width);
				_V2D2 = std::vector<double>(height * width);
				_V2D3 = std::vector<double>(height * width);
				fill(_V2D1.begin(), _V2D1.end(), 0);
				fill(_V2D2.begin(), _V2D2.end(), 0);
				fill(_V2D3.begin(), _V2D3.end(), 0);
				break;
			}
			default:
				throw std::runtime_error("input<T>::createBuffer: wrong number of planes");
		}
	}
}

template<typename T>
void Input<T>::import(const std::string& input_file)
{
	ifstream file(input_file);
	if(file.bad())
		throw std::ios_base::failure("input<T>::import: file " + input_file + " doesn't exist");

	TiXmlDocument doc(input_file.c_str());

	if(!doc.LoadFile())
		throw std::runtime_error("input<T>::import: unable to parse file " + input_file);
	else
	{
		TiXmlElement* root = doc.RootElement();
		if(!root->NoChildren())
		{
			_cycle_per_second = atof(root->FirstChildElement("parameters")->Attribute("cyclePerSecond"));
			_nb_frames = atoi(root->FirstChildElement("frame")->Attribute("TotNumberFrames"));
			_frame_width  = atoi(root->FirstChildElement("frame")->Attribute("ImageWidth"));
			_frame_height = atoi(root->FirstChildElement("frame")->Attribute("ImageHeight"));
			_nb_V2_planes  = atoi(root->FirstChildElement("frame")->Attribute("DepthPlanes"));
			_frame = new std::vector< SparseFrame<T> >(_nb_frames);

			create_buffer(_frame_height, _frame_width);

			auto pFrame = _frame->begin();

			TiXmlNode* child = nullptr;

			//Get first frame element
			child = root->FirstChildElement("parameters");

			while((child = root->IterateChildren(child)))
			{
				if(_import_motion == USE_MOTION)
				{
					//On motion channel
					TiXmlNode* onM = child->FirstChildElement("OnMOtion");

					//Get number of points and allocate proper space
					int nbOnPt = 0;
					onM->ToElement()->Attribute("NbPoints",&nbOnPt);
					pFrame->get_triplets(pFrame->on_motion(), onM, nbOnPt, _frame_width);

					//Off motion Channel
					TiXmlNode* offM = child->FirstChildElement("OffMOtion");

					//Get number of points and allocate proper space
					int nbOffPt = 0;
					offM->ToElement()->Attribute("NbPoints", &nbOffPt);

					pFrame->get_triplets(pFrame->off_motion(), offM, nbOffPt, _frame_width);
				}

				if(_import_V2 == USE_V2)
				{
					TiXmlNode* v2 = child->FirstChildElement("V2");

					TiXmlNode* v2Depth = nullptr;

					while((v2Depth = v2->IterateChildren(v2Depth)))
					{
						//For each depth plane
						int nbV2Pt = 0;
						v2Depth->ToElement()->Attribute("NbPoints",&nbV2Pt);

						int depth = 0;
						v2Depth->ToElement()->Attribute("Plane",&depth);

						string currentDepth;
						switch(depth)
						{
						case 1:
							currentDepth = "V2D1";
							break;
						case 2:
							currentDepth = "V2D2";
							break;
						case 3:
							currentDepth = "V2D3";
							break;
						default:
							throw std::runtime_error("input<T>::import: unsupported depth plane");
						}
						pFrame->get_triplets_double(pFrame->array_double(currentDepth.c_str()), v2Depth, nbV2Pt, _frame_width);
					}
				}

				if(_import_V1 == USE_V1)
				{
					//Access V1 data
				}

				if(_import_attention == USE_ATTEN)
				{
					TiXmlNode* att = child->FirstChildElement("att");

					if(att == nullptr)
					{
						//There's no top down attention in this stimulus or frame
					}
					else
					{
						//There's at least one att element
						TiXmlNode* attItem = nullptr;

						while((attItem = child->IterateChildren("att", attItem)))
						{
							int direction = 0;

							attItem->ToElement()->Attribute("d", &direction);
							std::vector<double>* attParam = pFrame->create_att_array((unsigned)direction);

							attItem->ToElement()->Attribute("s", &attParam->at(0));
							attItem->ToElement()->Attribute("sig", &attParam->at(1));
							attItem->ToElement()->Attribute("x", &attParam->at(2));
							attItem->ToElement()->Attribute("y", &attParam->at(3));
						}
					}
				}
				pFrame++;
			}

			//Here set current frame to nonexisting frame
			_cur_frame_nb = -1;
			return;
		}
		else
			throw std::runtime_error("input<T>::import: document root element has no children nodes");
	}
}

template<typename T>
std::vector<T> Input<T>::V2(unsigned int depth)
{
	switch(depth)
	{
	case 1:
		return _V2D1;
		break;
	case 2:
		return _V2D2;
		break;
	case 3:
		return _V2D3;
	default:
		throw std::runtime_error("input<T>::V2: invalid depth plane");
	}
}


template <typename T>
double* Input<T>::V2(unsigned int depth, double*& array)
{
	std::vector<double>* v = nullptr;
	switch(depth)
	{
	case 1:
		v = &_V2D1;
		break;
	case 2:
		v = &_V2D2;
		break;
	case 3:
		v = &_V2D3;
		break;
	default:
		throw std::runtime_error("input<T>::V2: invalid depth plane");
	}

	double* q = array;

	for(auto p = v->begin();p != v->end(); p++,q++)
		*q = *p;
	return array;
}

template<typename T>
void Input<T>::attention(unsigned int direction, std::vector<double>& array)
{
	switch(direction)
	{
	case 0: array = _attd0s1;
		break;
	case 1: array = _attd1s1;
		break;
	case 2: array = _attd2s1;
		break;
	case 3: array = _attd3s1;
		break;
	case 4: array = _attd4s1;
		break;
	case 5: array = _attd5s1;
		break;
	case 6: array = _attd6s1;
		break;
	case 7: array = _attd7s1;
		break;
	default:
		throw std::runtime_error("input<T>::attention: invalid direction");
	}
}

template<typename T>
void Input<T>::load()
{
	//Reset current input
	if(_import_motion == USE_MOTION)
	{
		fill(_on_motion.begin(), _on_motion.end(), 0);
		fill(_off_motion.begin(), _off_motion.end(), 0);

		//Transfer from sparsely coded frame to full input
		std::vector<T>* p = _frame->at(_cur_frame_nb).array("ONMOT");
		std::vector<T>* q = _frame->at(_cur_frame_nb).array("OFFMOT");

		decode(p,_on_motion);
		decode(q,_off_motion);
	}

	if(_import_V2 == USE_V2)
	{
		switch(_nb_V2_planes)
		{
		case 1:
			{
				fill(_V2D1.begin(),_V2D1.end(),0);
				std::vector<double>* p = _frame->at(_cur_frame_nb).array_double("V2D1");
				decode_double(p,_V2D1);
				break;
			}
		case 2:
			{
				fill(_V2D1.begin(),_V2D1.end(),0);
				fill(_V2D2.begin(),_V2D2.end(),0);
				std::vector<double>* p = _frame->at(_cur_frame_nb).array_double("V2D1");
				std::vector<double>* q = _frame->at(_cur_frame_nb).array_double("V2D2");
				decode_double(p,_V2D1);
				decode_double(q,_V2D2);
				break;
			}
		case 3:
			{
				fill(_V2D1.begin(),_V2D1.end(),0);
				fill(_V2D2.begin(),_V2D2.end(),0);
				fill(_V2D3.begin(),_V2D3.end(),0);
				std::vector<double>* p = _frame->at(_cur_frame_nb).array_double("V2D1");
				std::vector<double>* q = _frame->at(_cur_frame_nb).array_double("V2D2");
				std::vector<double>* r = _frame->at(_cur_frame_nb).array_double("V2D3");
				decode_double(p,_V2D1);
				decode_double(q,_V2D2);
				decode_double(r,_V2D3);
				break;
			}
		default:
			return;
		}
	}

	if(_import_attention == USE_ATTEN)
	{
		_attd0s1 = *_frame->at(_cur_frame_nb).attd0s1();
		_attd1s1 = *_frame->at(_cur_frame_nb).attd1s1();
		_attd2s1 = *_frame->at(_cur_frame_nb).attd2s1();
		_attd3s1 = *_frame->at(_cur_frame_nb).attd3s1();
		_attd4s1 = *_frame->at(_cur_frame_nb).attd4s1();
		_attd5s1 = *_frame->at(_cur_frame_nb).attd5s1();
		_attd6s1 = *_frame->at(_cur_frame_nb).attd6s1();
		_attd7s1 = *_frame->at(_cur_frame_nb).attd7s1();
	}
}

template<typename T>
bool Input<T>::attend_direction(unsigned int d)
{
	switch(d)
	{
		case 0: return !_attd0s1.empty();
		case 1: return !_attd1s1.empty();
		case 2: return !_attd2s1.empty();
		case 3: return !_attd3s1.empty();
		case 4: return !_attd4s1.empty();
		case 5: return !_attd5s1.empty();
		case 6: return !_attd6s1.empty();
		case 7: return !_attd7s1.empty();
	}
	return false;
}

template<typename T>
void Input<T>::decode(std::vector<T>*& source, std::vector<T>& end)
{
	auto q = end.begin();
	for(auto p = source->begin(); p != source->end(); p++)
	{
		q += *p++;	//Shift to current position
		*q = *p;
	}
}

template<typename T>
void Input<T>::decode_double(const std::vector<double>* source, std::vector<double>& end)
{
	auto q = end.begin();
	for(auto p = source->begin(); p != source->end(); p++)
	{
		q += (int)*p++;	//Shift to current position
		*q = *p;
	}
}
}
