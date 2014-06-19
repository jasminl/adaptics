/**
 * \file SparseFrame.h
 * \brief 3D formotion sparsely encoded input frame (decoded from XML)
 */
#pragma once

#include <vector>
#include <stdexcept>

namespace visio
{
/**
 *	\brief Sparsely encoded input frame
 */
template <typename T>
class SparseFrame
{
public:

	/**
	 *	Default constructor
	 */
	SparseFrame();

	/**
	 *	Destructor
	 */
	~SparseFrame();

	/**
	 *	Access currently buffered input by name
	 *	@param aName Can be either of "ONMOT", "OFFMOT", "V2D1", "V2D2", "V2D3"
	 */
	std::vector<T>*& array(const char* aName);

	/**
	 * Access double-precision depth plane by name
	 * @param aName Depth plane by name V2D1, V2D2, V2D3
	 */
	std::vector<double>*& array_double(const char* aName);

	/**	@name Direct buffer access member functions
		Used to access a buffer
	*/
	//@{
	/**
	 * ON channel motion input
	 */
	std::vector<T>*& on_motion(){return _on_motion;}

	/**
	 * OFF channel motion input
	 */
	std::vector<T>*& off_motion(){return _off_motion;}

	/**
	 * V2 depth 1 boundary input
	 */
	std::vector<T>*& v2d1(){return _V2D1;}

	/**
	 * V2 depth 2 boundary input
	 */
	std::vector<T>*& v2d2(){return _V2D2;}

	/** V2 depth 3 boundary input */
	std::vector<T>*& v2d3(){return _V2D3;}
	//@}

	/** @name Functions to access current attentional top-down
	*/
	//@{
	std::vector<double>* attd0s1(){return &_attd0s1;}
	std::vector<double>* attd1s1(){return &_attd1s1;}
	std::vector<double>* attd2s1(){return &_attd2s1;}
	std::vector<double>* attd3s1(){return &_attd3s1;}
	std::vector<double>* attd4s1(){return &_attd4s1;}
	std::vector<double>* attd5s1(){return &_attd5s1;}
	std::vector<double>* attd6s1(){return &_attd6s1;}
	std::vector<double>* attd7s1(){return &_attd7s1;}
	//@}

	/**
	 *	Creates array to store sparse input
	 *	\param pointer Reference to the pointer to the location where to allocate memory for an input array.
	 *	\param nb_points The number of points to allocate space for (2*nbPoints).
	 */
	void create_array(std::vector<T>*& pointer, unsigned int nb_points)
	{	//Allocate array such that each point corresponds to a triplet
		pointer = new std::vector<T>(nb_points * 2);
	}

	/**
	 * Creates array to store sparse input
	 * @param pointer Reference to the pointer to the location where to allocate memory for an input array.
	 * @param nb_points The number of points to allocate space for (2*nbPoints).
	 */
	void create_array(std::vector<double>*& pointer, unsigned int nb_points)
	{	//Allocate array such that each point corresponds to a triplet
		pointer = new std::vector<double>(nb_points * 2);
	}

	/**
	 *	Creates attentional array
	 *	@param d The motion direction ({0,1,2,3,4,5,6,7}) for which to create the 4-element attentional array.
	 *	@return The address of the allocated array.
	 */
	std::vector<double>* create_att_array(unsigned int d);

	/**
	 *	Indicates whether attention in a given direction is used for that time frame
	 *	@param d specifies which direction to attend.
	 */
	bool attend_direction(unsigned int d);

	/**
	 *	Cleans up memory space for a previously allocated array.
	 *	@param pointer Reference to a pointer to the address of the array to free.
	 */
	void delete_array(std::vector<T>*& pointer)
	{
		delete pointer;
	}

	/**
	 *	Cleans up memory space for a previously allocated array.
	 *	@param pointer Reference to a pointer to the address of the array to free.
	 */
	void delete_array(std::vector<double>*& pointer)
	{
		delete pointer;
	}

	/**
	 *	Reads xml triplets and stores them in array.
	 *	@param array The destination where to store the read triplets.
	 *	@param node A pointer to the current frame in the XML file object wherefrom to read triplets.
	 *	@param nbPts The number of points to read.
	 *	@param imageWidth The horizontal extent of the input frame.
	 */
	void get_triplets(std::vector<T>*& array, TiXmlNode* node, unsigned int nbPts, unsigned int imageWidth);

	/**
	 *	Reads xml triplets and stores them in a double-precision array. Necessary when the input data is not integer.
	 *	@param array The destination where to store the read triplets.
	 *	@param node A pointer to the current frame in the XML file object wherefrom to read triplets.
	 *	@param nbPts The number of points to read.
	 *	@param imageWidth The horizontal extent of the input frame.
	 */
	void get_triplets_double(std::vector<double>*& array, TiXmlNode* node, unsigned int nbPts, unsigned int imageWidth);

private:

	/**	@name Motion input
	*/
	std::vector<T>* _on_motion;			//!< Input to the ON motion channel
	std::vector<T>* _off_motion;		//!< Input to the OFF motion channel

	/**	@name V1 input
	*/
	//@{
	std::vector<T>* _on_V1;				//!< Input to the ON V1 channel
	std::vector<T>* _off_V1;			//!< Input to the OFF V1 channel
	//@}

	/** @name V2 boundary input
	*/
	//@{
	std::vector<double>* _V2D1;			//!< V2 boundary input in depth 1
	std::vector<double>* _V2D2;			//!< V2 boundary input in depth 2
	std::vector<double>* _V2D3;			//!< V2 boundary input in depth 3

	//@}

	/** @name Attentional top down
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
};

template<typename T>
SparseFrame<T>::SparseFrame()
: _on_motion(nullptr), _off_motion(nullptr), _on_V1(nullptr), _off_V1(nullptr),
	_V2D1(nullptr), _V2D2(nullptr), _V2D3(nullptr)
{
	_attd0s1.clear();
	_attd1s1.clear();
	_attd2s1.clear();
	_attd3s1.clear();
	_attd4s1.clear();
	_attd5s1.clear();
	_attd6s1.clear();
	_attd7s1.clear();
}

template<typename T>
SparseFrame<T>::~SparseFrame()
{
	 delete_array(_on_motion);
	 delete_array(_off_motion);
	 delete_array(_on_V1);
	 delete_array(_off_V1);

	 delete_array(_V2D1);
	 delete_array(_V2D2);
	 delete_array(_V2D3);
}

template<typename T>
std::vector<T>*& SparseFrame<T>::array(const char* aName)
{
	if(!strcmp(aName,"ONMOT"))	return _on_motion;
	if(!strcmp(aName,"OFFMOT"))	return _off_motion;

	throw std::runtime_error("SparseFrame<T>::array: invalid name " + std::string(aName));
}

template<typename T>
std::vector<double>*& SparseFrame<T>::array_double(const char* aName)
{
	if(!strcmp(aName,"V2D1"))	return _V2D1;
	if(!strcmp(aName,"V2D2"))	return _V2D2;
	if(!strcmp(aName,"V2D3"))	return _V2D3;

	throw std::runtime_error("SparseFrame<T>::array_double: invalid name " + std::string(aName));
}

template<typename T>
std::vector<double>* SparseFrame<T>::create_att_array(unsigned int d)
{
	switch(d)
	{
	case 0:
		_attd0s1 = std::vector<double>(4);
		return &_attd0s1;
	case 1:
		_attd1s1 = std::vector<double>(4);
		return &_attd1s1;
	case 2:
		_attd2s1 = std::vector<double>(4);
		return &_attd2s1;
	case 3:
		_attd3s1 = std::vector<double>(4);
		return &_attd3s1;
	case 4:
		_attd4s1 = std::vector<double>(4);
		return &_attd4s1;
	case 5:
		_attd5s1 = std::vector<double>(4);
		return &_attd5s1;
	case 6:
		_attd6s1 = std::vector<double>(4);
		return &_attd6s1;
	case 7:
		_attd7s1 = std::vector<double>(4);
		return &_attd7s1;
	}

	throw std::runtime_error("SparseFrame<T>::create_att_array: invalid direction");
}

template<typename T>
void SparseFrame<T>::get_triplets(std::vector<T>*& array, TiXmlNode* node, unsigned int nbPts, unsigned int imageWidth)
{
	if(array == nullptr)
		create_array(array, nbPts);

	TiXmlNode* child = nullptr;

	auto p = array->begin();

	int i = 0;
	unsigned int sum = 0;

	while((child = node->IterateChildren(child)))
	{
		if (i == 0)
		{
			//Absolute index
			*p++ = (T)(( atoi(child->ToElement()->Attribute("y") ) -1) * imageWidth
				    + (atoi(child->ToElement()->Attribute("x"))-1) );
		}
		else
		{
			//Difference index
			*p++ = (T)( (atoi(child->ToElement()->Attribute("y") )-1) * imageWidth
				    + (atoi(child->ToElement()->Attribute("x"))-1)  - sum) ;
		}

		sum += *(p-1);

		*p++ = (T)atoi(child->ToElement()->Attribute("value"));

		i++;
	}

}

template<typename T>
void SparseFrame<T>::get_triplets_double(std::vector<double>*& array, TiXmlNode* node, unsigned int nbPts, unsigned int imageWidth)
{
	if(array == nullptr)
		create_array(array,nbPts);

	TiXmlNode* child = nullptr;

	auto p = array->begin();

	int i = 0;
	double sum = 0;

	while((child = node->IterateChildren(child)))
	{
		if(i == 0)
		{
			//Absolute index
			*p++ = (double)(( atof(child->ToElement()->Attribute("y") ) -1) * imageWidth
				    + (atof(child->ToElement()->Attribute("x"))-1) );
		}
		else
		{
			//Difference index
			*p++ = (double)( (atof(child->ToElement()->Attribute("y") )-1) * imageWidth
				    + (atof(child->ToElement()->Attribute("x"))-1)  - sum) ;
		}

		sum += *(p-1);

 		*p++ = (double)atof(child->ToElement()->Attribute("value"));
 		i++;
	}
}

template<typename T>
bool SparseFrame<T>::attend_direction(unsigned int d)
{
	switch (d)
	{
	case 0:
		return !_attd0s1.empty();
	case 1:
		return !_attd1s1.empty();
	case 2:
		return !_attd2s1.empty();
	case 3:
		return !_attd3s1.empty();
	case 4:
		return !_attd4s1.empty();
	case 5:
		return !_attd5s1.empty();
	case 6:
		return !_attd6s1.empty();
	case 7:
		return !_attd7s1.empty();
	}
	return false;
}
}
