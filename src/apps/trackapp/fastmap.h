#pragma once

#include <vector>
#include <string.h>

/**
	std::map where key is not const
*/
template<typename T, typename Q> class fastmap
{
protected:

	std::vector<T> m_key;
	std::vector< Q* > m_data;

	std::vector<unsigned int> m_size;			/*<<Number of entries in the array */

	std::vector<unsigned int> m_maxSize;			/*<<Actual space occupied */


public:

	/**
		Constructor
	*/
	fastmap<T,Q>(std::vector<T> key):m_key(key){}

	/**
		Default constructor
	*/
	fastmap<T,Q>(){}

	~fastmap()
	{
		clear();
	}

	/**
		Add space for one element
	*/
	void add(T key, unsigned int size)
	{
		m_key.push_back(key);
		m_data.push_back(new Q[size]);
		m_size.push_back(size);

		m_maxSize.push_back(size);
	}

	/**
		Operator[]
	*/
	std::pair<Q*,unsigned int> operator[](T index)
	{
		typename std::vector<Q*>::iterator q = m_data.begin();
		std::vector<unsigned int>::iterator r = m_size.begin();

		for(auto p =m_key.begin(); p != m_key.end(); p++, q++, r++)
		{
			if(*p == index)
			{
				return std::pair<Q*, unsigned int>(*q,*r);
			}

		}
		return std::pair<Q*,unsigned int>(nullptr, 0);	//Index not found
	}

	/**
		Operator () accesses by address
	*/
	std::pair<Q*, unsigned int > operator()(unsigned int index)
	{
		return std::pair<Q*, unsigned int >(m_data[index],m_size[index]);
	}

	/**
		Gives max size for a given scale
	*/
	unsigned& maxSize(T index)
	{
		std::vector<unsigned int>::iterator q = m_maxSize.begin();
		for(auto p = m_key.begin() ; p != m_key.end(); p++,q++)
		{
			if(index == *p)
			{
				return *q;
			}
		}
	}

	/**
		Reallocates space and/or update index
	*/
	Q* realloc(T index, unsigned int minSize)
	{
		typename std::vector<Q*>::iterator p = m_data.begin();
		for(auto q = m_key.begin() ; q != m_key.end() ; q++,p++)
		{
			if(index == *q)
			{
				if(maxSize(index) < minSize)
				{
					delete *p;

					//Reallocate only if is smaller
					*p = new Q[minSize];

					//Set to 0
					Q* arr = *p;
					memset(arr,0,sizeof(Q)*minSize);

					maxSize(index) = minSize;
					actualSize(index) = minSize;
				}

				return *p;	//Always exit function with current array, whether or not it was reassigned
			}
		}

		//Not found so add it
		m_key.push_back(index);
		m_size.push_back(minSize);
		m_maxSize.push_back(minSize);
		m_data.push_back(new Q[minSize]);

		return m_data[m_data.size()-1];

	}

	/**
		Change scale index
	*/
	void changeIndex(unsigned int index,T value)
	{
		if(index>=m_key.size())
		{
			m_key.push_back(value);

			m_size.push_back(0);
			m_maxSize.push_back(0);
			m_data.push_back(nullptr);
		}

		m_key[index] = value;
	}

	/**
		Nb of elements
	*/
	unsigned int size() const
	{
		return (unsigned int)m_key.size();
	}

	/**
		Empty
	*/
	void clear()
	{
		m_key.empty();
		m_size.empty();
		for(auto p = m_data.begin() ; p != m_data.end() ; p++)
		{
			delete *p;
		}
	}

	std::pair<Q*,unsigned int> centerScale()
	{
		unsigned int index = static_cast<unsigned int>(m_key.size()/ 2);
		return std::pair<Q*,unsigned int> (m_data[index],m_size[index]);

	}

	unsigned int& actualSize(T index)
	{
		std::vector<unsigned int>::iterator q = m_size.begin();
		for(auto p = m_key.begin() ; p != m_key.end() ; p++, q++)
		{
			if(*p == index)
			{
				return *q;
			}
		}
	}

};
