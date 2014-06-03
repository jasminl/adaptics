#pragma once

#include <list>
#include <string>
 
/**
 *	Manages simulation directories
 */
class CDirManager
{
public:
	//typedefs, etc.
	static const int DEFAULT_START;				//!<Indicates to start numbering output folders from 1
	enum{USE_MANAGER = true, NO_USE_MANAGER = false};	//!<Whether to use the directory manager or not

public:
	/**
	 *	\param output_folder The root folder where all results are stored
	 *	\param startFolder  The folder number from which to start numbering.
	 */
	CDirManager(const std::string& output_folder, int startFolder);

	/**
	 *	Destructor
	 */
	~CDirManager();

	/**
	 *	Indicates whether to use the directory manager or not
	 *	\return Boolean value indicating whether to use or not the manager
	 */
	bool use() const;

	/** 
	 *	Obtains the name of the current directory
	 *	\return The name of the current directory
	 */
	std::string currentDir();

	/**
	 *	Gets the working directory
	 *	\return The name of the working directory
	 */
	const char* wkdir() const
	{
		return m_cwdir.c_str();
	}

	/**
	 *	Creates new directory, based on existing list
	 */
	std::string createNewDir();

	/**
	 *	Indicates whether a certain directory exists
	 *	\param dname Directory of which to check the existence.
	 *	\return Boolean value indicating the existence of the directory.
	 */
	bool exist(const char* dname) const;
	
	/**
	 *	Moves a given xml file to destination directory
	 *	\param Fname Name of the file to move.
	 *	\param dname Name of the directory where to move it.
	 */
	void moveFile(std::string Fname, const char* dname);

	/**
	 * Loads a file and returns a list of strings
	 */
	std::list<std::string> load_file(const std::string& file_path);

	/**
	 * Creates a new current parameter directory
	 */
	std::string mk_par_dir();

	/**
	 * Creates a new current stimulus directory
	 */
	std::string mk_stim_dir(const std::string& stim_file);

	/**
	 * Retrieve the current parameter directory
	 */
	std::string par_dir() const;

	/**
	 * Retrieve the current stimulus directory
	 */
	std::string stim_dir() const;

	/**
	 * Retrieve the root directory
	 */
	std::string root_dir() const;

protected: //Data member
	std::list<int> m_existingDir;		//!< List of existing directories
	std::string   m_cwdir;				//!< Current working directory
	unsigned m_bufferLength;			//!< Maximum length of directory names
	bool     m_use;						//!< Indicates whether to use the CDirManager or not
	unsigned m_currentDirNumber;		//!< The number specifying the next directory to be created
	std::string _cur_p_dir;				//!< Current parameter directory
	std::string _cur_s_dir;				//!< Current stimulus directory
};

inline
bool CDirManager::use() const
{
	return m_use;
}

inline
std::string CDirManager::par_dir() const
{
	return _cur_p_dir;
}

inline
std::string CDirManager::stim_dir() const
{
	return _cur_s_dir;
}

inline
std::string CDirManager::root_dir() const
{
	return m_cwdir;
}
