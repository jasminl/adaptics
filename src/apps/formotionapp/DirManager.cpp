#include "tinyxml/tinyxml.h"
#include "visio/Layer.h"
#include "DirManager.h"
#include <time.h>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <stdexcept>

#pragma GCC diagnostic ignored "-Wunused-result"	//Since some of the IO functions generate that warning

using namespace std;

const int CDirManager::DEFAULT_START = -1;				//!<Indicates to start numbering output folders from 1

CDirManager::CDirManager(const string& output_folder, int startFolder)
: m_bufferLength(200), m_use(true)
{
	if(startFolder == DEFAULT_START)
		m_currentDirNumber = 1;
	else
		m_currentDirNumber = startFolder;

	//setCDir();
	m_cwdir = output_folder;
}

CDirManager::~CDirManager()
{

}

string CDirManager::currentDir()
{
	if(!use()) return string();

	char* buffer = new char[m_bufferLength];
	memset(buffer,0,sizeof(char)*m_bufferLength);

	stringstream str;
	str<<m_currentDirNumber;
	//_itoa_s(m_currentDirNumber,buffer,m_bufferLength,10);

	string dname = str.str();//string(buffer);
	delete[] buffer;

	return dname;
}

bool CDirManager::exist(const char* dname) const
{
	if(!use()) return false;

	char* buffer = new char[m_bufferLength];
	memset(buffer,0,sizeof(char)*m_bufferLength);

	getcwd(buffer, m_bufferLength);	//Get current directory to return to it later

	bool answer = false;

	if(0 == chdir(dname))
	{
		chdir(buffer);
		answer = true;
	}

	delete[] buffer;

	return answer;
}

string CDirManager::createNewDir()
{

	string p_dir;
	while(true)
	{
		stringstream str;
		str<<m_currentDirNumber;
		if(0 == chdir((m_cwdir + str.str()).c_str()))
		{
			//Indicates that this directory already exist, therefore we don't want it
			m_currentDirNumber++;		//Increase number
			chdir(m_cwdir.c_str());		//Go back to working directory
		}
		else
		{
			p_dir = str.str();
			break;
		}
	}

	if(-1 == mkdir((m_cwdir + "/" + p_dir).c_str(), S_IRWXU|S_IRGRP|S_IXGRP))
		throw runtime_error("Cannot create directory " + m_cwdir + "/" + p_dir);

	cout<<"Created directory "<<p_dir<<endl;
	m_existingDir.push_back(m_currentDirNumber);
	return p_dir;

}

void CDirManager::moveFile(string Fname, const char* dname)
{
	if(true == exist(dname))
	{
		TiXmlDocument doc(Fname.c_str());
		if(!doc.LoadFile())
		{
			//XML file can't be loaded
			cout<<"Cannot open XML file: "<<Fname<<", in 'moveFile'"<<endl;
			return;
		}
		else
		{
			int base_length = Fname.size() - Fname.rfind("/");
			string base_name = Fname.substr(Fname.rfind("/"), base_length);
			string output_file = string(dname) + "/" + base_name;

			if(false == doc.SaveFile(output_file.c_str()))
				cout<<"Couldn't create copy of "<<Fname<<", in 'moveFile'.\n";
		}
	}
	else
	{
		//Directory does not exist
		cout<<"Invalid directory: "<<dname<<", in 'moveFile'\n.";
		return;
	}
}

list<string> CDirManager::load_file(const string& file_path)
{
	ifstream file(file_path);
	if(file.bad())
		throw runtime_error("CDirManager::load_file: file does not exist: " + file_path);

	list<string> out;
	if(string::npos != file_path.find(".xml",file_path.size()-4))
		out.push_back(file_path); //The file is a single xml file: insert into m_pfile list
	else if(string::npos != file_path.find(".txt",file_path.size()-4))
	{	//This is a bunch of xml files
		string str;
		ifstream txt_file(file_path);
		while(txt_file >> str)
		{	//Read each line as an element
			out.push_back(str);
			str.clear();
		}
	}

	return out;
}

string CDirManager::mk_par_dir()
{
	_cur_p_dir = createNewDir();
	return _cur_p_dir;
}

string CDirManager::mk_stim_dir(const string& stim_file)
{
	int base_length = stim_file.size() - stim_file.rfind("/");
	string base_name = stim_file.substr(stim_file.rfind("/"), base_length);
	_cur_s_dir = m_cwdir + "/" + _cur_p_dir + "/" + base_name;

	if(-1 == mkdir(_cur_s_dir.c_str(), S_IRWXU|S_IRGRP|S_IXGRP))
		throw runtime_error("CDirManager::mk_stim_dir: cannot create directory " + _cur_s_dir);

	return _cur_s_dir;
}
