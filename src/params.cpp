/*! 
  @file params.cpp
  @brief Contains functions for reading input files.
  
  File containing the needed functions for reading input files.
  The main function read a file, and create a map with pairs of
  the parameter and his value.
*/

#include "params.hpp"

/*!
  @function input
  @brief Class constructor 
  
  @params[in] string filename
  @params[out] map readmap
*/

input::input(std::string filename)
{
  infile.open(filename,std::ifstream::in);
  if(!infile.is_open())
    cout << "Unable to open input file";
  
}

/*!
  @function read_input
  @brief Function for reading input files. 
  A map is created with the loaded parameters.
  
  @params[in] string filename
  @params[out] map readmap
*/

inputMap input::read_input(std::string filename)
{
  std::string key, line, eqstr;
  std::string strvalue;
  inputMap readmap;
  
  // Read the file
  while(std::getline(infile,line))
    {
      if( (line.length() != 0) && (line[0] != '#') && line.find('='))
	{
	  std::istringstream fileline(line);
	  while(fileline >> key >> eqstr >> strvalue )
	    readmap[key] = strvalue;
	}
      
    }
  return readmap;	
}
/*!
  @function read_integer
  @brief Read integer from input file
  
  @params[in] std::string mykey Wanted parameter to read from
  @params[in] int default_value Default value if wanted value not found
  @params[out] int Integer value 
*/
int input::read_integer(std::string mykey, int default_value)
{
  std::string key, line, eqstr;
  int         intvalue;
  
  while(std::getline(infile,line))
    if( (line.length() != 0) && (line[0] != '#') && line.find('=') )
      {
	std::istringstream fileline(line);
	fileline >> key >> eqstr >> intvalue;
	if(key==mykey)
	  {
	    infile.clear();
	    infile.seekg(0, ios::beg);
	    return intvalue;
	  }
      }
  infile.clear();
  infile.seekg(0, ios::beg);
  return default_value;
}
/*!
  @function read_double
  @brief Read double from input file
  
  @params[in] std::string mykey Wanted parameter to read from
  @params[in] double default_value Default value if wanted value not found
  @params[out] int Double value 
*/
double input::read_double(std::string mykey, double default_value)
{
  std::string key, line, eqstr;
  double      dvalue;
  
  while(std::getline(infile,line))
    if( (line.length() != 0) && (line[0] != '#') && line.find('=') )
      {
	std::istringstream fileline(line);
	fileline >> key >> eqstr >> dvalue;
	if(key==mykey)
	  {
	    infile.clear();
	    infile.seekg(0, ios::beg);
	    return dvalue;
	  }
      }
	
  infile.clear();
  infile.seekg(0, ios::beg);
  return default_value;
}
/*!
  @function read_string
  @brief Read an string from input file
  
  @params[in] std::string mykey Wanted parameter to read from
  @params[in] std::string default_value Default value if wanted value not found
  @params[out] std::string Integer value 
*/
std::string input::read_string(std::string mykey, std::string default_value)
{
  std::string key, line, eqstr, strvalue;
  
  while(std::getline(infile,line))
    if( (line.length() != 0) && (line[0] != '#') && line.find('=') )
      {
	std::istringstream fileline(line);
	fileline >> key >> eqstr >> strvalue;
	if(key==mykey)
	  {
	    infile.clear();
	    infile.seekg(0, ios::beg);
	    return strvalue;
	  }
      }
  infile.clear();
  infile.seekg(0, ios::beg);
  return default_value;
}

/*!
  @function read_boolean
  @brief Read an string from input file
  
  @params[in] std::string mykey Wanted parameter to read from
  @params[in] bool default_value Default value if wanted value not found
  @params[out] std::string Integer value 
*/
bool input::read_boolean(std::string mykey, bool default_value)
{
  std::string key, line, eqstr;
  bool        boolvalue;
  
  while(std::getline(infile,line))
    if( (line.length() != 0) && (line[0] != '#') && line.find('=') )
      {
	std::istringstream fileline(line);
	fileline >> key >> eqstr >> boolvalue;
	if(key==mykey)
	  {
	    infile.clear();
	    infile.seekg(0, ios::beg);
	    return boolvalue;
	  }
      }
  
  infile.clear();
  infile.seekg(0, ios::beg);
  return default_value;
}
