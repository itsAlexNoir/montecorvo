/*!
  @file params.hpp
*/

#ifndef ____PARAMS__
#define ____PARAMS__

#include <iostream>
#include "constants.hpp"
#include "boost/variant.hpp"

typedef boost::variant<std::string, int, double> input_type;
typedef std::map<std::string, input_type> inputMap;

class input{
private:
  std::ifstream infile;

public:
  input(std::string filename);
  /*!
    @function ~input
    \brief Class destructor
  */
  ~input(){infile.close();}

  /// Class methods
  inputMap read_input(std::string filename);
  int read_integer(std::string mykey, int default_value = 1);
  double read_double(std::string mykey, double default_value = 1.0);
  std::string read_string(std::string mykey, std::string default_value = "None");
  bool read_boolean(std::string mykey, bool default_value = false);
};


#endif // ____PARAMS__
