
/** 
	Exeption for parameter being in wrong range.
*/
class ex_wrong_range: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Parameter is in wrong range.";
  }
} x_wrong_range;