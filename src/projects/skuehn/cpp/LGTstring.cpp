#include "LGTstring.h"

using namespace std;


bool operator==(string_t str1,string_t str2)
{
  return (str1.pos==str2.pos)&&(str1.leng==str2.leng);
}

bool operator!=(string_t str1,string_t str2)
{
  return (str1.pos!=str2.pos)||(str1.leng!=str2.leng);
}

#include <iostream>
ostream& operator<<(ostream& os,string_t str)
{
  os<< "(" << str.pos << "," << str.leng << ")";
  return os;
}


