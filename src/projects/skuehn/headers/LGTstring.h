#ifndef LGTSTRING_H
#define LGTSTRING_H

#include <fstream>

/** Encapsulate a string as starting position and length, how it acutally looks like is dependent on the model */
typedef struct string_t
{
  int pos;
  int leng;
} string_t;

bool operator==(string_t str1,string_t str2);
bool operator!=(string_t str1,string_t str2);
std::ostream& operator<<(std::ostream& os,string_t str);

#endif // LGTSTRING_H
