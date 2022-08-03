#ifndef MISC_H
#define MISC_H

#include <string>

/** Return true if the file with the given name exists */

bool file_exists(const std::string filename);
bool file_exists(const char* filename);


int gcd(int a,int b);
int lcm(int a,int b);


#endif // MISC_H
