#ifndef __TOKENIZE_H__
#define __TOKENIZE_H__


#include <vector>
#include <string>
#include <stdlib.h>
using namespace std;

string MyTokenize(const char *input, const char *delimiter, int Number);

void MyTokenize(vector<string> &result, const char *input, const char *delimiter, int Number);

string FindTokenWithPrefix(const char *input,const char *delimiter, string CheckPrefix);



#endif
