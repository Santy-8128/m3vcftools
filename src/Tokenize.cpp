#include "Tokenize.h"
using namespace std;

string MyTokenize(const char *input, const char *delimiter, int Number)
{
    string result;
    size_t wordCount = 1;
    result.clear();
    std::string *word = &result;


    while (*input)
    {
        if (*input==*delimiter)
        {
            // we got a delimeter, and since an empty word following
            // a delimeter still counts as a word, we allocate it here
            wordCount++;

            if((int)wordCount>Number)
                return result;

            result.clear();
            word = &result;
        }
        else
        {
            word->push_back(*input);
        }
        input++;
    }
    return result;

}

void MyTokenize(vector<string> &result, const char *input, const char *delimiter, int Number)
{

    size_t wordCount = 1;
    result[0].clear();
    std::string *word = &result[0];


    while (*input)
    {
        if (*input==*delimiter)
        {
            // we got a delimeter, and since an empty word following
            // a delimeter still counts as a word, we allocate it here
            wordCount++;

            if((int)wordCount>Number)
                return;

            result[wordCount-1].clear();
            word = &result[wordCount-1];
        }
        else
        {
            word->push_back(*input);
        }
        input++;
    }

}

string FindTokenWithPrefix(const char *input,const char *delimiter, string CheckPrefix)
{

    std::string word = "";
    int Size = (int)CheckPrefix.size();
    while (*input)
    {
        if (*input==*delimiter)
        {
            int Index=0;
            ++input;
            while(*input)
            {
                word=word + (*input);
                if(Index<Size && *input!=CheckPrefix[Index++])
                    break;

                ++input;
                if(*input==*delimiter || *input=='\0')
                    return word;
            }
        }
        input++;
        word="";
    }
    return word;

}

