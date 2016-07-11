#ifndef ANALYSE_H_INCLUDED
#define ANALYSE_H_INCLUDED

#include "HaplotypeSet.h"
//#include <unordered_set>
#include<map>

using namespace std;




class Analyse
    {

        public:

        bool BgZip;
        String OutputFileName;
		vector<double> AlleleFreq;
		vector<int> AlleleCount;
        vector<bool> major;

        float m3vcfVersion;

//
//		vector<vector<double> > HapR2;


        String Hapr2PosFile;
        int LdWindow,LdWindowMin,LdWindowBp,LdWindowBpMin;
        double MinHapR2;


            Analyse(bool bgzip,String &filename)
            {
                BgZip=bgzip;
                OutputFileName=filename;

            }


            void writem3vcfFile         (HaplotypeSet &ThisData);
            void PrintFreqOutput        (HaplotypeSet &ThisData);
            void CalculateFreq          (HaplotypeSet &ThisData);
            void PrintCountOutput       (HaplotypeSet &ThisData);
            void CalculatCount          (HaplotypeSet &ThisData);
            void PrintHapR2Output       (vector<vector<double> > &HapR2, HaplotypeSet &ThisData, int Block1, int Block2);
            void CalculateHapCrossBlock (vector<vector<double> > &HapR2, HaplotypeSet &ThisData, int Block1, int Block2);
            void CalculateHapSelfBlock  (vector<vector<double> > &HapR2, HaplotypeSet &ThisData, int Block1);
            void CalculatHapR2          (HaplotypeSet &ThisData);

            double CalculateDPrime(double &ThisD, double &iFreq, double &jFreq);



            void    InitializeLDParameters(String &hapr2PosFile,int &ldWindow,int &ldWindowMin,int &ldWindowBp,
                                           int &ldWindowBpMin,double &minHapR2)

            {

//                String &hapr2PosFile,int &ldWindow,int &ldWindowMin,int &ldWindowBp,
//                                               int &ldWindowBpMin,float &minHapR2


                Hapr2PosFile=hapr2PosFile.c_str();
                LdWindow=ldWindow;
                LdWindowMin=ldWindowMin;
                LdWindowBp=ldWindowBp;
                LdWindowBpMin=ldWindowBpMin;
                MinHapR2=minHapR2;


cout<<LdWindow<<endl;
            }
//MinHapR2




};




#endif // ANALYSE_H_INCLUDED
