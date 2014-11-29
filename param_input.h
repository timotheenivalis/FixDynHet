#ifndef PARAM_INPUT_H_INCLUDED
#define PARAM_INPUT_H_INCLUDED

#include <vector>

int cmp_nocase(const std::string& s, const std::string& s2);
void rtrim(std::string *s);
int evaluateBool(bool &boolean, std::string buf);
int seeks_settings_file_name(const std::string cmdlinefilename,std::string& settingsfilename);
int read_settings_file(const std::string filename);


extern unsigned int PoissonSimulNumber;
extern unsigned int NeutralSimulNumber;

extern double AdultSurvival;
extern double JuvenileSurvival;
extern std::vector<long double> IndVarReproS;
extern std::vector<long double> IndVarSurvS;
extern std::vector<long double> CorReproSurvS;
extern std::vector<long double> MarkovS;
extern long double VarianceCutOff;

extern std::vector<int> MeanReproS;
extern unsigned int AgeMax;
extern unsigned int StudyLength;
extern unsigned int IndInCohort;

extern unsigned long int _ptSamplingSeed;
extern bool cinGetOnError;
extern bool pauseGP;
extern bool DisplayProgression;
#endif // PARAM_INPUT_H_INCLUDED
