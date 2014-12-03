//    DynFixSim simulate animal populations with various levels of fixed and dynamic heterogeneity in life history
//    and extract statistics for further analysis.
//
//	  Copyright (C) 2014  Timothée Bonnet - timothee.bonnet@ieu.uzh.ch
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

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
