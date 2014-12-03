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


#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <vector>
#include <map>
#include "CIndividu.h"
#include "CCohort.h"
#include "CStats.h"

int FCorrectPoisson(long double& CorrectionPoisson, long double& CorrectionBinomial);
CIndividu FPoissonSimul(unsigned int& IndKey,long double const& CorrectionPoisson, long double const& CorrectionBinomial, unsigned int const& Year);
CIndividu FNeutralSimul(unsigned int& NIndKey, unsigned int const& NYear,double const& RealizedJuvenileSurv,double const& RealizedAdultSurv,vector<long double> const& RealizedARS, unsigned int const& PoissonMaxLRS);
int FDemographer(std::vector<CCohort> const& PopulationS,long double& H, long double& P);
std::vector<std::vector<long double> > FEcologist(double& RealizedAdultSurv,double& RealizedJuvenileSurv, std::vector<CCohort>& PopulationS);
std::vector<long double> FNeutralEcologist(double& NeutralJuvenileSurv,double& NeutralAdultSurv,std::vector<CCohort>& PopulationN);
int FAverageLRS(std::vector<long double>& MeanNeutralLRS,std::vector<long double>& NeutralLRS, long double& M, long double& V);
long double FMean(std::vector<long double> const& LRSDistri);
long double FVariance(std::vector<long double> const& LRSDistri, long double const& Mean);
long double FEntropy(std::map<unsigned int, std::map<unsigned int, long double> >& Transitions, std::map<unsigned int, long double>& Collector);
long double FPersistence(std::map<unsigned int, std::map<unsigned int, long double> >& Transitions, std::map<unsigned int, long double>& Collector);
int FWritingLRSFiles(std::vector<std::vector<long double> >& AllMeanNeutralLRS,std::vector<std::vector<long double> >& AllRealizedLRS,unsigned int const& AbsoluteMaxLRS);
int FWritingStatsFiles(CStats& HPMVcontainer);
int FWritingMMFiles(std::vector<std::vector<CCohort> >& AllPopulationS);
int FWritingRunInfo();
#endif // MAIN_H_INCLUDED
