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
