//    DynFixSim simulate animal populations with various levels of fixed and dynamic heterogeneity in life history
//    and extract statistics for further analysis.
//
//	  Copyright (C) 2014-2015  Timothée Bonnet - timothee.bonnet@ieu.uzh.ch
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

//instructions preprocesseur
#include <iostream>
#include <sstream> //stringstream
#include <fstream> // ofstream
#include <iomanip>//setprecision
#include <vector> // for vector
#include <map> //for map
#include <cstdlib> //for folder creation
//#include <direct.h> //for mkdir
#include <dirent.h> //alternative for mkdir which should be cross platform
#include <ctime>//clock


//#include <sys/stat.h>
//#include <sys/types.h>


using namespace std;

#include "main.h"
#include "MersenneTwister.h"
#include "param_input.h"
#include "CIndividu.h"
#include "CCohort.h"
#include "GaussianGenerator.h"
#include "CStats.h"

unsigned int PoissonSimulNumber=10;
unsigned int NeutralSimulNumber=10;
double AdultSurvival=0.2;
double JuvenileSurvival=0.4;

vector<long double> IndVarSurvS;
vector<long double> IndVarReproS;
vector<long double> CorReproSurvS;
vector<long double> MarkovS;

long double IndVarRepro=0.;
long double IndVarSurv=0.;
long double CorReproSurv=0.;
long double Markov=0.;
long double VarianceCutOff=3; // how much of the variance in qualities should we discard?

vector<int> MeanReproS;
int MeanRepro=3;
unsigned int AgeMax=5;
unsigned int StudyLength=10;
unsigned int IndInCohort=100;


bool cinGetOnError;
bool pauseGP;
bool DisplayProgression=false;



MTRand alea;
unsigned long int _ptSamplingSeed=67144630;

string cmdlinefilename="cmdlineArguments.txt";
string settingsfilename="DynFixSimParam.txt";//input file

int main(int argc, char *argv[])
{
    clock_t start,end;
    double temps_ecoule(0.0);
    start=clock();
    if (argc>1)
     {
        // to give inline the name of the file in which command line is written

        string buf(argv[1]);
        string::size_type pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
        string var=buf.substr(0,pos).c_str();
        if(cmp_nocase(var,"CmdlineFileName")==0) cmdlinefilename=buf.substr(pos+1);
        ofstream cmdline(cmdlinefilename.c_str(),ios::out);
        for (int it=1;it<=argc;it++) cmdline<<argv[it]<<endl;
        cmdline<<endl;
        cmdline.close();
        // seeks optional SettingsFile in cmdline
        seeks_settings_file_name(cmdlinefilename,settingsfilename);
     }
    read_settings_file(settingsfilename); //cf migraine.cpp... READS and SETS...
    if (argc>1) read_settings_file(cmdlinefilename);
    //remove(cmdlinefilename.c_str());// on le supprime si on prevoit de faire plein de trucs dans le meme dossier

    alea.seed(_ptSamplingSeed);

    for (unsigned int RV(0);RV<IndVarReproS.size();RV++)
    for (unsigned int SV(0);SV<IndVarSurvS.size();SV++)
    for (unsigned int CRS(0);CRS<CorReproSurvS.size();CRS++)
    for (unsigned int MR(0);MR<MeanReproS.size();MR++)
    for (unsigned int MaK(0);MaK<MarkovS.size();MaK++)
        {
            IndVarRepro=IndVarReproS[RV];
            IndVarSurv=IndVarSurvS[SV];
            CorReproSurv=CorReproSurvS[CRS];
            MeanRepro=MeanReproS[MR];
            Markov=MarkovS[MaK];

            vector<vector<long double > > AllMeanNeutralLRS;
            vector<vector<long double> > AllRealizedLRS;
            unsigned int AbsoluteMaxLRS(0);
            unsigned int PoissonMaxLRS(0);
            long double H(0.);
            long double P(0.);
            long double M(0.);
            long double V(0.);
            vector<long double> HlistSelection;
            vector<long double> PlistSelection;
            vector<long double> MlistSelection;
            vector<long double> VlistSelection;
            vector<long double> ASlistSelection;
            vector<long double> JSlistSelection;

            vector<vector<long double> > HlistNeutral;
            vector<vector<long double> > PlistNeutral;
            vector<vector<long double> > MlistNeutral;
            vector<vector<long double> > VlistNeutral;
            vector<vector<long double> > ASlistNeutral;
            vector<vector<long double> > JSlistNeutral;

            vector<vector<CCohort> > AllPopulationS;

            long double CorrectionPoisson(1.);
            long double CorrectionBinomial(0.);
            FCorrectPoisson(CorrectionPoisson,CorrectionBinomial);// aims at correcting the mean reproduction for the asymmetrical variance introduced
            for (unsigned int nbSPoisson(0);nbSPoisson<PoissonSimulNumber;nbSPoisson++)
                {
                    vector<CCohort> PopulationS;
                    unsigned int IndKey(0);
                    for (unsigned int Year(0);Year<StudyLength;Year++)
                        {
                            CCohort cohort;
                            cohort.CohortBirth=Year;
                            for (unsigned int Ind(0);Ind<IndInCohort;Ind++)
                                {
                                    CIndividu individual=FPoissonSimul(IndKey,CorrectionPoisson,CorrectionBinomial,Year);
                                    cohort.Individus.push_back(individual);
                                }
                            PopulationS.push_back(cohort);
                        }
                    AllPopulationS.push_back(PopulationS);

                    double RealizedAdultSurv(0.);
                    double RealizedJuvenileSurv(0.);

                    FDemographer(PopulationS,H,P);
                    HlistSelection.push_back(H);
                    PlistSelection.push_back(P);

                    vector<vector<long double> > Contain=FEcologist(RealizedAdultSurv,RealizedJuvenileSurv,PopulationS);//estimates realized survival and Annual Reproductive Success from Selective Simulations
                    ASlistSelection.push_back(RealizedAdultSurv);
                    JSlistSelection.push_back(RealizedJuvenileSurv);
                    vector<long double> RealizedARS;
                    RealizedARS.resize(Contain[0].size());
                    RealizedARS=Contain[0];
                    vector<long double> RealizedLRS;
                    RealizedLRS.resize(Contain[1].size());
                    RealizedLRS=Contain[1];
                    AllRealizedLRS.push_back(RealizedLRS);

                    M=FMean(RealizedLRS);
                    V=FVariance(RealizedLRS,M);
                    MlistSelection.push_back(M);
                    VlistSelection.push_back(V);

                    if(PoissonMaxLRS<Contain[1].size())
                        {
                            PoissonMaxLRS=Contain[1].size();//we need it for latter writting, and we will update the value during neutral simulations
                        }
                    unsigned int NStudyLength=StudyLength;
                    unsigned int NIndInCohort=IndInCohort;

                    vector<long double> MeanNeutralLRS;//alternative to the storage
                    vector<long double> RunNeutralH;//entropy
                    vector<long double> RunNeutralP;//persistence
                    vector<long double> RunNeutralM;//mean LRS
                    vector<long double> RunNeutralV;//variance in LRS
                    vector<long double> RunNeutralAS;//adult survival
                    vector<long double> RunNeutralJS;//juvenile survival

                    for (unsigned int NSim(0);NSim<NeutralSimulNumber;NSim++)
                        {
                        if (DisplayProgression==true)  cerr<<"\r nbSPoisson "<<nbSPoisson<<" nbNeutral "<<NSim<<" starting"<<flush;
                            vector<CCohort> PopulationN;
                            unsigned int NIndKey(0);
                            for (unsigned int NYear(0);NYear<NStudyLength;NYear++)//neutral simulations
                                {
                                    CCohort Ncohort;
                                    Ncohort.CohortBirth=NYear;
                                    for (unsigned int NInd(0);NInd<NIndInCohort;NInd++)
                                        {
                                            CIndividu Nindividual=FNeutralSimul(NIndKey,NYear,RealizedJuvenileSurv,RealizedAdultSurv,RealizedARS,PoissonMaxLRS);
                                            Ncohort.Individus.push_back(Nindividual);
                                        }
                                    PopulationN.push_back(Ncohort);
                                }

                            FDemographer(PopulationN,H,P);
                            RunNeutralH.push_back(H);
                            RunNeutralP.push_back(P);

                            double NeutralAdultSurv(0.);
                            double NeutralJuvenileSurv(0.);

                            vector<long double> NeutralLRS=FNeutralEcologist(NeutralJuvenileSurv,NeutralAdultSurv,PopulationN);
                            RunNeutralAS.push_back(NeutralAdultSurv);
                            RunNeutralJS.push_back(NeutralJuvenileSurv);
                            FAverageLRS(MeanNeutralLRS,NeutralLRS,M,V);
                            RunNeutralM.push_back(M);
                            RunNeutralV.push_back(V);

                            if(AbsoluteMaxLRS<(NeutralLRS.size())) AbsoluteMaxLRS=(NeutralLRS.size());//update the max size for latter file writting

                        }
                    if(AbsoluteMaxLRS<PoissonMaxLRS) AbsoluteMaxLRS=PoissonMaxLRS;

                    AllMeanNeutralLRS.push_back(MeanNeutralLRS);
                    HlistNeutral.push_back(RunNeutralH);
                    PlistNeutral.push_back(RunNeutralP);
                    MlistNeutral.push_back(RunNeutralM);
                    VlistNeutral.push_back(RunNeutralV);
                    ASlistNeutral.push_back(RunNeutralAS);
                    JSlistNeutral.push_back(RunNeutralJS);

                }//end for (unsigned int nbSPoisson(1);nbSPoisson<=PoissonSimulNuber;nbSPoisson++)


            CStats HPMVcontainer;
            HPMVcontainer.NH=HlistNeutral;
            HPMVcontainer.NP=PlistNeutral;
            HPMVcontainer.NM=MlistNeutral;
            HPMVcontainer.NV=VlistNeutral;
            HPMVcontainer.NAS=ASlistNeutral;
            HPMVcontainer.NJS=JSlistNeutral;
            HPMVcontainer.SH=HlistSelection;
            HPMVcontainer.SP=PlistSelection;
            HPMVcontainer.SM=MlistSelection;
            HPMVcontainer.SV=VlistSelection;
            HPMVcontainer.SAS=ASlistSelection;
            HPMVcontainer.SJS=JSlistSelection;

            stringstream stst;
            stst << IndVarRepro; //transforms integer in strstr
            stst<<IndVarSurv;
            stst<<CorReproSurv;
            stst<<MeanRepro;
            stst<<Markov;

            string folderName("OutputDFS_");
            folderName+=stst.str();

            cout<<endl;
            string command("mkdir ");
            system((command+folderName).c_str());

            chdir(folderName.c_str());

            FWritingLRSFiles(AllMeanNeutralLRS,AllRealizedLRS,AbsoluteMaxLRS);

            FWritingStatsFiles(HPMVcontainer);


            FWritingMMFiles(AllPopulationS);
            FWritingRunInfo();
            string stepback("\..");

            chdir(stepback.c_str());

            if (DisplayProgression==true)
                {
                    end=clock();
                    temps_ecoule=double(end-start)/1000;
                    cout<<endl<<"time spent="<<temps_ecoule<<"s"<<endl;
                }
        }//end     for (unsigned int RV(0);RV<IndVarReproS.size();RV++)
    return 0;
}//end int main(int argc,char *argv[])

/***********************************************************************************************************************************************************************/
//**************************************************************/SECONDARY FUNCTIONS/**********************************************************************************//
/***********************************************************************************************************************************************************************/
int FCorrectPoisson(long double& CorrectionPoisson, long double& CorrectionBinomial)
{
    long double CorrectlogscaleLambda(0.);
    long double sumRepro(0.);
    long double nbcorrection(10000);

    if(IndVarRepro!=0.)
        {
            for (unsigned int i(0);i<nbcorrection;i++)
                {

                    CorrectlogscaleLambda=exp(log(MeanRepro)+sqrt(IndVarRepro)*normsinv(alea()));
                    long double L=exp(-double(CorrectlogscaleLambda));//generates Poisson distribution
                    long unsigned int k=0;
                    double p=1;
                    do{
                    k++;
                    p*=alea();
                    }while(p>L);
                    k--;
                    sumRepro+=k;
                }
            sumRepro/=nbcorrection;
            CorrectionPoisson=sumRepro/MeanRepro;
        }

    long double sumSurv(0.);
    if(IndVarSurv!=0.)
        {
            long double nbstep(0.);
            for (unsigned int i(0);i<nbcorrection;i++)
                {

                    long double Alive(1.);
                    long double intercept=log(AdultSurvival/(1-AdultSurvival));
                    long double logitvalue=intercept+sqrt(IndVarSurv)*normsinv(alea());
                    long double Phi=exp(logitvalue)/(1+exp(logitvalue));
                    unsigned int age(1);
                    while ((Alive==1.) & (age<=AgeMax))
                        {
                            long double fatum=alea();
                            if(fatum>Phi)
                                {
                                    Alive=0.;
                                }
                            sumSurv+=Alive;
                            age++;
                            nbstep++;
                        }
                }
            sumSurv/=nbstep;
            CorrectionBinomial=log(sumSurv/(1-sumSurv))-log(AdultSurvival/(1-AdultSurvival));
        }
    return 0;
}//end FCorrectPoisson()
/***************************************************************************************************************************/

CIndividu FPoissonSimul(unsigned int& IndKey,long double const& CorrectionPoisson, long double const& CorrectionBinomial, unsigned int const& Year)
{
    bool Alive(1);
    unsigned int age(1);
    int LHS=1;// 1 means Juvenile, 2 means Adult
    unsigned int repro(0);
    long double pastRepro(MeanRepro/CorrectionPoisson); // starts with mean repro so that when there is markovian effect the first reproduction is completely random.
    CIndividu Individu;
    CLifeStage LHind;
    Individu.IndividuKey=IndKey;
    long double ranRepro(0.);
    long double ranSurv(0.);
    bool OK(0);
    while(OK==0)// these two while loops ensure that we do not get too extreme qualities, most importantly avoiding huge ARS which are memory consuming. The variance should be minored
    {
        ranRepro=normsinv(alea());
        if (abs(ranRepro)<VarianceCutOff) OK=1;
    }
    OK=0;
    while(OK==0)
    {
        ranSurv=normsinv(alea());
        if (abs(ranSurv)<VarianceCutOff) OK=1;
    }

    Individu.ReproductionQuality=sqrt(IndVarRepro)*ranRepro;
    Individu.SurvivalQuality=sqrt(IndVarSurv)*ranSurv;


    while(Alive==1)
        {

            if (LHS==1)// Juvenile
                {
                    repro=0;
                    long double fatum=alea();
                    long double intercept=log(JuvenileSurvival/(1-JuvenileSurvival));
                    long double logitvalue=intercept+Individu.SurvivalQuality;
                    long double Phi=exp(logitvalue)/(1+exp(logitvalue));
                    if(fatum>Phi)
                        {
                            Alive=0;
                        }
                }
            if (LHS==2) // Adult
                {
                    long double logLambda=log(MeanRepro/CorrectionPoisson+Markov*(pastRepro-MeanRepro/CorrectionPoisson))+Individu.ReproductionQuality;
                    long double Lambda=exp(logLambda);
                    double L=exp(-double(Lambda));//generates Poisson distribution
                    int k(0);
                    double p(1.);
                    do{
                    k++;
                    p*=alea();
                    }while(p>L);
                    k--;
                    if(k<0.) {cerr<<endl<<"PATATRAF POISSON, k should not <0, k="<<k<<endl<<flush;k=0.;}
                    repro=k;
                    long double fatum=alea();
                    long double intercept=log(AdultSurvival/(1-AdultSurvival));
                    long double logitvalue=intercept-CorrectionBinomial+(repro-MeanRepro)*CorReproSurv+Individu.SurvivalQuality;
                    long double Phi=exp(logitvalue)/(1+exp(logitvalue));

                    if( ((fatum>Phi) || ((age+Year)>=StudyLength)) || (age>=AgeMax))
                        {
                            Alive=0;
                        }
                    if( ((age+Year)>=StudyLength)) //are we out the study period?
                        {
                            LHS=3;
                        }
                    pastRepro=repro;
                }
            LHind.repro=repro;
            LHind.age=age;
            LHind.survival=Alive;
            LHind.stage=LHS;

            Individu.LifeHistory.push_back(LHind);

            LHS=2;
            age++;
        }
        IndKey++;
    return(Individu);
}//end CIndividu FPoissonSimul(unsigned int& IndKey,long double const& CorrectionPoisson)

/***************************************************************************************************************************/

CIndividu FNeutralSimul(unsigned int& NIndKey, unsigned int const& NYear,double const& RealizedJuvenileSurv,double const& RealizedAdultSurv,vector<long double> const& RealizedARS,unsigned int const& PoissonMaxLRS)
{
    bool Alive(1);
    unsigned int age(1);// means period from age 0 to age 1
    int LHS(1);// 1 means Juvenile, 2 means Adult
    unsigned int repro(0);
    CIndividu Individu;
    CLifeStage LHind;
    Individu.IndividuKey=NIndKey;
    unsigned int LRS(0);
    bool Believe(false);// avoid extrem LRS. Note that we have to cut LRS distribution extreme-tail for KS test

    while(Believe==false)
        {
            LRS=0;
            while(Alive==1)
                {

                    if (LHS==1)// Juvenile
                        {
                            repro=0;
                            long double fatum=alea();
                            if(fatum>RealizedJuvenileSurv)
                                {
                                    Alive=0;
                                }
                        }
                    if (LHS==2) // Adult
                        {
                            long double r=alea();
                            if (r==1.) { r=r-0.00001;}
                            //if (r==0.) { r=r+0.00001;}// very important to avoid cases of k<0, and DO NOT try to solve it by some kind of casting!
                            long double cumul(0.);
                            int k(0);
                            while(cumul<r)
                            {
                                cumul+=RealizedARS[k];
                                k++;
                            }
                            k--;
                            if(k<0) {cerr<<endl<<"ET PATATRAF NEUTRAL, k should not <0, k="<<k<<endl<<flush; k=0.;}
                            repro=k;

                            long double fatum=alea();

                            if( ((fatum>RealizedAdultSurv) || ((age+NYear)>=StudyLength)) || (age>=AgeMax))
                                {
                                    Alive=0;
                                }
                            if (((age+NYear)>=StudyLength))// out study period
                                {
                                    LHS=3;
                                }
                            LRS+=repro;
                        }
                    LHind.repro=repro;
                    LHind.age=age;
                    LHind.survival=Alive;
                    LHind.stage=LHS;
                    Individu.LifeHistory.push_back(LHind);

                    LHS=2;
                    age++;
                }
            if(LRS<=(PoissonMaxLRS*1.1)) {Believe=true; }
            else {Individu.LifeHistory.resize(0); }
        }
    NIndKey++;
    return(Individu);
}//end FNeutralSimul(unsigned int& NIndKey, unsigned int const& NYear)

/***************************************************************************************************************************/

int FDemographer(vector<CCohort> const& Population,long double& H, long double& P)// should estimate entropie and persistence
{
    map<unsigned int, long double> Collector;
    map<unsigned int, map<unsigned int, long double> > Transitions;

    unsigned int MaxARS(1);

    const CCohort *FocalCohort(0);
    const CIndividu *FocalIndividu(0);

    for (unsigned int cohort(0);cohort<Population.size();cohort++)
        {
            FocalCohort=&Population[cohort];
            for (unsigned int ind(0);ind<(*FocalCohort).Individus.size();ind++)
                {
                    FocalIndividu=&(*FocalCohort).Individus[ind];
                    for (unsigned int LS(0);LS<(*FocalIndividu).LifeHistory.size();LS++)//We do not start from 1 because we need anyway a test to exclude juveniles only
                        {
                            if (((*FocalIndividu).LifeHistory[LS].stage==2)&&((*FocalIndividu).LifeHistory.size()>(LS+1)))
                                {
                                    if (MaxARS<(*FocalIndividu).LifeHistory[LS].repro) MaxARS=(*FocalIndividu).LifeHistory[LS].repro;
                                    Collector[(*FocalIndividu).LifeHistory[LS].repro]++;// count the number of transitions starting there
                                    Transitions[(*FocalIndividu).LifeHistory[LS].repro][(*FocalIndividu).LifeHistory[LS+1].repro]++;// count the number of transition going there
                                }
                        }
                }
        }

if(MaxARS>(exp(log(double(MeanRepro))+4*pow(MeanRepro+IndVarRepro,0.5)))) {cerr<<"A problem might have occured in FDemographer, MaxARS="<<MaxARS<<" while maximal prediction is"<<(exp(log(double(MeanRepro))+4*pow(MeanRepro+IndVarRepro,0.5)))<<endl;}

    for (unsigned int i(0);i<=MaxARS;i++)//this creates map cells even when transition does not occur
        {
            Collector[i]++;
            Collector[i]--;
            if(Transitions[i].empty()==false)
                {
                    for (unsigned int j(0);j<=MaxARS;j++)
                        {
                            Transitions[i][j]++;
                            Transitions[i][j]--;
                        }
                }
        }

    long double SumCollector(0.);
    map<unsigned int, long double>::iterator it;

    for (it=Collector.begin();it!=Collector.end();it++)
        {
            map<unsigned int, long double>::iterator it2;
            for (it2=Transitions[(*it).first].begin();it2!=Transitions[(*it).first].end();it2++)
                {
                    if((*it).second!=0)
                        {
                            (*it2).second/=(*it).second;
                        }
                    else // in case a transition smaller than the maximal ARS has 0 occurence (avoids nan)
                        {
                            (*it2).second=0;
                        }
                }
            SumCollector+=(*it).second;
        }

    for (it=Collector.begin();it!=Collector.end();it++)
        {
            (*it).second/=SumCollector;
        }

    H=FEntropy(Transitions,Collector);
    P=FPersistence(Transitions,Collector);
    return 0;
}// end FDemographer(vector<CCohort>& Population,long double& H, long double& P)

/***************************************************************************************************************************/

vector<vector<long double> > FEcologist(double& RealizedAdultSurv,double& RealizedJuvenileSurv, vector<CCohort>& PopulationS)
{
    CCohort *FocalCohort(0);
    CIndividu *FocalIndividu(0);
    double AdNumerator(0.);
    double AdDenumerator(0.);
    double JuvNumerator(0.);
    double JuvDenumerator(0.);

    vector<long double> RealizedARS;
    vector<long double> RealizedLRS;
    vector<vector<long double> > Contain;

    unsigned int maxARS(0);
    unsigned int maxLRS(0);
    unsigned int LRSind(0);
    vector<unsigned int> CountARS;
    vector<unsigned int> CountLRS;
    unsigned int CollectorARS(0);
    unsigned int CollectorLRS(0);

    for (unsigned int cohort(0);cohort<PopulationS.size();cohort++)//frist I need to know the maximal reproductive values
        {
            FocalCohort=&PopulationS[cohort];
            for (unsigned int ind(0);ind<(*FocalCohort).Individus.size();ind++ )
                {
                    FocalIndividu=&(*FocalCohort).Individus[ind];
                    LRSind=0;
                    for (unsigned int LStage(1);LStage<(*FocalIndividu).LifeHistory.size();LStage++ )//we start on the second age (1)
                        {
                            if (maxARS<(*FocalIndividu).LifeHistory[LStage].repro) maxARS=(*FocalIndividu).LifeHistory[LStage].repro;
                            LRSind+=(*FocalIndividu).LifeHistory[LStage].repro;
                        }
                    if(maxLRS<LRSind)
                        {
                            maxLRS=LRSind;
                        }
                }
        }

    RealizedARS.resize(maxARS+1);
    CountARS.resize(maxARS+1);
    RealizedLRS.resize(maxLRS+1);
    CountLRS.resize(maxLRS+1);

    for (unsigned int cohort(0);cohort<PopulationS.size();cohort++)
        {
            FocalCohort=&PopulationS[cohort];
            for (unsigned int ind(0);ind<(*FocalCohort).Individus.size();ind++ )
                {
                    FocalIndividu=&(*FocalCohort).Individus[ind];
                    LRSind=0;
                    for (unsigned int LStage(0);LStage<(*FocalIndividu).LifeHistory.size();LStage++ )
                        {
                            if ( (*FocalIndividu).LifeHistory[LStage].stage==1)//Juvenile
                                {
                                    JuvDenumerator++;
                                    if((*FocalIndividu).LifeHistory[LStage].survival==1)
                                        {
                                            JuvNumerator++;
                                        }
                                }
                            if ( (*FocalIndividu).LifeHistory[LStage].stage==2)//Adult NB: could be 3, if we are out the study period, in which case we don't know survival
                                {
                                    AdDenumerator++;
                                    if((*FocalIndividu).LifeHistory[LStage].survival==1)
                                        {
                                            AdNumerator++;
                                        }
                                    CountARS[(*FocalIndividu).LifeHistory[LStage].repro]++;
                                    CollectorARS++;
                                    LRSind+=(*FocalIndividu).LifeHistory[LStage].repro;
                                }
                            if ( (*FocalIndividu).LifeHistory[LStage].stage==3)// we don't know about survival
                                {
                                    CountARS[(*FocalIndividu).LifeHistory[LStage].repro]++;
                                    CollectorARS++;
                                    LRSind+=(*FocalIndividu).LifeHistory[LStage].repro;
                                }
                        }
                    CountLRS[LRSind]++;
                    CollectorLRS++;
                }
        }

    for (unsigned int i(0);i<=maxARS;i++)
        {
            RealizedARS[i]=double(CountARS[i])/double(CollectorARS);
        }
   for (unsigned int i(0);i<=maxLRS;i++)
        {
            RealizedLRS[i]=double(CountLRS[i])/double(CollectorLRS);
        }
    RealizedAdultSurv=(AdNumerator/AdDenumerator);
    RealizedJuvenileSurv=(JuvNumerator/JuvDenumerator);

    Contain.push_back(RealizedARS);
    Contain.push_back(RealizedLRS);

    return(Contain);
}//end FEcologist(long double& RealizedAdultSurv,long double& RealizedJuvenileSurv,vector<long double>& RealizedARS)

/***************************************************************************************************************************/

vector<long double> FNeutralEcologist(double& NeutralJuvenileSurv,double& NeutralAdultSurv,vector<CCohort>& PopulationN)
{
    vector<long double> LRS;
    CCohort *FocalCohort(0);
    CIndividu *FocalIndividu(0);
    unsigned int maxLRS(0);
    vector<long double> CountLRS;
    long double CollectorLRS(0);
    unsigned int LRSind(0);
    double AdNumerator(0.);
    double AdDenumerator(0.);
    double JuvNumerator(0.);
    double JuvDenumerator(0.);

    for (unsigned int cohort(0);cohort<PopulationN.size();cohort++)//first I need to know the maximal reproductive values
        {
            FocalCohort=&PopulationN[cohort];
            for (unsigned int ind(0);ind<(*FocalCohort).Individus.size();ind++ )
                {
                    FocalIndividu=&(*FocalCohort).Individus[ind];
                    LRSind=0;
                    for (unsigned int LStage(1);LStage<(*FocalIndividu).LifeHistory.size();LStage++ )//we start on the second age (1)
                        {
                            LRSind+=(*FocalIndividu).LifeHistory[LStage].repro;
                        }
                    if(maxLRS<LRSind) {maxLRS=LRSind;}
                }
        }
    LRS.resize(maxLRS+1);
    CountLRS.resize(maxLRS+1);

    for (unsigned int cohort(0);cohort<PopulationN.size();cohort++)//first I need to know the maximal reproductive values
        {
            FocalCohort=&PopulationN[cohort];
            for (unsigned int ind(0);ind<(*FocalCohort).Individus.size();ind++ )
                {
                    FocalIndividu=&(*FocalCohort).Individus[ind];
                    LRSind=0;
                    for (unsigned int LStage(0);LStage<(*FocalIndividu).LifeHistory.size();LStage++ )//we start on the second age (1)
                        {
                             if ( (*FocalIndividu).LifeHistory[LStage].stage==1)//Juvenile
                                {
                                    JuvDenumerator++;
                                    if((*FocalIndividu).LifeHistory[LStage].survival==1)
                                        {
                                            JuvNumerator++;
                                        }
                                }
                            if ( (*FocalIndividu).LifeHistory[LStage].stage==2)//Adult NB: could be 3, if we are out the study period, in which case we don't know survival
                                {
                                    AdDenumerator++;
                                    if((*FocalIndividu).LifeHistory[LStage].survival==1)
                                        {
                                            AdNumerator++;
                                        }
                                    LRSind+=(*FocalIndividu).LifeHistory[LStage].repro;
                                }
                            if ( (*FocalIndividu).LifeHistory[LStage].stage==3)// we don't know about survival
                                {
                                    LRSind+=(*FocalIndividu).LifeHistory[LStage].repro;
                                }
                        }
                    CountLRS[LRSind]++;
                    CollectorLRS++;
                }
        }

   for (unsigned int i(0);i<=maxLRS;i++)
        {
            LRS[i]=CountLRS[i]/CollectorLRS;
        }
    NeutralAdultSurv=(AdNumerator/AdDenumerator);
    NeutralJuvenileSurv=(JuvNumerator/JuvDenumerator);

    return(LRS);
}//end FNeutralEcologist(vector<CCohort> PopulationN)
/***************************************************************************************************************************/

int FAverageLRS(vector<long double>& MeanNeutralLRS,vector<long double>& NeutralLRS, long double& M, long double& V)
{
    M=FMean(NeutralLRS);
    V=FVariance(NeutralLRS,M);
    int diffsize(NeutralLRS.size()-MeanNeutralLRS.size());
    if (diffsize>0)
        {
            for (int i(0);i<(diffsize);i++)
                {
                        MeanNeutralLRS.push_back(0.);
                }
        }
    for (unsigned int i(0);i<NeutralLRS.size();i++)
        {
            MeanNeutralLRS[i]+=(NeutralLRS[i])/NeutralSimulNumber;

        }
    return 0;
}//end FAverageLRS(vector<long double>& MeanNeutralLRS,vector<long double>& NeutralLRS)

/***************************************************************************************************************************/

long double FMean(vector<long double> const& LRSDistri)
{
    long double Mean(0.);
    for (unsigned int i(0);i<LRSDistri.size();i++)
        {
            Mean+=i*LRSDistri[i];
        }
    return(Mean);
}//end FMean(vector<long double> const&LRSDistri)

/***************************************************************************************************************************/
long double FVariance(vector<long double> const& LRSDistri, long double const& Mean)
{
    long double Variance(0.);
    for (unsigned int i(0);i<LRSDistri.size();i++)
        {
            Variance+=pow((Mean-i),2)*LRSDistri[i];
        }
    return(Variance);
}//end FMean(vector<long double> const&LRSDistri)
/***************************************************************************************************************************/

long double FEntropy(map<unsigned int, map<unsigned int, long double> >& Transitions, map<unsigned int, long double>& Collector)//Collector represent the weights and Transition the probabilities
{
    long double H(0.);

    map<unsigned int, long double>::iterator it;
    map<unsigned int, long double>::iterator it2;

    for (it=Collector.begin();it!=Collector.end();it++)
        {

            for (it2=Transitions[(*it).first].begin();it2!=Transitions[(*it).first].end();it2++)
                {
                    if ((*it2).second!=0)
                        {
                            H+=((*it).second)*(-((*it2).second)*(log((*it2).second)))/log(Collector.size());
                        }
                }
        }
    return (H);
}//end FEntropy(map<unsigned int, map<unsigned int, long double> > Transitions)
/***************************************************************************************************************************/
long double FPersistence(map<unsigned int, map<unsigned int, long double> >& Transitions, map<unsigned int, long double>& Collector)
{
    long double P(0.);
    map<unsigned int, long double>::iterator it;
    for (it=Collector.begin();it!=Collector.end();it++)
        {
            P+=(Transitions[(*it).first][(*it).first])*(*it).second;
        }
    return(P);
}//end FPersistence(map<unsigned int, map<unsigned int, long double> >& Transitions, map<unsigned int, long double>& Collector)

/***************************************************************************************************************************/
int FWritingLRSFiles(vector<vector<long double > >& AllMeanNeutralLRS,vector<vector<long double> >& AllRealizedLRS,unsigned int const& AbsoluteMaxLRS)
{
    ofstream LRSFile("LRSdistri.txt");
    LRSFile<<setprecision(5);

    LRSFile<<"#SamplingSeed="<<_ptSamplingSeed<<"\tPoissonSimulNumber="<<PoissonSimulNumber<<"\tNeutralSimulNumber="<<NeutralSimulNumber<<endl;
    LRSFile<<"#AdultSurvival="<<AdultSurvival<<"\tJuvenileSurvival="<<JuvenileSurvival<<"\tIndVarSurv="<<IndVarSurv<<"\tIndVarRepro="<<IndVarRepro<<"\tCorReproSurv="<<CorReproSurv<<endl;
    LRSFile<<"#MeanRepro="<<MeanRepro<<"\tAgeMax="<<AgeMax<<"\tStudyLength="<<StudyLength<<"\tIndInCohort="<<IndInCohort<<endl;

    LRSFile<<"PoissonNB\tNOrS\t";
    for (unsigned int i(0);i<AbsoluteMaxLRS;i++)
        {
                LRSFile<<i<<"\t";
        }
    LRSFile<<endl;

    for (unsigned int Poi(0);Poi<AllRealizedLRS.size();Poi++)
        {
            LRSFile<<Poi<<"\tS\t";
            for (unsigned int SSim(0);SSim<AllRealizedLRS[Poi].size();SSim++)
                {
                    LRSFile<<AllRealizedLRS[Poi][SSim]<<"\t";
                }
            for (unsigned int zeros(AllRealizedLRS[Poi].size());zeros<AbsoluteMaxLRS;zeros++)
                {
                    LRSFile<<0<<"\t";
                }
            LRSFile<<endl;

            LRSFile<<Poi<<"\tN\t";
            for (unsigned int NSim(0);NSim<AllMeanNeutralLRS[Poi].size();NSim++)
                {
                    LRSFile<<AllMeanNeutralLRS[Poi][NSim]<<"\t";
                }
            for (unsigned int zeros(AllMeanNeutralLRS[Poi].size());zeros<AbsoluteMaxLRS;zeros++)
                {
                    LRSFile<<0<<"\t";
                }

            LRSFile<<endl;
        }
    LRSFile.close();
    return 0;
}//end FWrittingLRSFiles(vector<vector<vector<long double> > > AllNeutralLRS,vector<vector<long double> > AllRealizedLRS)

/***************************************************************************************************************************/

int FWritingStatsFiles(CStats& HPMVcontainer)// write file with entropy, persistence, mean LRS and variance in LRS for neutral and Poisson simulations
{
    ofstream StatsFile("Stats.txt");
    StatsFile<<setprecision(5);

    StatsFile<<"#SamplingSeed="<<_ptSamplingSeed<<"\tPoissonSimulNumber="<<PoissonSimulNumber<<"\tNeutralSimulNumber="<<NeutralSimulNumber<<endl;
    StatsFile<<"#AdultSurvival="<<AdultSurvival<<"\tJuvenileSurvival="<<JuvenileSurvival<<"\tIndVarSurv="<<IndVarSurv<<"\tIndVarRepro="<<IndVarRepro<<"\tCorReproSurv="<<CorReproSurv<<endl;
    StatsFile<<"#MeanRepro="<<MeanRepro<<"\tAgeMax="<<AgeMax<<"\tStudyLength="<<StudyLength<<"\tIndInCohort="<<IndInCohort<<endl;

    StatsFile<<"PoissonNB\tSimulNB\tNOrS\tMean\tVariance\tEntropy\tPersistence\tAdultSurvival\tJuvSurvival";
    StatsFile<<endl;
    for (unsigned int Poi(0);Poi<HPMVcontainer.SM.size();Poi++)
        {
            StatsFile<<Poi<<"\t"<<Poi<<"\tS\t";
            StatsFile<<HPMVcontainer.SM[Poi]<<"\t";
            StatsFile<<HPMVcontainer.SV[Poi]<<"\t";
            StatsFile<<HPMVcontainer.SH[Poi]<<"\t";
            StatsFile<<HPMVcontainer.SP[Poi]<<"\t";
            StatsFile<<HPMVcontainer.SAS[Poi]<<"\t";
            StatsFile<<HPMVcontainer.SJS[Poi]<<"\t";

            StatsFile<<endl;
            for (unsigned int Neutral(0);Neutral<HPMVcontainer.NM[Poi].size();Neutral++)
                {
                    StatsFile<<Poi<<"\t"<<Neutral<<"\tN\t";
                    StatsFile<<HPMVcontainer.NM[Poi][Neutral]<<"\t";
                    StatsFile<<HPMVcontainer.NV[Poi][Neutral]<<"\t";
                    StatsFile<<HPMVcontainer.NH[Poi][Neutral]<<"\t";
                    StatsFile<<HPMVcontainer.NP[Poi][Neutral]<<"\t";
                    StatsFile<<HPMVcontainer.NAS[Poi][Neutral]<<"\t";
                    StatsFile<<HPMVcontainer.NJS[Poi][Neutral]<<"\t";
                    StatsFile<<endl;
                }
        }
    StatsFile.close();
    return 0;
}//end WritingStatsFiles(CStats& HPMVcontainer)


/***************************************************************************************************************************/
int FWritingMMFiles(vector<vector<CCohort> >& AllPopulationS)
{
    CCohort *FocalCohort(0);
    CIndividu *FocalIndividu(0);

    ofstream MMFile("MMFile.txt");
    MMFile<<setprecision(5);

    MMFile<<"#SamplingSeed="<<_ptSamplingSeed<<"\tPoissonSimulNumber="<<PoissonSimulNumber<<"\tNeutralSimulNumber="<<NeutralSimulNumber<<endl;
    MMFile<<"#AdultSurvival="<<AdultSurvival<<"\tJuvenileSurvival="<<JuvenileSurvival<<"\tIndVarSurv="<<IndVarSurv<<"\tIndVarRepro="<<IndVarRepro<<"\tCorReproSurv="<<CorReproSurv<<endl;
    MMFile<<"#MeanRepro="<<MeanRepro<<"\tAgeMax="<<AgeMax<<"\tStudyLength="<<StudyLength<<"\tIndInCohort="<<IndInCohort<<endl;


    MMFile<<"Run\tCohort\tIndividual\tAge\tStage\tRepro\tSurvival\tReproQ\tSurvQ"<<endl;

    for (unsigned int run(0);run<AllPopulationS.size();run++)
        {
            for (unsigned int coh(0);coh<AllPopulationS[run].size();coh++)
                {
                    FocalCohort=&AllPopulationS[run][coh];
                    for (unsigned int ind(0);ind<(*FocalCohort).Individus.size();ind++)
                        {
                            FocalIndividu=&(*FocalCohort).Individus[ind];
                            for (unsigned int ls(0);ls<(*FocalIndividu).LifeHistory.size();ls++)
                                {
                                    MMFile<<run<<"\t"<<coh<<"\t"<<(*FocalIndividu).IndividuKey<<"\t"<<ls<<"\t";
                                    MMFile<<(*FocalIndividu).LifeHistory[ls].stage<<"\t";
                                    MMFile<<(*FocalIndividu).LifeHistory[ls].repro<<"\t";
                                    MMFile<<(*FocalIndividu).LifeHistory[ls].survival<<"\t";
                                    MMFile<<(*FocalIndividu).ReproductionQuality<<"\t";
                                    MMFile<<(*FocalIndividu).SurvivalQuality<<"\n";
                                }
                        }
                }
        }
    MMFile.close();
    return 0;
}//end FWritingMMFiles()

/***************************************************************************************************************************/

int FWritingRunInfo()
{
    ofstream RunP ("RunParameters.txt");
    RunP<<"SamplingSeed\tPoissonSimulNumber\tNeutralSimulNumber\tAdultSurvival\tJuvenileSurvival\tIndVarSurv\tIndVarRepro\tCorReproSurv\tMarkov\tMeanRepro\tAgeMax\tStudyLength\tIndInCohort"<<endl;
    RunP<<_ptSamplingSeed<<"\t"<<PoissonSimulNumber<<"\t"<<NeutralSimulNumber<<"\t"<<AdultSurvival<<"\t"<<JuvenileSurvival<<"\t"<<IndVarSurv<<"\t"<<IndVarRepro<<"\t"<<CorReproSurv<<"\t"<<Markov<<"\t"<<MeanRepro<<"\t"<<AgeMax<<"\t"<<StudyLength<<"\t"<<IndInCohort<<endl;
    RunP.close();

    return 0;
}//end int FWritingRunInfo()

/***************************************************************************************************************************/
