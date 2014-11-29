#ifndef CSTATS_H_INCLUDED
#define CSTATS_H_INCLUDED

#include <vector>

class CStats
{
    public:
    CStats();
    std::vector<std::vector<long double> > NH;
    std::vector<std::vector<long double> > NP;
    std::vector<std::vector<long double> > NM;
    std::vector<std::vector<long double> > NV;
    std::vector<std::vector<long double> > NAS;// adult survival
    std::vector<std::vector<long double> > NJS;// juvenile survival
    std::vector<long double> SH;
    std::vector<long double> SP;
    std::vector<long double> SM;
    std::vector<long double> SV;
    std::vector<long double> SAS;
    std::vector<long double> SJS;
};


#endif // CSTATS_H_INCLUDED
