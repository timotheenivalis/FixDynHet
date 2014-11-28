#ifndef CINDIVIDU_H_INCLUDED
#define CINDIVIDU_H_INCLUDED

#include <vector>
#include "CLifeStage.h"

class CIndividu
{
    public:
    CIndividu();
    unsigned int IndividuKey;
    long double ReproductionQuality;
    long double SurvivalQuality;
    std::vector<CLifeStage> LifeHistory;
};
#endif // CINDIVIDU_H_INCLUDED
