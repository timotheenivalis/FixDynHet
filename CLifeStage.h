#ifndef CLIFESTAGE_H_INCLUDED
#define CLIFESTAGE_H_INCLUDED

#include <cstring>
#include <string>

class CLifeStage
{
    public:
    CLifeStage();
    unsigned int age;
    bool survival;
    unsigned int repro;
    int stage;// 1 means Juvenile and 2 Adult. We do not code it by a boolean because it is succeptible to be change to more than 2 stages
};

#endif // CLIFESTAGE_H_INCLUDED
