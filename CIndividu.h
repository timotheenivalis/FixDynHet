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
