//    DynFixSim simulate animal populations with various levels of fixed and dynamic heterogeneity in life history
//    and extract statistics for further analysis.
//
//	  Copyright (C) 2014  Timoth�e Bonnet - timothee.bonnet@ieu.uzh.ch
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
