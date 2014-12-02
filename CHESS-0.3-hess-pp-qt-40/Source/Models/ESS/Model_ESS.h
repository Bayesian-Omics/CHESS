/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * CHESS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CHESS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CHESS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MODEL_ESS_H
#define MODEL_ESS_H

#include "../../General_Classes/Model_Generic.h"
#include "Start_Up_ESS.h"

#include <algorithm>
#include <string>
#include <vector>

#include "../../Kernel/Routines/struc.h"
#include "../../Kernel/Routines/dyn_name.h"
#include "../../Kernel/Routines/matrix_handling.h"
#include "../../Kernel/Routines/rand.h"
#include "../../Kernel/Routines/moves.h"
#include "../../Kernel/Routines/regression.h"
#include "../../Kernel/Routines/cond_post.h"
#include "../../Kernel/Routines/post_processing.h"
#include "../../Kernel/Routines/xml_file_read.h"
#include "../../Kernel/Classes/String_Matrices.h"
#include "../../Kernel/Classes/Int_Matrices.h"
#include "../../Kernel/Classes/Int_Matrices_var_dim.h"
#include "../../Kernel/Classes/Double_Matrices.h"
#include "../../Kernel/Classes/Double_Matrices_cont.h"
#include "../../Kernel/Classes/Prior_param.h"
#include "../../Kernel/Classes/Temperatures.h"
#include "../../Kernel/Classes/Move_monitor.h"
#include "../../Kernel/Classes/AdMH.h"
#include "../../Kernel/Classes/DR.h"
#include "../../Kernel/Classes/CM.h"


using namespace std;

class Model_ESS : public Model_Generic
{
public:
    Model_ESS();
    //Inherited constructor:
    Model_ESS(Command_Line My_Command_Line) : Model_Generic(My_Command_Line){}


    /*_______________
        Overloaded:
      ---------------*/
    void Run();


};


#endif // MODEL_ESS_H
