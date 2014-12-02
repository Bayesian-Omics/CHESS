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

#include <iostream>

#include "../General_Classes/Kernel_Single_Gamma.h"
#include "../General_Classes/Command_Line.h"
#include "../General_Classes/Model_Generic.h"
#include "../Models/ESS/Model_ESS.h"
#include "../Models/HESS/Model_HESS.h"

//#define GSL_RANGE_CHECK_OFF


using namespace std;


int main(int argc, char *  argv[])
{

    Command_Line My_Command_Line;
    My_Command_Line.Process_Comm_Line(argc,argv);


    //Choice the model and the regression type
    cout << "Model: " << My_Command_Line.Model_Tag << endl;
    cout << "Regression type: " << My_Command_Line.Regression_Model << endl << endl;

    if (My_Command_Line.Model_Tag=="ESS")
    {
            Model_ESS My_Model(My_Command_Line);
            My_Model.Run();
    }
    else if (My_Command_Line.Model_Tag=="HESS")
    {
            Model_HESS My_Model(My_Command_Line);
            My_Model.Run();
    }
    else
    {
        cout << "Unknown Model" << endl;
        exit(1);
    }

    return 0;


}

