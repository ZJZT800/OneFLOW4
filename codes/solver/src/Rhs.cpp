/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Rhs.h"
#include "UNsBcSolver.h"
#include "Zone.h"
#include "DataBase.h"
#include "Iteration.h"
#include "NsCom.h"
#include "UCom.h"
#include "UNsCom.h"
#include "UnsGrid.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "UNsInvFlux.h"
#include "UNsVisFlux.h"
#include "UNsUnsteady.h"
#include "Ctrl.h"

//#include "UINsCorrectPress.h"
//#include "UINsCorrectSpeed.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "UINsInvterm.h"
#include "UINsVisterm.h"
//#include "UINsUnsteady.h"
#include "UINsBcSolver.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Rhs::Rhs()
{
    ;
}

Rhs::~Rhs()
{
    ;
}

void Rhs::UpdateNsResiduals()
{
    NsCmpRHS();
}

void NsCmpBc()
{
    UNsBcSolver * uNsBcSolver = new UNsBcSolver();
    uNsBcSolver->Init();
    uNsBcSolver->CmpBc();
    delete uNsBcSolver;
}

void NSCmpGamaT( int flag )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    ug.Init();
    unsf.Init();
    ug.SetStEd( flag );

    if ( nscom.chemModel == 1 )
    {
    }
    else
    {
        //if ( ZoneState::zid == 0 )
        //{
        //    cout << " ug.ist = " << ug.ist  << " ug.ied = " << ug.ied << "\n";
        //    int kkk = 1;
        //}
        Real oamw = one;
        for ( int cId = ug.ist; cId < ug.ied; ++ cId )
        {
            Real & density  = ( * unsf.q )[ IDX::IR ][ cId ];
            Real & pressure = ( * unsf.q )[ IDX::IP ][ cId ];
            ( * unsf.gama )[ 0 ][ cId ] = nscom.gama_ref;
            ( * unsf.tempr )[ IDX::ITT ][ cId ] = pressure / ( nscom.statecoef * density * oamw );
        }
    }
}

void NsCmpRHS()
{
    NsCmpInvFlux();

    NsCmpVisFlux();

    NsCmpSrcFlux();
}


void NsCmpInvFlux()
{
    UNsInvFlux * uNsInvFlux = new UNsInvFlux();
    uNsInvFlux->CmpFlux();
    delete uNsInvFlux;
}

void NsCmpVisFlux()
{
    UNsVisFlux * uNsVisFlux = new UNsVisFlux();
    uNsVisFlux->CmpFlux();
    delete uNsVisFlux;
}

void NsCmpSrcFlux()
{
    if ( nscom.chemModel == 1 )
    {
        NsCmpChemSrc();
    }

    NsCmpTurbEnergy();

    //dual time step source
    if ( ctrl.idualtime == 1 )
    {
        NsCmpDualTimeStepSrc();
    }
}

void NsCmpChemSrc()
{
    ;
}

void NsCmpTurbEnergy()
{
    ;
}

void NsCmpDualTimeStepSrc()
{
    UNsUnsteady * unsUnsteady = new UNsUnsteady();
    unsUnsteady->CmpDualTimeSrc();
    delete unsUnsteady;
}

void Rhs::UpdateINsResiduals()
{
	//INsCmpRHS();
}

void INsCmpBc()
{
	UINsBcSolver * uINsBcSolver = new UINsBcSolver();
	uINsBcSolver->Init();
	uINsBcSolver->CmpBc();
	delete uINsBcSolver;
}

void INSCmpGamaT(int flag)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	ug.Init();
	uinsf.Init();
	ug.SetStEd(flag);

	if (inscom.chemModel == 1)
	{
	}
	else
	{
		//if ( ZoneState::zid == 0 )
		//{
		//    cout << " ug.ist = " << ug.ist  << " ug.ied = " << ug.ied << "\n";
		//    int kkk = 1;
		//}
		Real oamw = one;
		for (int cId = ug.ist; cId < ug.ied; ++cId)
		{
			Real & density = ( * uinsf.q )[ IIDX::IIR ][ cId ];
			Real & pressure = ( * uinsf.q )[ IIDX::IIP ][ cId ];
		}
	}
}

void Rhs::FieldInit()
{
	INsPreflux();      
}

void Rhs::TrainsAssign()
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.u0[cId] = (*uinsf.q)[IIDX::IIU][cId];
		iinv.v0[cId] = (*uinsf.q)[IIDX::IIV][cId];
		iinv.w0[cId] = (*uinsf.q)[IIDX::IIW][cId];
	}
}

void Rhs::UINsSolver()
{

	    INsCmpVis();

		INsCmpInv(); 

		//INsCmpUnstead(); 

		INsCmpSrc(); 

		DifEqua();    

		Relaxation();

		INsMomPre();

		INsCmpFaceflux();

		INsCorrectPresscoef(); 

		INsCmpPressCorrectEquandUpdatePress();  

		INsCmpSpeedCorrectandUpdateSpeed();  

		INsUpdateFaceflux();   

		INsUpdateRes();

}

void INsCmpTimestep()
{
	
UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpINsTimestep();
	delete uINsInvterm;
}

void INsPreflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpINsPreflux();
	delete uINsInvterm;
}

void INsCmpInv()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpInvcoff();
	delete uINsInvterm;
}

void INsCmpVis()
{
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CmpViscoff();
	delete uINsVisterm;
}

void INsCmpUnstead()
{
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CmpUnsteadcoff();
	delete uINsVisterm;
}

void INsCmpSrc()
{
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CmpINsSrc();
	delete uINsVisterm;
}

void DifEqua()
{
	UINsVisterm* uINsVisterm = new UINsVisterm();
	uINsVisterm->DifEquaMom();
	delete uINsVisterm;
}

void Relaxation()
{
	UINsVisterm* uINsVisterm = new UINsVisterm();
	uINsVisterm->RelaxMom(0.5);
	delete uINsVisterm;
}

void INsMomPre()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->MomPre();
	delete uINsInvterm;
}

void INsCmpFaceflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpFaceflux();
	delete uINsInvterm;
}

void INsCorrectPresscoef()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpCorrectPresscoef();
	delete uINsInvterm;
}

void INsCmpPressCorrectEquandUpdatePress()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpPressCorrectEqu();
	delete uINsInvterm;
}

void INsUpdateFaceflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateFaceflux();
	delete uINsInvterm;
}

void INsCmpSpeedCorrectandUpdateSpeed()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateSpeed();
	delete uINsInvterm;
}

void INsUpdateRes()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateINsRes();
	delete uINsInvterm;
}

//void INsCorrectSpeed()
//{
//	UINsInvterm * uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpCorrectSpeed();
//	delete uINsInvterm;
//}



void INsCmpChemSrc()
{
	;
}

void INsCmpTurbEnergy()
{
	;
}

//void INsCmpDualTimeStepSrc()
//{
//	UINsUnsteady * uinsUnsteady = new UINsUnsteady();
//	uinsUnsteady->CmpDualTimeSrc();
//	delete uinsUnsteady;
//}


EndNameSpace