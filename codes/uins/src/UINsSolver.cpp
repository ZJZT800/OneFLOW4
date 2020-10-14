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

#include "UINsSolver.h"
#include "Mesh.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "SolverInfo.h"
#include "SolverState.h"
#include "FaceTopo.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "DataBase.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "UINsLusgs.h"
#include "INsInvterm.h"
#include "UINsVisterm.h"
#include "UINsInvterm.h"
#include "UCom.h"
#include "Zone.h"
#include "UINsCom.h"
#include <iostream>

using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( UINsSolver )

UINsSolver::UINsSolver()
{
}

UINsSolver::~UINsSolver()
{
}

void UINsSolver::Init()
{
}

void UINsSolver::StaticInit()
{
    INsSolver::StaticInit();
    LusgsState::AddSolver( this->sid, this->gridType, new UINsLusgs() );
}

void UINsSolver::Run()
{
}

//void UINsSolver::uINsSolver()
//{
//	//INsCmpTimestep();
//
//	INsPreflux();
//
//	INsCmpInv(); //���������
//
//	INsCmpVis(); //������ɢ��
//
//	INsCmpUnstead(); //�������̬��
//
//	INsCmpSrc(); //����Դ��Ͷ�������ϵ��
//
//	INsMomPre(); //��⶯������
//
//	INsCmpFaceflux(); //�����������
//
//	INsCorrectPresscoef(); //����ѹ����������ϵ��
//
//	INsCmpPressCorrectEquandUpdatePress();  //��Ҫ��ѹ�����������飬���赥Ԫ����ѹ��δ֪��
//
//	INsCmpSpeedCorrectandUpdateSpeed();  //��Ҫ��������������ٶ�δ֪�����������,���µ�Ԫ�ٶȺ�ѹ��
//
//	INsUpdateFaceflux();   //���½�������
//
//	INsUpdateRes();
//
//}

//void INsPreflux()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpINsPreflux();
//	delete uINsInvterm;
//}
//
//void INsCmpInv()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpInvcoff();
//	delete uINsInvterm;
//}
//
//void INsCmpVis()
//{
//	UINsVisterm* uINsVisterm = new UINsVisterm();
//	uINsVisterm->CmpViscoff();
//	delete uINsVisterm;
//}
//
//void INsCmpUnstead()
//{
//	UINsVisterm* uINsVisterm = new UINsVisterm();
//	uINsVisterm->CmpUnsteadcoff();
//	delete uINsVisterm;
//}
//
//void INsCmpSrc()
//{
//	UINsVisterm* uINsVisterm = new UINsVisterm();
//	uINsVisterm->CmpINsSrc();
//	delete uINsVisterm;
//}
//
//void INsMomPre()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->MomPre();
//	delete uINsInvterm;
//}
//
//void INsCmpFaceflux()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpFaceflux();
//	delete uINsInvterm;
//}
//
//void INsCorrectPresscoef()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpCorrectPresscoef();
//	delete uINsInvterm;
//}
//
//void INsCmpPressCorrectEquandUpdatePress()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpPressCorrectEqu();
//	delete uINsInvterm;
//}
//
//void INsUpdateFaceflux()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->UpdateFaceflux();
//	delete uINsInvterm;
//}
//
//void INsCmpSpeedCorrectandUpdateSpeed()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->UpdateSpeed();
//	delete uINsInvterm;
//}
//
//void INsUpdateRes()
//{
//	UINsInvterm* uINsInvterm = new UINsInvterm();
//	uINsInvterm->UpdateINsRes();
//	delete uINsInvterm;
//}

EndNameSpace