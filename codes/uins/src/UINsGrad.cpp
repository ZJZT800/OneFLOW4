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

#include "UINsGrad.h"
#include "UCom.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "INsInvterm.h"
#include "INsVisterm.h"
#include "HXMath.h"
#include "DataBase.h"
#include "FieldImp.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "INsIdx.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"


BeginNameSpace( ONEFLOW )

UINsGrad uins_grad;
UITGrad uit_grad;


UINsGrad::UINsGrad()
{
    ;
}

UINsGrad::~UINsGrad()
{
    ;
}

void UINsGrad::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    //name  = "qf";
    namex = "dqdx";
    namey = "dqdy";
    namez = "dqdz";

	//qf   = GetFieldPointer< MRField > ( grid, name  );
    dqdx = GetFieldPointer< MRField > ( grid, namex );
    dqdy = GetFieldPointer< MRField > ( grid, namey );
    dqdz = GetFieldPointer< MRField > ( grid, namez );

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		(*uinsf.qf)[IIDX::IIU][ug.fId] = iinv.uf[ug.fId];
		(*uinsf.qf)[IIDX::IIV][ug.fId] = iinv.vf[ug.fId];
		(*uinsf.qf)[IIDX::IIW][ug.fId] = iinv.wf[ug.fId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		(*uinsf.qf)[IIDX::IIU][ug.fId] = cl * (*uinsf.q)[IIDX::IIU][ug.lc] + cr * (*uinsf.q)[IIDX::IIU][ug.rc];;
		(*uinsf.qf)[IIDX::IIV][ug.fId] = cl * (*uinsf.q)[IIDX::IIV][ug.lc] + cr * (*uinsf.q)[IIDX::IIV][ug.rc];;
		(*uinsf.qf)[IIDX::IIW][ug.fId] = cl * (*uinsf.q)[IIDX::IIW][ug.lc] + cr * (*uinsf.q)[IIDX::IIW][ug.rc];;
	}
    //bdqdx = GetFieldPointer< MRField > ( grid, "bcdqdx" );
    //bdqdy = GetFieldPointer< MRField > ( grid, "bcdqdy" );
    //bdqdz = GetFieldPointer< MRField > ( grid, "bcdqdz" );

    //this->nEqu = inscom.nTEqu;

    //this->istore = 1;
}

UITGrad::UITGrad()
{
	;
}

UITGrad::~UITGrad()
{
	;
}

void UITGrad::Init()
{
	UnsGrid * grid = Zone::GetUnsGrid();

	name = "tempr";
	namex = "dtdx";
	namey = "dtdy";
	namez = "dtdz";

	q = GetFieldPointer< MRField >(grid, name);
	dqdx = GetFieldPointer< MRField >(grid, namex);
	dqdy = GetFieldPointer< MRField >(grid, namey);
	dqdz = GetFieldPointer< MRField >(grid, namez);

	this->nEqu = inscom.nTModel;

	this->istore = 0;
}





EndNameSpace