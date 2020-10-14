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

#include "UINsVisterm.h"
#include "INsInvterm.h"
#include "UINsInvterm.h"
#include "INsVisterm.h"
#include "Iteration.h"
#include "HeatFlux.h"
#include "UGrad.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "INsCtrl.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "VisGrad.h"
#include "UINsGrad.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "BcData.h"
#include "ULimiter.h"
#include "UINsLimiter.h"
#include "FieldImp.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include <iostream>
using namespace std;

BeginNameSpace(ONEFLOW)


UINsVisterm::UINsVisterm()
{
	;
}

UINsVisterm::~UINsVisterm()
{
	;
}


void UINsVisterm::CmpViscoff()
{
	if (vis_model.vismodel == 0) return;

	this->CmpVisterm();

}

void UINsVisterm::Alloc()
{
	uinsf.qf = new MRField(inscom.nEqu, ug.nFace);
}

void UINsVisterm::DeAlloc()
{
	delete uinsf.qf;
}

void UINsVisterm::PrepareField()
{

	ONEFLOW::CmpINsGrad((*uinsf.q)[IIDX::IIU], (*uinsf.dqdx)[IIDX::IIU], (*uinsf.dqdy)[IIDX::IIU], (*uinsf.dqdz)[IIDX::IIU]);
	ONEFLOW::CmpINsGrad((*uinsf.q)[IIDX::IIV], (*uinsf.dqdx)[IIDX::IIV], (*uinsf.dqdy)[IIDX::IIV], (*uinsf.dqdz)[IIDX::IIV]);
	ONEFLOW::CmpINsGrad((*uinsf.q)[IIDX::IIW], (*uinsf.dqdx)[IIDX::IIW], (*uinsf.dqdy)[IIDX::IIW], (*uinsf.dqdz)[IIDX::IIW]);

}

void UINsVisterm::CmpPreandVisGrad()
{
	(*uinsf.dqdx)[IIDX::IIR] = 0;
	(*uinsf.dqdy)[IIDX::IIR] = 0;
	(*uinsf.dqdz)[IIDX::IIR] = 0;

	(*uinsf.dqdx)[IIDX::IIU] = 0;
	(*uinsf.dqdy)[IIDX::IIU] = 0;
	(*uinsf.dqdz)[IIDX::IIU] = 0;

	(*uinsf.dqdx)[IIDX::IIV] = 0;
	(*uinsf.dqdy)[IIDX::IIV] = 0;
	(*uinsf.dqdz)[IIDX::IIV] = 0;

	(*uinsf.dqdx)[IIDX::IIW] = 0;
	(*uinsf.dqdy)[IIDX::IIW] = 0;
	(*uinsf.dqdz)[IIDX::IIW] = 0;

	(*uinsf.dqdx)[IIDX::IIP] = 0;
	(*uinsf.dqdy)[IIDX::IIP] = 0;
	(*uinsf.dqdz)[IIDX::IIP] = 0;


	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		if (fId == 432)
		{
			int kkk = 1;
		}
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

		Real value1 = cl * (*uinsf.q)[IIDX::IIU][ug.lc] + cr * (*uinsf.q)[IIDX::IIU][ug.rc];
		Real value2 = cl * (*uinsf.q)[IIDX::IIV][ug.lc] + cr * (*uinsf.q)[IIDX::IIV][ug.rc];
		Real value3 = cl * (*uinsf.q)[IIDX::IIW][ug.lc] + cr * (*uinsf.q)[IIDX::IIW][ug.rc];
		Real value4 = cl * (*uinsf.q)[IIDX::IIP][ug.lc] + cr * (*uinsf.q)[IIDX::IIP][ug.rc];

		/*Real value1 = iinv.uf[ug.fId];
		Real value2 = iinv.vf[ug.fId];
		Real value3 = iinv.wf[ug.fId];
		Real value4 = iinv.pf[ug.fId];*/

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		(*uinsf.dqdx)[IIDX::IIU][ug.lc] += fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.lc] += fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.lc] += fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.lc] += fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.lc] += fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.lc] += fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.lc] += fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.lc] += fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.lc] += fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.lc] += fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.lc] += fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.lc] += fnza * value4;

		//if (ug.fId < ug.nBFace) continue;

		(*uinsf.dqdx)[IIDX::IIU][ug.rc] += -fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.rc] += -fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.rc] += -fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.rc] += -fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.rc] += -fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.rc] += -fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.rc] += -fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.rc] += -fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.rc] += -fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] += -fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.rc] += -fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.rc] += -fnza * value4;
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		if (fId == 432)
		{
			int kkk = 1;
		}
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real value1 = iinv.uf[ug.fId];
		Real value2 = iinv.vf[ug.fId];
		Real value3 = iinv.wf[ug.fId];
		Real value4 = iinv.pf[ug.fId];

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		(*uinsf.dqdx)[IIDX::IIU][ug.lc] += fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.lc] += fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.lc] += fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.lc] += fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.lc] += fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.lc] += fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.lc] += fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.lc] += fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.lc] += fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.lc] += fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.lc] += fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.lc] += fnza * value4;
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		Real ovol = one / (*ug.cvol)[ug.cId];
		(*uinsf.dqdx)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIP][ug.cId] *= ovol;
	}

}


void UINsVisterm::CmpVisterm()
{
	RealField dudx, dudy, dudz;
	RealField dvdx, dvdy, dvdz;
	RealField dwdx, dwdy, dwdz;

	dudx.resize(ug.nCell);
	dudy.resize(ug.nCell);
	dudz.resize(ug.nCell);
	dvdx.resize(ug.nCell);
	dvdy.resize(ug.nCell);
	dvdz.resize(ug.nCell);
	dwdx.resize(ug.nCell);
	dwdy.resize(ug.nCell);
	dwdz.resize(ug.nCell);

	ONEFLOW::CmpINsGrad(iinv.uf, dudx, dudy, dudz);
	ONEFLOW::CmpINsGrad(iinv.vf, dvdx, dvdy, dvdz);
	ONEFLOW::CmpINsGrad(iinv.wf, dwdx, dwdy, dwdz);

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[fId];
		ug.rc = (*ug.rcf)[fId];
        
		this->CmpFaceVisterm(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);

	}

	//Direchlet Boundary Condition
	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			ug.bcfId = ibc;

			BcInfo* bcInfo = ug.bcRecord->bcInfo;

			ug.fId = bcInfo->bcFace[ug.ir][ibc];
			ug.bcr = bcInfo->bcRegion[ug.ir][ibc];
			ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

			if (ug.bcr == -1) return; //interface
			int dd = bcdata.r2d[ug.bcr];
			if (dd != -1)
			{
				ug.bcdtkey = 1;
				//inscom.bcflow = &bcdata.dataList[dd];
			}
			//if (ug.bcdtkey == 1)
			//{
			ug.lc = (*ug.lcf)[ug.fId];
			this->CmpBcFaceVisterm(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
			//}
		}
	}
	

}

void UINsVisterm::CmpVisCoef(Real& vis)
{
	if (ug.fId < ug.nBFace)
	{
		iinv.rf = (*uinsf.q)[IIDX::IIR][ug.lc] * (*ug.fl)[ug.fId] + (*uinsf.q)[IIDX::IIR][ug.rc] * (1 - (*ug.fl)[ug.fId]);
		Real v = sqrt(iinv.uf[ug.fId] * iinv.uf[ug.fId] + iinv.vf[ug.fId] * iinv.vf[ug.fId] + iinv.wf[ug.fId] * iinv.wf[ug.fId]);
		vis = iinv.rf * v * 1 / inscom.reynolds;
	}
	else if (ug.fId > ug.nBFace - 1)
	{
		iinv.rf = (*uinsf.q)[IIDX::IIR][ug.lc];
		Real v = sqrt(iinv.uf[ug.fId] * iinv.uf[ug.fId] + iinv.vf[ug.fId] * iinv.vf[ug.fId] + iinv.wf[ug.fId] * iinv.wf[ug.fId]);
		vis = iinv.rf * v * 1 / inscom.reynolds;
	}
}

void UINsVisterm::CmpFaceVisterm(RealField & dudx, RealField & dudy, RealField & dudz, RealField & dvdx, RealField & dvdy, RealField & dvdz, RealField & dwdx, RealField& dwdy, RealField& dwdz)
{

	Real l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];  
	Real l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	Real l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

	Real vis = 1.0 / inscom.reynolds;
	//CmpVisCoef(vis);

	Real dist = (*ug.a1)[ug.fId] * l2rdx + (*ug.a2)[ug.fId] * l2rdy + (*ug.a3)[ug.fId] * l2rdz;

	Real Fn = (*ug.a1)[ug.fId] * (*ug.a1)[ug.fId] + (*ug.a2)[ug.fId] * (*ug.a2)[ug.fId] + (*ug.a3)[ug.fId] * (*ug.a3)[ug.fId];

	Fn = Fn / dist;

	Real T1 = (*ug.a1)[ug.fId] - l2rdx * Fn;
	Real T2 = (*ug.a2)[ug.fId] - l2rdy * Fn;
	Real T3 = (*ug.a3)[ug.fId] - l2rdz * Fn;

	Real fdudx = (*ug.fl)[ug.fId] * dudx[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dudx[ug.rc];
	Real fdudy = (*ug.fl)[ug.fId] * dudy[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dudy[ug.rc];
	Real fdudz = (*ug.fl)[ug.fId] * dudz[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dudz[ug.rc];
	Real fdvdx = (*ug.fl)[ug.fId] * dvdx[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dvdx[ug.rc];
	Real fdvdy = (*ug.fl)[ug.fId] * dvdy[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dvdy[ug.rc];
	Real fdvdz = (*ug.fl)[ug.fId] * dvdz[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dvdz[ug.rc];
	Real fdwdx = (*ug.fl)[ug.fId] * dwdx[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dwdx[ug.rc];
	Real fdwdy = (*ug.fl)[ug.fId] * dwdy[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dwdy[ug.rc];
	Real fdwdz = (*ug.fl)[ug.fId] * dwdz[ug.lc] + (1 - (*ug.fl)[ug.fId]) * dwdz[ug.rc];

	iinv.ai[ug.fId][0] += vis * Fn;
	iinv.ai[ug.fId][1] += vis * Fn;

	iinv.buc[ug.lc] += vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);
	iinv.buc[ug.rc] -= vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);

	iinv.bvc[ug.lc] += vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
	iinv.bvc[ug.rc] -= vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);

	iinv.bwc[ug.lc] += vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
	iinv.bwc[ug.rc] -= vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
}

void UINsVisterm::CmpBcFaceVisterm(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz)
{
	Real l2rdx = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];  
	Real l2rdy = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real l2rdz = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dist = (*ug.a1)[ug.fId] * l2rdx + (*ug.a2)[ug.fId] * l2rdy + (*ug.a3)[ug.fId] * l2rdz;

	Real Fn = (*ug.a1)[ug.fId] * (*ug.a1)[ug.fId] + (*ug.a2)[ug.fId] * (*ug.a2)[ug.fId] + (*ug.a3)[ug.fId] * (*ug.a3)[ug.fId];

	Fn = Fn / dist;

	Real vis = 1.0 / inscom.reynolds;
	//CmpVisCoef(vis);

	Real T1 = (*ug.a1)[ug.fId] - l2rdx * Fn;
	Real T2 = (*ug.a2)[ug.fId] - l2rdy * Fn;
	Real T3 = (*ug.a3)[ug.fId] - l2rdz * Fn;

	Real fdudx = dudx[ug.lc];
	Real fdudy = dudy[ug.lc];
	Real fdudz = dudz[ug.lc];
	Real fdvdx = dvdx[ug.lc];
	Real fdvdy = dvdy[ug.lc];
	Real fdvdz = dvdz[ug.lc];
	Real fdwdx = dwdx[ug.lc];
	Real fdwdy = dwdy[ug.lc];
	Real fdwdz = dwdz[ug.lc];

	iinv.rf = iinv.rl;

	iinv.spc[ug.lc] += vis * Fn;

	iinv.buc[ug.lc] += vis * Fn * iinv.uf[ug.fId] + vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);
	iinv.bvc[ug.lc] += vis * Fn * iinv.vf[ug.fId] + vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
	iinv.bwc[ug.lc] += vis * Fn * iinv.wf[ug.fId] + vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);

}


void UINsVisterm::CmpUnsteadcoff()
{
	iinv.timestep = GetDataValue< Real >("global_dt");

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		iinv.spc[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][cId] / iinv.timestep;  

		iinv.buc[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][cId] * iinv.u0[cId] / iinv.timestep;
		iinv.bvc[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][cId] * iinv.v0[cId] / iinv.timestep;
		iinv.bwc[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][cId] * iinv.w0[cId] / iinv.timestep;
	}

}


void UINsVisterm::CmpINsSrc()
{

	RealField dpdx, dpdy, dpdz;
	dpdx.resize(ug.nCell);
	dpdy.resize(ug.nCell);
	dpdz.resize(ug.nCell);
	ONEFLOW::CmpINsGrad(iinv.pf, dpdx, dpdy, dpdz);

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Real vol = (*ug.cvol)[cId];
		iinv.buc[cId] = iinv.buc[cId] - vol * dpdx[cId];
		iinv.bvc[cId] = iinv.bvc[cId] - vol * dpdy[cId];
		iinv.bwc[cId] = iinv.bwc[cId] - vol * dpdz[cId];
		//iinv.spc[cId] += iinv.spt[cId];
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.spc[lc] += iinv.ai[fId][0];
		iinv.spc[rc] += iinv.ai[fId][1];
	}
}

void UINsVisterm::DifEquaMom()
{
	iinv.remax_up = 0;
	iinv.remax_vp = 0;
	iinv.remax_wp = 0;
	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.buc[lc] += iinv.ai[fId][0] * (*uinsf.q)[IIDX::IIU][rc];
		iinv.bvc[lc] += iinv.ai[fId][0] * (*uinsf.q)[IIDX::IIV][rc];
		iinv.bwc[lc] += iinv.ai[fId][0] * (*uinsf.q)[IIDX::IIW][rc];

		iinv.buc[rc] += iinv.ai[fId][1] * (*uinsf.q)[IIDX::IIU][lc];
		iinv.bvc[rc] += iinv.ai[fId][1] * (*uinsf.q)[IIDX::IIV][lc];
		iinv.bwc[rc] += iinv.ai[fId][1] * (*uinsf.q)[IIDX::IIW][lc];
			
	}

	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.buc[cId] -= iinv.spc[cId] * (*uinsf.q)[IIDX::IIU][cId];
		iinv.bvc[cId] -= iinv.spc[cId] * (*uinsf.q)[IIDX::IIV][cId];
		iinv.bwc[cId] -= iinv.spc[cId] * (*uinsf.q)[IIDX::IIW][cId];
		iinv.remax_up = MAX(abs(iinv.remax_up), abs(iinv.buc[cId]));
		iinv.remax_vp = MAX(abs(iinv.remax_vp), abs(iinv.bvc[cId]));
		iinv.remax_wp = MAX(abs(iinv.remax_wp), abs(iinv.bwc[cId]));
	}
}

void UINsVisterm::RelaxMom(Real a)
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.spc[cId] = iinv.spc[cId] * (1 + a);
	}
}

EndNameSpace

