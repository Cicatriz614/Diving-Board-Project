
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>

//#include <Engine.h>

//#pragma comment (lib, "libmat.lib")
//#pragma comment (lib, "libmx.lib")
//#pragma comment (lib, "libmex.lib")
//#pragma comment (lib, "libeng.lib")


//#pragma once

class CirMain
{
	// operations
	public:
		CirMain(void);
		CirMain(int nodes, double E);
		~CirMain(void);

		void Stiff_Multip204	();
		void Stiff_Single204	(int i_node, double d_Length, double d_CrArea, double d_MomInt, double d_GloAgl);

		void Stiff_MultipHut	();
		void Stiff_SingleHut	(int i_node, double d_Length, double d_CrArea, double d_MomInt, double d_GloAgl);

		void Trans_PosIF ();
		void Trans_PosAg ();
		void Trans_PosLg ();

		void Trans_Multip	();
		void Trans_Single	(int i_node, double d_GloAgl);

		double Trans_StsTMtx	();

		void Cin_PosIMtx(double* inipos);
		void Cin_DisMtx(double* inidis);

        void Cin_AreMtx(double* area);
        void Cout_AreMtx();

        void Cin_AMIMtx(double* areamoment);
        void Cout_AMIMtx();

		void Cout_PosIMtx();
		void Cout_PosFMtx();

		void Cout_PosLMtx();
		void Cout_PosAMtx();
		void Cout_DisMtx();
		void Cout_ForMtx();

		void Cout_DisLMtx();
		void Cout_ForLMtx();

		void Cout_StsAMtx();
		void Cout_StsBMtx();
		void Cout_StsTMtx();

	// variables
	public:
		int		i_nodenm;
		double	d_ElaMod;

		double* d_PosIMtx;
		double* d_PosFMtx;

		double* d_PosLIMtx;
		double* d_PosLFMtx;
		double* d_PosAMtx;

        double* d_AreMtx;
        double* d_AMIMtx;

		double*	d_DisMtx;
		double*	d_ForMtx;

		double*	d_DisLMtx;
		double*	d_ForLMtx;

		double*	d_StsAMtx;
		double*	d_StsBMtx;
		double*	d_StsTMtx;

	// accessors
	public:
		int		a_nodenm(void)	{return i_nodenm;}
		bool	a_ElaMod(void) 	{return d_ElaMod;}
};

///*/
