#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "matrixops.h"
#include "node.h"
#include "CirMain.h"
/*
#define DTau (2.5e-5)
#define DTSQUARED (6.25e-10)
#define DAMPING 1
#define STIFFNESS 100
*/

using namespace std;

//Global Variables
double DTau;
double TFinal;
double DTauSqd;
double TwiceDTau;
double *CurDisplacement;
double *LastDisplacement;
double *NextDisplacement;
double ModulusofElasticity = 7000000000; //About right for 7000 series alloy
const double Bull = 5;
const double dampfac = 1.26044887;
const double density = 2800; //kg/m3

//select chooses 0 (stiffness) 1 (dampness) or 2 (mass/moment)
void localwbending(double**k, Node*list, int first, int sec)
{
    //int number = i_node;
	double l = list[first].get_init_length(list[sec]);
	double a = list[first].GetConArea(0);
	double i = list[first].GetAreaMoment();
	double c = list[first].getFrac(list[sec],0,1);
	double s = list[first].getFrac(list[sec],1,0);

	double al2 = a*l*l;
	double i12 = 12*i;
	double e13;

    e13 = ModulusofElasticity/l/l/l;
	//double el3 = d_ElaMod/l/l/l;

	double c2 = c*c;
	double s2 = s*s;
	double cs = c*s;

//	printf("%.3f, %.3f", c,s);
/*
	k =
	{
//		{			0,				1,			2,				3,					4,		5,	}
		{  i12*s2 +al2*c2,	 (al2-i12)*cs, -6*i*l*s, -i12*s2 -al2*c2,	  (i12-al2)*cs, -6*i*l*s},	// row 0
		{	 (al2-i12)*cs, i12*c2 +al2*s2,  6*i*l*c,	(i12-al2)*cs,  -i12*c2 -al2*s2,  6*i*l*c},	// row 1
		{		 -6*i*l*s,		  6*i*l*c,  4*i*l*l,		 6*i*l*s,		  -6*i*l*c,  2*i*l*l},	// row 2
		{ -i12*s2 -al2*c2,	 (i12-al2)*cs,	6*i*l*s,  i12*s2 +al2*c2,	  (al2-i12)*cs,  6*i*l*s},	// row 3
		{	 (i12-al2)*cs,-i12*c2 -al2*s2, -6*i*l*c,	(al2-i12)*cs,	i12*c2 +al2*s2, -6*i*l*c},	// row 4
		{		 -6*i*l*s,		  6*i*l*c,	2*i*l*l,		 6*i*l*s,		  -6*i*l*c,  4*i*l*l},	// row 5
	};
*/
    k[0][0] = i12*s2 +al2*c2;
    k[0][1] = (al2-i12)*cs;
    k[0][2] = -6*i*l*s;
    k[0][3] = -i12*s2 -al2*c2;
    k[0][4] = (i12-al2)*cs;
    k[0][5] = -6*i*l*s;
    k[1][0] = (al2-i12)*cs;
    k[1][1] = i12*c2 +al2*s2;
    k[1][2] = 6*i*l*c;
    k[1][3] = (i12-al2)*cs;
    k[1][4] = -i12*c2 -al2*s2;
    k[1][5] = 6*i*l*c;
    k[2][0] =  -6*i*l*s;
    k[2][1] = 6*i*l*c;
    k[2][2] = 4*i*l*l;
    k[2][3] = 6*i*l*s;
    k[2][4] = -6*i*l*c;
    k[2][5] = 2*i*l*l;
    k[3][0] = -i12*s2 -al2*c2;
    k[3][1] = (i12-al2)*cs;
    k[3][2] = 6*i*l*s;
    k[3][3] = i12*s2 +al2*c2;
    k[3][4] = (al2-i12)*cs;
    k[3][5] = 6*i*l*s;
    k[4][0] = (i12-al2)*cs;
    k[4][1] = -i12*c2 -al2*s2;
    k[4][2] = -6*i*l*c;
    k[4][3] = (al2-i12)*cs;
    k[4][4] = i12*c2 +al2*s2;
    k[4][5] = -6*i*l*c;
    k[5][0] = -6*i*l*s;
    k[5][1] = 6*i*l*c;
    k[5][2] = 2*i*l*l;
    k[5][3] = 6*i*l*s;
    k[5][4] = -6*i*l*c;
    k[5][5] = 4*i*l*l;

	for (int f = 0; f < 6; f++)
	{
		for (int g = 0; g < 6; g++)
		{
			k[f][g] *= e13;
		}

	}
    //MPrint(k,size,size);
}

void KAssemble(double**k,Node*list,int NCNT)
{
    double**local = CreateMatrix(6);
    for(int i=0;i<NCNT;i++)
    {
        for(int j=0;j<list[i].getconn();j++)
        {
            localwbending(local,list,i,list[i].connto(0));
            addLocToGlo(k,local,i,list[i].connto(0),3);
        }
    }
    DeleteMatrix(6,local);
}

void localc(double**c,Node*list,int first,int sec)
{
    double a = list[first].getFrac(list[sec],0,1);
	double b = list[first].getFrac(list[sec],1,0);
	double D = dampfac*list[first].getconnlength(0);

/*
	c =
	{
//		{			0,				1,			2,				3,					4,		5,	}
		{  Da2,	 Dab, 0, -Da2,	  -Dab, 0},	// row 0
		{  Dab,	 Db2, 0, -Dab,	  -Db2, 0},	// row 1
		{  0  ,    0, B,    0,       0,-B},	// row 2
		{ -Da2,	-Dab, 0,  Da2,	   Dab, 0},	// row 3
		{ -Dab,	-Db2, 0,  Dab,	   Db2, 0},	// row 4
		{  0  ,    0,-B,    0,       0, B},	// row 5
	};
*/

    c[0][0] = D*a*a;
    c[0][1] = D*a*b;
    c[0][2] = 0;
    c[0][3] = -D*a*a;
    c[0][4] = -D*a*b;
    c[0][5] = 0;
    c[1][0] = D*a*b;
    c[1][1] = D*b*b;
    c[1][2] = 0;
    c[1][3] = -D*a*b;
    c[1][4] = -D*b*b;
    c[1][5] = 0;
    c[2][0] = 0;
    c[2][1] = 0;
    c[2][2] = Bull;
    c[2][3] = 0;
    c[2][4] = 0;
    c[2][5] = -Bull;
    c[3][0] = -D*a*a;
    c[3][1] = -D*a*b;
    c[3][2] = 0;
    c[3][3] = D*a*a;
    c[3][4] = D*a*b;
    c[3][5] = 0;
    c[4][0] = -D*a*b;
    c[4][1] = -D*b*b;
    c[4][2] = 0;
    c[4][3] = D*a*b;
    c[4][4] = D*b*b;
    c[4][5] = 0;
    c[5][0] = 0;
    c[5][1] = 0;
    c[5][2] = -Bull;
    c[5][3] = 0;
    c[5][4] = 0;
    c[5][5] = Bull;
}

void CAssemble(double**c,Node*list,int NCNT)
{
    double**local = CreateMatrix(6);
    for(int i=0;i<NCNT;i++)
    {
        for(int j=0;j<list[i].getconn();j++)
        {
            localc(local,list,i,list[i].connto(0));
            addLocToGlo(c,local,i,list[i].connto(0),3);
        }
    }
    DeleteMatrix(6,local);
}

void localm(double**k,Node*list,int first,int sec)
{
    double factor = (density*(list[first].GetConArea(0))*list[first].getconnlength(0))/420;
    double c = list[first].getFrac(list[sec],0,1);
    double s = list[first].getFrac(list[sec],1,0);

    k[0][0] = factor*(140*c*c+15*s*s);
    k[0][1] = -16*c*s*factor;
    k[0][2] = -22*s*list[first].getconnlength(0)*factor;
    k[0][3] = factor*(70*c*c+54*s*s);
    k[0][4] = 16*c*s*factor;
    k[0][5] = 13*s*list[first].getconnlength(0)*factor;
    k[1][0] = -16*c*s*factor;
    k[1][1] = (140*s*s+156*c*c)*factor;
    k[1][2] = factor*(22*c*list[first].getconnlength(0));
    k[1][3] = factor*16*c*s;
    k[1][4] = factor*(70*s*s*54*c*c);
    k[1][5] = -13*(c)*factor*list[first].getconnlength(0);
    k[2][0] = -22*s*list[first].getconnlength(0)*factor;
    k[2][1] = factor*(22*c*list[first].getconnlength(0));
    k[2][2] = 4*(list[first].getconnlength(0))*(list[first].getconnlength(0))*factor;
    k[2][3] = -13*s*list[first].getconnlength(0)*factor;
    k[2][4] = 13*c*(list[first].getconnlength(0))*factor;
    k[2][5] = -3*(list[first].getconnlength(0))*(list[first].getconnlength(0))*factor;
    k[3][0] = factor*(70*c*c+54*s*s);
    k[3][1] = factor*16*c*s;
    k[3][2] = -13*s*list[first].getconnlength(0)*factor;
    k[3][3] = (140*c*c+156*s*s)*factor;
    k[3][4] = -16*c*s*factor;
    k[3][5] = 22*s*list[first].getconnlength(0)*factor;
    k[4][0] = 16*c*s*factor;
    k[4][1] = factor*(70*s*s*54*c*c);
    k[4][2] = 13*c*(list[first].getconnlength(0))*factor;
    k[4][3] = -16*c*s*factor;
    k[4][4] = (140*s*s+156*c*c)*factor;
    k[4][5] = -22*c*(list[first].getconnlength(0))*factor;
    k[5][0] = 13*s*list[first].getconnlength(0)*factor;
    k[5][1] = -13*(c)*factor*list[first].getconnlength(0);
    k[5][2] = -3*(list[first].getconnlength(0))*(list[first].getconnlength(0))*factor;
    k[5][3] = 22*s*list[first].getconnlength(0)*factor;
    k[5][4] = -22*c*(list[first].getconnlength(0))*factor;
    k[5][5] = 4*(list[first].getconnlength(0))*(list[first].getconnlength(0))*factor;
}

void MAssemble(double**m,Node*list,int NCNT)
{
    double**local = CreateMatrix(6);
    for(int i=0;i<NCNT;i++)
    {
        for(int j=0;j<list[i].getconn();j++)
        {
            localm(local,list,i,list[i].connto(0));
            addLocToGlo(m,local,i,list[i].connto(0),3);
        }
    }
    DeleteMatrix(6,local);
}
//Select is used to choose double derivative, derivative, or proportional (2,1, or 0)
//It is assumed that the dependant variable (ex. force) is only transferred axially
//transformation matrix is passed as a parameter so that it can be used multiple
//times without recalculating


// matrices named for M, C, and K for visualization purposes
void GAssemble(int size,double**M,double**C,double**K,double*G,double*forces)
{
    for(int i=0;i<size;i++)
    {
        G[i] = forces[i];
        for(int j=0;j<size;j++)
        {
            G[i] += ((2*M[i][j]/DTauSqd)-K[i][j])*CurDisplacement[j] + ((C[i][j]/TwiceDTau) - (M[i][j]/DTauSqd))*LastDisplacement[j];
        }
    }
}
void AAssemble(int size,double**M,double**C,double**A)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            A[i][j] = (M[i][j]/DTauSqd+C[i][j]/TwiceDTau);
        }
    }
}

void modExforce(double* exforce, int noofnodes, double curtau, Node*nodes, int noofjumps)
{
     double modcurtau = curtau;
     for(int i = 0; i < noofnodes*3; i++)
     {
          exforce[i] = 0;
     }
     for(int i = 0; i < noofnodes*3; i++)
     {
          if(i % 3 == 1)
          {
               exforce[i] -= nodes[(i-1)/3].getselffactor(2)*9.81;
          }
     }
     if(modcurtau < 1)
     {
          exforce[noofnodes-2] -= 700.0;
          return;
     }
     if(modcurtau < 1.1)
     {
          exforce[noofnodes-2] -= -5000*modcurtau + 5700;
          return;
     }
     if(modcurtau > (1.6 + noofjumps*1.5))
     {
         return;
     }
     if(modcurtau >= 1.1)
     {
          while(modcurtau > 2.6)
          {
               modcurtau -= 1.5;
          }
     }
     if(modcurtau < 1.4)
     {
          exforce[noofnodes-2] -= 5000*modcurtau - 5300;
          return;
     }
     if(modcurtau < 1.5)
     {
          exforce[noofnodes-2] -= -1000*modcurtau + 3100;
          return;
     }
     if(modcurtau < 1.6)
     {
          exforce[noofnodes-2] -= -16000*modcurtau + 25600;
          return;
     }
     if(modcurtau < 2.2)
     {
          return;
     }
     if(modcurtau < 2.3)
     {
          exforce[noofnodes-2] -= 20000*modcurtau - 44000;
          return;
     }
     if(modcurtau < 2.6)
     {
          exforce[noofnodes-2] -= -(1800/0.3)*modcurtau + 15800;
          return;
     }
     return;
}
int main()
{
    //Initialization.
    int NCNT, DoF = 3, numFixed = 4;
    int nodesfixed[4];
    nodesfixed[0] = 0;
    nodesfixed[1] = 1;
    nodesfixed[2] = 2;
    Node*nodes = NULL;
    bool grav = false;
    int size;
    ofstream output ("Output.txt");

    //set up timestep values
    cout << "Enter the timestep in seconds." << endl;
    cin >> DTau;
    cout << "Enter last time to solve for." << endl;
    cin >> TFinal;
    cout << "Enter 1 for gravity, 0 for no gravity"
        << endl << "(Gravity is assumed to act in the negative y direction)" << endl;
    cin >> grav;

    cout << "Enter the number of nodes you would like to use" << endl;
    cin >> NCNT;

    nodes = new Node[NCNT];
    //6 ribs, rib thickness 3mm, outer 6mm, 2.4384m half length, 0.034925m to 0.0508m to 0.022225m thick
    //0.498475m wide
    nodes[0].initProb(nodes,6,0.002,0.002,2.4384,2.4384,0.034925,0.0508,0.022225,0.498475,NCNT,2.4384,1,1,density,nodesfixed);
    size = DoF*NCNT;

    output << NCNT << endl;
    //allocate and zero out
    double *exforce = new double[size];
    CurDisplacement = new double[size];
    LastDisplacement = new double[size];
    NextDisplacement = new double[size];
    double*initialpos = new double[NCNT*2];
    double*areamoment = new double[NCNT-1];
    double*area = new double[NCNT-1];
    double *G = new double[size];
    double **A = CreateMatrix(size);
    double **K = CreateMatrix(size);
    double **C = CreateMatrix(size);
    double **M = CreateMatrix(size);
    double updatevector[DoF];
    double **submatrix = CreateMatrix(size-(numFixed));
    double *subG = new double[size-(numFixed)];
    double *subnextdis = new double[size-(numFixed)];
    double *velocity = new double[size];
    double *acceleration = new double[size];
    double *amplitude = new double[size];
    double *mean = new double[size];
    double *peak = new double[size];
    int ony = 0;

    double maxcomp = 0;
    double maxtens = 0;

    for(int i=0;i<size;i++)
    {
        CurDisplacement[i] = 0;
        LastDisplacement[i] = 0;
        NextDisplacement[i] = 0;
        velocity[i] = 0;
        acceleration[i] = 0;
        amplitude[i] = 0;
        mean[i] = 0;
        peak[i] = -1e20;
        exforce[i] = 0;
//        if(grav && (ony == 1))
//        {
//            exforce[i] -= 9.80665*nodes[(i/DoF)].getselffactor(2); //exploiting integer division here
//        }
//        if((ony++) >1)
//        {
//            ony = 0;
//        }
//        if(i == 28)
//        {
//            exforce[i] += -600;
//        }
    }

    int fixindx;
    int lj;

    //set up time
    long int cycles = ceil(TFinal/DTau)+1;
    double frames = ceil(TFinal*60);
    int plotinterval = cycles/frames; //plotting millions of points is insane

    if(plotinterval<1)
    {
        plotinterval = 1;
    }
    DTauSqd = DTau*DTau; //Precalculate repeatedly used values
    TwiceDTau = DTau*2;

    for(int i=0;i<NCNT;i++)
    {
        if(i>0)
        {
            areamoment[i-1] = nodes[i].GetAreaMoment();
            area[i-1] = nodes[i].GetConArea(0);
        }
        initialpos[i*2] = nodes[i].getInitialstate(0);
        initialpos[i*2+1] = nodes[i].getInitialstate(1);
        output << nodes[i].getInitialstate(0) << " "
               << nodes[i].getInitialstate(1) << " ";
    }
    output << endl;

    for (int i = 0;i<cycles;i++)
    {
        //now the coupled elements

        MAssemble(M,nodes,NCNT);
        KAssemble(K,nodes,NCNT);
        CAssemble(C,nodes,NCNT);

        //MPrint(K,size,size);
        //MPrint(C,size,size);
        //MPrint(M,size,size);
        fixindx = 0;

        modExforce(exforce,NCNT,i*DTau,nodes,2);
        GAssemble(size,M,C,K,G,exforce);
        AAssemble(size,M,C,A);

        //MPrint(A,size,size);
        //PrintV(G,size);

        createsubmatrix(A,submatrix,size,nodesfixed,numFixed);
        createsubG(G,subG,size,nodesfixed,numFixed);
        lud(submatrix,subG,(size-(numFixed)),subnextdis);

        //put next displacements into main vector
        fixindx = 0;
        lj = 0;

        for(int j=0;j<(size-(numFixed));j++)
        {
            while((nodesfixed[fixindx] == lj) && (fixindx < numFixed))
            {
                lj++;
                fixindx++;
                if(lj >= size)
                {
                    continue;
                }
            }
            NextDisplacement[lj] = subnextdis[j];
            lj++;
        }
        //update velocity and acceleration
        for(int j=0;j<size;j++)
        {
            velocity[j] = (NextDisplacement[j]-LastDisplacement[j])/TwiceDTau;
            acceleration[j] = ((NextDisplacement[j]-2*CurDisplacement[j]+LastDisplacement[j])/DTauSqd);
        }

        if((i%plotinterval) == 0)
        {
            for(int h=0;h<NCNT;h++)
            {
                for(int d=0;d<DoF;d++)
                    output << CurDisplacement[h*DoF + d] << " " ;
            }
            output << endl;
        }
        //Shift Indices
        for(int j=0;j<NCNT;j++)
        {
            for(int d=0;d<DoF;d++)
            {
                updatevector[d] = NextDisplacement[(DoF*j+d)];
                LastDisplacement[(j*DoF + d)] = CurDisplacement[(j*DoF + d)];
                CurDisplacement[(j*DoF + d)] = NextDisplacement[(j*DoF + d)];
            }
            nodes[j].UpdateDelta(updatevector);
        }



    CirMain* c_SSolve = new CirMain(NCNT,ModulusofElasticity);

	c_SSolve->Cin_PosIMtx(initialpos);
//	c_SSolve->Cout_PosIMtx();

	c_SSolve->Cin_DisMtx(CurDisplacement);
//	c_SSolve->Cout_DisMtx();

	c_SSolve->Cin_AreMtx(area);
//	c_SSolve->Cout_AreMtx();

    c_SSolve->Cin_AMIMtx(areamoment);
//	c_SSolve->Cout_AMIMtx();

	c_SSolve->Trans_PosIF();
//	c_SSolve->Cout_PosFMtx();

	c_SSolve->Trans_PosLg();
//	c_SSolve->Cout_PosLMtx();

	c_SSolve->Trans_PosAg();
//	c_SSolve->Cout_PosAMtx();

	c_SSolve->Trans_Multip();
//	c_SSolve->Cout_DisLMtx();

	c_SSolve->Stiff_MultipHut();
//	c_SSolve->Cout_ForLMtx();

//	c_SSolve->Cout_StsAMtx();
//	c_SSolve->Cout_StsBMtx();

	double tmptens = c_SSolve->Trans_StsTMtx();

	if(tmptens<0)
    {
        if(maxtens<(tmptens*-1))
        {
            maxtens = tmptens;
        }
    }
    else
    {
        if(maxtens<tmptens)
        {
            maxtens = tmptens;
        }
    }

	//c_SSolve->Cout_StsTMtx();

    delete c_SSolve;

    }
    cout << maxtens << endl;
    cout << "no seg faults" << endl;
    system("PAUSE");

    return 0;
}
