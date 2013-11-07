#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "matrixops.h"
#include "node.h"
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
double ModulusofElasticity = 70000000000; //
double ModulusofDamping;

struct ForceSinusoid
{
    public:
    double frequency;
    double amplitude;
};

struct Force
{
    public:
    int type; //0 for constant force, 1 for sinusoidal, 2 for ramp
    ForceSinusoid FS;
    double curval;
    double multiplier;
};

void ForcePrint(Force*array,int num)
{
    cout << endl;
    for(int i = 0;i<num;i++)
    {
        cout << array[i].curval << endl;
    }
}
//select chooses 0 (stiffness) 1 (dampness) or 2 (mass/moment)
void localwbending(double**k,int size, Node*list, int first, int sec, int select)
{
    //int number = i_node;
	double l = list[first].get_init_length(list[sec],0,1);
	double a = list[first].GetConArea(0);
	double i = list[first].GetAreaMoment();
	double c = list[first].getFrac(list[sec],0,1);
	double s = list[first].getFrac(list[sec],1,0);

	double al2 = a*l*l;
	double i12 = 12*i;
	double e13;
	if(select == 0)
    {
        e13 = ModulusofElasticity/l/l/l;
    }
    else if(select == 1)
    {
        e13 = ModulusofDamping/l/l/l;
    }
	//double el3 = d_ElaMod/l/l/l;

	double c2 = c*c;
	double s2 = s*s;
	double cs = c*s;

//	printf("%.3f, %.3f", c,s);
/*
	k =
	{
//		{			0,				1,			2,				3,					4,		5,	}
		{  i12*s2 +al2*c2,	 (al2-i12)*cs, -6*i*l*s, -i12*s2 -al2*c2,	  (i12-al2)*cs,  6*i*l*s},	// row 0
		{	 (al2-i12)*cs, i12*c2 +al2*s2,  6*i*l*c,	(i12-al2)*cs,  -i12*c2 -al2*s2,  6*i*l*c},	// row 1
		{		 -6*i*l*s,		  6*i*l*c,  4*i*l*l,		 6*i*l*s,		  -6*i*l*c,  2*i*l*l},	// row 2
		{ -i12*s2 -al2*c2,	 (i12-al2)*cs,	6*i*l*s,  i12*s2 +al2*c2,	  (al2-i12)*cs,  6*i*l*s},	// row 3
		{	 (i12-al2)*cs,-i12*c2 -al2*s2, -6*i*l*c,	(al2-i12)*cs,	i12*c2 +al2*s2, -6*i*l*c},	// row 4
		{		  6*i*l*s,		  6*i*l*c,	2*i*l*l,		 6*i*l*s,		  -6*i*l*c,  4*i*l*l},	// row 5
	};
*/
    k[0][0] = i12*s2 +al2*c2;
    k[0][1] = (al2-i12)*cs;
    k[0][2] = -6*i*l*s;
    k[0][3] = -i12*s2 -al2*c2;
    k[0][4] = (i12-al2)*cs;
    k[0][5] = 6*i*l*s;
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
    k[5][0] = 6*i*l*s;
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
    ////MPrint(k,size,size);
}

//Select is used to choose double derivative, derivative, or proportional (2,1, or 0)
//Note: This function has been basterdized for the diving board project and no longer works for general applications
void selfrefAssemble(Node*list,double**m,int select, int NCNT, int size, int DOF)
{
    zero(m,size);
    for(int i=0;i<NCNT;i++)
    {
        if(!(list[i].getselffactor(select) <= 1e-5))
        {
            for(int j=0; j<DOF;j++)
            {
                if(j != 2)
                    m[i*DOF+j][i*DOF+j]=list[i].getselffactor(select);
                else
                    m[i*DOF+j][i*DOF+j]=list[i].GetMoment();
            }
        }

    }
}

//Select is used to choose double derivative, derivative, or proportional (2,1, or 0)
//It is assumed that the dependant variable (ex. force) is only transferred axially
//transformation matrix is passed as a parameter so that it can be used multiple
//times without recalculating
void coupledAssemble(Node*list,double**c,int select, int NCNT, int DoF)
{
    zero(c,DoF*NCNT);
    double**transformation = CreateMatrix(DoF*2);
    if(select == 0)
    {
        for(int cnt=0;cnt<NCNT;cnt++)
        {
            for(int cnum=0;cnum<list[cnt].getconn();cnum++)
            {
                //list[cnt].gettransformation(list[list[cnt].connto(cnum)], transformation);
                localwbending(transformation,DoF*2,list,cnt,list[cnt].connto(cnum),select);
                MPrint(transformation,6,6);
                addLocToGlo(c,transformation,cnt,list[cnt].connto(cnum),DoF);
            }
        }
    }
    else if(select == 1)
    {
        for(int cnt=0;cnt<NCNT;cnt++)
        {
            for(int cnum=0;cnum<list[cnt].getconn();cnum++)
            {
                //list[cnt].gettransformation(list[list[cnt].connto(cnum)], transformation);
                list[cnt].formDampMatrix(transformation,CurDisplacement,list[list[cnt].connto(cnum)],cnt,0.497475,2);
                MPrint(transformation,6,6);
                addLocToGlo(c,transformation,cnt,list[cnt].connto(cnum),DoF);
            }
        }
    }
    DeleteMatrix(DoF*2,transformation);
}

// matrices named for M, C, and K for visualization purposes
void GAssemble(int size,double**M,double**C,double**K,double*G,Force*forces)
{
    for(int i=0;i<size;i++)
    {
        G[i] = forces[i].curval;
        for(int j=0;j<size;j++)
        {
            G[i] += ((-2*M[i][j]/DTauSqd)+K[i][j])*CurDisplacement[j] + (-(C[i][j]/TwiceDTau) + (M[i][j]/DTauSqd))*LastDisplacement[j];
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

int main()
{
    //Initialization.
    int NCNT, DoF = 3, numFixed = 3;
    int nodesfixed[3];
    nodesfixed[0] = 0;
    nodesfixed[1] = 1;
    Node*nodes = NULL;
    Force*exforce = NULL;
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
//    cout << "Enter the Modulus of Elasticity in GPa" << endl;
//    cin >> ModulusofElasticity ;
//    ModulusofElasticity *= 1e9;
    cout << "Enter the Modulus of Damping in whatever units" << endl;
    cin >> ModulusofDamping;
    cout << "Enter the number of nodes you would like to use" << endl;
    cin >> NCNT;

    nodes = new Node[NCNT];
    //6 ribs, rib thickness 3mm, outer 6mm, 2.4384m half length, 0.034925m to 0.0508m to 0.022225m thick
    //0.498475m wide
    nodes[0].initProb(nodes,6,0.006,0.003,2.4384,2.4384,0.034925,0.0508,0.022225,0.498475,NCNT,2.4384,1,1,2850,nodesfixed);
    size = DoF*NCNT;

    output << size << endl;
    //allocate and zero out
    exforce = new Force[size];
    CurDisplacement = new double[size];
    LastDisplacement = new double[size];
    NextDisplacement = new double[size];
    double *G = new double[size];
    double **A = CreateMatrix(size);
    double **Ks = CreateMatrix(size);
    double **Cs = CreateMatrix(size);
    double **Ms = CreateMatrix(size);
    double **Kc = CreateMatrix(size);
    double **Cc = CreateMatrix(size);
    double **Mc = CreateMatrix(size);
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
        exforce[i].curval = 0;
        if(grav && ((ony == 1)))
        {
            exforce[i].curval -= 9.80665*nodes[(i/DoF)].getselffactor(2); //exploiting integer division here
        }
        if((ony++) >1)
        {
            ony = 0;
        }
    }

    //first the selfreferential elements
    //they don't change
    //selfrefAssemble(nodes,Ks,0,NCNT,size,DoF);
    //selfrefAssemble(nodes,Cs,1,NCNT,size,DoF);
    //selfrefAssemble(nodes,Ms,2,NCNT,size,DoF);

    int fixindx;
    int lj;

    //set up time
    int cycles = ceil(TFinal/DTau)+1;
    int plotinterval = cycles/5000; //plotting millions of points is insane
    if(plotinterval<1)
    {
        plotinterval = 1;
    }
    DTauSqd = DTau*DTau; //Precalculate repeatedly used values
    TwiceDTau = DTau*2;
    for (int i = 0;i<cycles;i++)
    {
        //now the coupled elements

        if(1) //this should be replaced by criteria for updating the matrices
        {
            coupledAssemble(nodes,Kc,0,NCNT,DoF);
            coupledAssemble(nodes,Cc,1,NCNT,DoF);
            //addm(Kc,Ks,Kc,size,size);
            //addm(Cc,Cs,Cc,size,size);
            //coupledAssemble(nodes,Mc,2,NCNT,DoF);
        }

        MPrint(Kc,size,size);
        MPrint(Cc,size,size);
        MPrint(Ms,size,size);
        fixindx = 0;

        GAssemble(size,Ms,Cc,Kc,G,exforce);
        AAssemble(size,Ms,Cc,A);

        MPrint(A,size,size);
        PrintV(G,size);

        createsubmatrix(A, submatrix,size,nodesfixed,numFixed);
        createsubG(G,subG,size,nodesfixed,numFixed);

        MPrint(submatrix,(size-numFixed),(size-numFixed));

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

        //code to see if system if stable.
        //it's not (very bad)
        for(int d=0;d<size;d++)
        {
            exforce[d].curval = 0;
        }
    }

    cout << "no seg faults" << endl;
    cin.get();
    return 0;
}
