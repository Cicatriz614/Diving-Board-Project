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
double ModulusofElasticity = 70000000000;
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
    double cosine = list[first].getFrac(list[list[first].connto(sec)],0,1);
    double sine = list[first].getFrac(list[list[first].connto(sec)],1,0);
    double length = list[first].get_2D_Distance(list[list[first].connto(sec)],0,1); //can hopefully replace this with something more efficient since node are all adjacent to one another
    double factor = 0;

    if(select == 0)
    {
        factor = ModulusofElasticity;
    }
    else if(select == 1)
    {
        factor = ModulusofDamping;
    }
    int compare;
    for(int i = 0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            compare = i%3 + j%3;
            k[i][j] = factor/(length*length*length);
            switch(compare)
            {
            case 0:
                k[i][j] *= 12*(list[first].GetMoment())*sine*sine + (list[first].GetConArea(sec))*length*length*cosine*cosine;
                break;
            case 1:
                k[i][j] *= cosine*sine*((list[first].GetConArea(sec))*length*length - 12*(list[first].GetMoment()));
                break;
            case 2:
                if(i%3 != j%3)
                {
                    k[i][j] *= -6*(list[first].GetMoment())*length*sine;
                }
                else
                {
                    k[i][j] *= -6*(list[first].GetMoment())*length*sine;
                }
            case 3:
                k[i][j] *= 6*(list[first].GetMoment())*length*cosine;
            case 4:
                k[i][j] *= 4*(list[first].GetMoment())*length*length;
            }
            if((i>2 && j<2) || (i<2 && j>2)) //top right or bottom left quadrants of matrix
            {
                k[i][j] *= -1;
            }
        }
    }
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
    for(int cnt=0;cnt<NCNT;cnt++)
    {
        for(int cnum=0;cnum<list[cnt].getconn();cnum++)
        {
            //list[cnt].gettransformation(list[list[cnt].connto(cnum)], transformation);
            localwbending(transformation,DoF*2,list,cnt,cnum,select);
            addLocToGlo(c,transformation,cnt,list[cnt].connto(cnum),DoF,list[cnt].getconnfactor(cnum,select));
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

int main()
{
    //Initialization.
    int NCNT, DoF = 3, numFixed = 0;
    int *nodesfixed;
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

    //allocate and zero out
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
    double **submatrix = CreateMatrix(size-(DoF*numFixed));
    double *subG = new double[size-(DoF*numFixed)];
    double *subnextdis = new double[size-(DoF*numFixed)];
    double *velocity = new double[size];
    double *acceleration = new double[size];
    double *amplitude = new double[size];
    double *mean = new double[size];
    double *peak = new double[size];

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
        if(grav && ((i&1) == 1))
        {
            exforce[i].curval -= 9.80665*nodes[(i/DoF)].getselffactor(2); //exploiting integer division here
        }
    }

    //first the selfreferential elements
    //they don't change
    selfrefAssemble(nodes,Ks,0,NCNT,size,DoF);
    selfrefAssemble(nodes,Cs,1,NCNT,size,DoF);
    selfrefAssemble(nodes,Ms,2,NCNT,size,DoF);

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
            coupledAssemble(nodes,Mc,2,NCNT,DoF);
        }

        addm(Kc,Ks,Kc,size,size);
        addm(Cc,Cs,Cc,size,size);
        addm(Mc,Ms,Mc,size,size);

        fixindx = 0;
        for(int j=0;j<size;j++)
        {
            if(exforce[j].type != 0) //if the force is external varying
            {
                if(exforce[j].type == 1) //sinusoid
                {
                    if(j&1)
                    {
                        exforce[j].curval = (exforce[j].FS.amplitude)*sin((exforce[j].FS.frequency)*((i)*DTau)) - grav*(9.80665*nodes[(j/DoF)].getselffactor(2));
                    }
                    else
                    {
                        exforce[j].curval = (exforce[j].FS.amplitude)*sin((exforce[j].FS.frequency)*((i)*DTau));
                    }
                }
                else //ramp
                {
                    if(j&1)
                    {
                        exforce[j].curval = (exforce[j].multiplier)*((i)*DTau) - grav*(9.80665*nodes[(j/DoF)].getselffactor(2));
                    }
                    else
                    {
                        exforce[j].curval = (exforce[j].multiplier)*((i)*DTau);
                    }
                }
            }
        }
        GAssemble(size,Mc,Cc,Kc,G,exforce);
        AAssemble(size,Mc,Cc,A);
        createsubmatrix(A, submatrix,size,nodesfixed);
        createsubG(G,subG,size,nodesfixed);
        lud(submatrix,subG,(size-(DoF*numFixed)),subnextdis);

        //put next displacements into main vector
        fixindx = 0;
        lj = 0;

        for(int j=0;j<(size-(DoF*numFixed));j++)
        {
            while(nodesfixed[fixindx] == lj)
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
    }

    cout << "no seg faults" << endl;
    cin.get();
    return 0;
}
