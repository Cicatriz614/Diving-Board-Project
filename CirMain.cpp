

#include "stdafx.h"
#include "CirMain.h"
///*

#include <algorithm> // for std::swap
#include <cstddef>
#include <cassert>
 
////////////////////////////////////////////////////////////////////////////////////
 ////////////begin code borrowed

// Matrix traits: This describes how a matrix is accessed. By
// externalizing this information into a traits class, the same code
// can be used both with native arrays and matrix classes. To use the
// dfault implementation of the traits class, a matrix type has to
// provide the following definitions as members:
//
// * typedef ... index_type;
//   - The type used for indexing (e.g. size_t)
// * typedef ... value_type;
//   - The element type of the matrix (e.g. double)
// * index_type min_row() const;
//   - returns the minimal allowed row index
// * index_type max_row() const;
//   - returns the maximal allowed row index
// * index_type min_column() const;
//   - returns the minimal allowed column index
// * index_type max_column() const;
//   - returns the maximal allowed column index
// * value_type& operator()(index_type i, index_type k)
//   - returns a reference to the element i,k, where
//     min_row() <= i <= max_row()
//     min_column() <= k <= max_column()
// * value_type operator()(index_type i, index_type k) const
//   - returns the value of element i,k
//
// Note that the functions are all inline and simple, so the compiler
// should completely optimize them away.
/*
template<typename MatrixType> 

struct matrix_traits
{
  typedef typename MatrixType::index_type index_type;
  typedef typename MatrixType::value_typ value_type;
  static index_type min_row(MatrixType const& A)
  { return A.min_row(); }
  static index_type max_row(MatrixType const& A)
  { return A.max_row(); }
  static index_type min_column(MatrixType const& A)
  { return A.min_column(); }
  static index_type max_column(MatrixType const& A)
  { return A.max_column(); }
  static value_type& element(MatrixType& A, index_type i, index_type k)
  { return A(i,k); }
  static value_type element(MatrixType const& A, index_type i, index_type k)
  { return A(i,k); }
};
 
// specialization of the matrix traits for built-in two-dimensional arrays
template<typename T, std::size_t rows, std::size_t columns>
 struct matrix_traits<T[rows][columns]>
{
  typedef std::size_t index_type;
  typedef T value_type;
  static index_type min_row(T const (&)[rows][columns])
  { return 0; }
  static index_type max_row(T const (&)[rows][columns])
  { return rows-1; }
  static index_type min_column(T const (&)[rows][columns])
  { return 0; }
  static index_type max_column(T const (&)[rows][columns])
  { return columns-1; }
  static value_type& element(T (&A)[rows][columns],
                             index_type i, index_type k)
  { return A[i][k]; }
  static value_type element(T const (&A)[rows][columns],
                            index_type i, index_type k)
  { return A[i][k]; }
};
 
// Swap rows i and k of a matrix A
// Note that due to the reference, both dimensions are preserved for
// built-in arrays
template<typename MatrixType>
 void swap_rows(MatrixType& A,
                 typename matrix_traits<MatrixType>::index_type i,
                 typename matrix_traits<MatrixType>::index_type k)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;
 
  // check indices
  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));
 
  assert(mt.min_row(A) <= k);
  assert(k <= mt.max_row(A));
 
  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    std::swap(mt.element(A, i, col), mt.element(A, k, col));
}
 
// divide row i of matrix A by v
template<typename MatrixType>
 void divide_row(MatrixType& A,
                  typename matrix_traits<MatrixType>::index_type i,
                  typename matrix_traits<MatrixType>::value_type v)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;
 
  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));
 
  assert(v != 0);
 
  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    mt.element(A, i, col) /= v;
}
 
// in matrix A, add v times row k to row i
template<typename MatrixType>
 void add_multiple_row(MatrixType& A,
                  typename matrix_traits<MatrixType>::index_type i,
                  typename matrix_traits<MatrixType>::index_type k,
                  typename matrix_traits<MatrixType>::value_type v)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;
 
  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));
 
  assert(mt.min_row(A) <= k);
  assert(k <= mt.max_row(A));
 
  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    mt.element(A, i, col) += v * mt.element(A, k, col);
}
 
// convert A to reduced row echelon form
template<typename MatrixType>
 void to_reduced_row_echelon_form(MatrixType& A)
{
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;
 
  index_type lead = mt.min_row(A);
 
  for (index_type row = mt.min_row(A); row <= mt.max_row(A); ++row)
  {
    if (lead > mt.max_column(A))
      return;
    index_type i = row;
    while (mt.element(A, i, lead) == 0)
    {
      ++i;
      if (i > mt.max_row(A))
      {
        i = row;
        ++lead;
        if (lead > mt.max_column(A))
          return;
      }
    }
    swap_rows(A, i, row);
    divide_row(A, row, mt.element(A, row, lead));
    for (i = mt.min_row(A); i <= mt.max_row(A); ++i)
    {
      if (i != row)
        add_multiple_row(A, i, row, -mt.element(A, i, lead));
    }
  }
}
 
 */
 ////////////end code borrowed
 //////////////////////////////////////////////////////////////////////////////




CirMain::CirMain(void)
{
	i_nodenm = 0;
	d_ElaMod = 0;
}

CirMain::CirMain(int nodes, double ElasticMod)
{
	i_nodenm = nodes;
	d_ElaMod = ElasticMod;

	d_PosIMtx = (double*)calloc(nodes*2, sizeof (double));	//Initial Node Postions, in x, y, with respect to First Node
	d_PosFMtx = (double*)calloc(nodes*2, sizeof (double));	//Initial Node Postions, in x, y, with respect to First Node
	
	d_PosLIMtx = (double*)calloc(nodes*1, sizeof (double));	//Initial Element Length
	d_PosLFMtx = (double*)calloc(nodes*1, sizeof (double));	//Final Element Length
	d_PosAMtx = (double*)calloc(nodes*1, sizeof (double));	//Strain Element Length, d_PosLFMtx - d_PosLIMtx

	d_DisMtx = (double*)calloc(nodes*3, sizeof (double));	//Global Node Displacement, in x, y, w, with respect to Respective Node
	d_ForMtx = (double*)calloc(nodes*3, sizeof (double));	//Initial Node Postions, with respect to First Node

	d_DisLMtx = (double*)calloc(nodes*6, sizeof (double));	//Local Element Displacement, in x1, y1, w1, x2, y2, w2, with respect to Respective Node
	d_ForLMtx = (double*)calloc(nodes*6, sizeof (double));	//Local Element Forces, in u1, v1, o1, u2, v2, o2, with respect to Respective Node

	d_StsAMtx = (double*)calloc(nodes*1, sizeof (double));	//Local Element Axial Stress
	d_StsBMtx = (double*)calloc(nodes*1, sizeof (double));	//Local Element Bending
	d_StsTMtx = (double*)calloc(nodes*2, sizeof (double));	//Local Element Bending
}


CirMain::~CirMain(void)
{
	delete d_DisMtx;
	delete d_ForMtx;
}


void CirMain::Trans_PosIF() // angle in rads
{	
	for (int i = 0; i < i_nodenm; i++)
	{		
		d_PosFMtx[2*i+0] = d_PosIMtx[2*i+0]+ d_DisMtx[3*i+0];
		d_PosFMtx[2*i+1] = d_PosIMtx[2*i+1]+ d_DisMtx[3*i+1];
	}
}
void CirMain::Trans_PosLg() // angle in rads
{	
	double dx = 1;
	double dy = 1;
	double ddx = 1;
	double ddy = 1;

	for (int i = 0; i < i_nodenm-1; i++)
	{		
		dx = (d_PosIMtx[2*(i+1) +0] - d_PosIMtx[2*(i+0) +0]);
		dy = (d_PosIMtx[2*(i+1) +1] - d_PosIMtx[2*(i+0) +1]);

		d_PosLIMtx[i] = sqrt(dy*dy+dx*dx);

		ddx = (d_PosFMtx[2*(i+1) +0] - d_PosFMtx[2*(i+0) +0]);
		ddy = (d_PosFMtx[2*(i+1) +1] - d_PosFMtx[2*(i+0) +1]);

		d_PosLFMtx[i] = sqrt(ddy*ddy+ddx*ddx);
	}
}
void CirMain::Trans_PosAg() // angle in rads
{	
	double dx = 1;
	double dy = 1;
	for (int i = 0; i < i_nodenm-1; i++)
	{		
		dx = (d_PosFMtx[2*(i+1) +0] - d_PosFMtx[2*(i+0) +0]);
		dy = (d_PosFMtx[2*(i+1) +1] - d_PosFMtx[2*(i+0) +1]);

		d_PosAMtx[i] = atan(dy/dx);
	}
}


void CirMain::Stiff_Multip204() // angle in rads
{
	for (int i = 0; i < i_nodenm-1; i++)
	{
		Stiff_Single204(i, d_PosLFMtx[i], 0.01, 0.0000000001, d_PosAMtx[i]);
	}
}
void CirMain::Stiff_Single204(int i_node, double d_Length, double d_CrArea, double d_MomInt, double d_GloAgl) // angle in rads
{
	int number = i_node;
	double l = d_Length;
	double a = d_CrArea;
	double i = d_MomInt;
	double c = cos(d_GloAgl);
	double s = sin(d_GloAgl);

	double al2 = a*l*l;
	double i12 = 12*i;
	double el3 = d_ElaMod/l/l/l;
	
	double c2 = c*c;
	double s2 = s*s;
	double cs = c*s;

//	printf("%.3f, %.3f", c,s);

	double stiffarray[6][6] =
	{
//		{			1,				2,			3,				4,					5,		6,	}
		{  i12*s2 +al2*c2,	 (al2-i12)*cs, -6*i*l*s, -i12*s2 -al2*c2,	  (i12-al2)*cs, -6*i*l*s},	// row 0
		{	 (al2-i12)*cs, i12*c2 +al2*s2,  6*i*l*c,	(i12-al2)*cs,  -i12*c2 -al2*s2,  6*i*l*c},	// row 1
		{		 -6*i*l*s,		  6*i*l*c,  4*i*l*l,		 6*i*l*s,		  -6*i*l*c,  2*i*l*l},	// row 2
		{ -i12*s2 -al2*c2,	 (i12-al2)*cs,	6*i*l*s,  i12*s2 +al2*c2,	  (al2-i12)*cs,  6*i*l*s},	// row 3
		{	 (i12-al2)*cs,-i12*c2 -al2*s2, -6*i*l*c,	(al2-i12)*cs,	i12*c2 +al2*s2, -6*i*l*c},	// row 4
		{		 -6*i*l*s,		  6*i*l*c,	2*i*l*l,		 6*i*l*s,		  -6*i*l*c,  4*i*l*l},	// row 5
	};
	
	for (int f = 0; f < 6; f++)
	{
		for (int g = 0; g < 6; g++)
		{
			d_ForMtx[number*6+f] += el3*stiffarray[f][g]*d_DisMtx[number*3+g];
		}
	}
}


void CirMain::Stiff_MultipHut() // angle in rads
{
	for (int i = 0; i < i_nodenm-1; i++)
	{
		Stiff_SingleHut(i, d_PosLFMtx[i], 0.01, 0.00000001, d_PosAMtx[i]);
	}
}
void CirMain::Stiff_SingleHut(int i_node, double d_Length, double d_CrArea, double d_MomInt, double d_GloAgl) // angle in rads
{
	int number = i_node;
	double l = d_Length;
	double a = d_CrArea;
	double i = d_MomInt;
//	double c = cos(d_GloAgl);
//	double s = sin(d_GloAgl);

	double al2 = a*l*l;
	double i12 = 12*i;
	double el3 = d_ElaMod/l/l/l;

//	double c2 = c*c;
//	double s2 = s*s;
//	double cs = c*s;

//	printf("%.3f, %.3f", c,s);

	double stiffarray[6][6] =
	{
//		{	1,		2,		3,		4,		5,		6}
		{  al2,		0,		0,	 -al2,		0,		0},	// row 0
		{	 0,	  i12,  6*i*l,		0,	 -i12,  6*i*l},	// row 1
		{	 0,	6*i*l,4*i*l*l,		0, -6*i*l,2*i*l*l},	// row 0
		{ -al2,		0,		0,	  al2,		0,		0},	// row 0
		{	 0,	 -i12, -6*i*l,		0,	  i12, -6*i*l},	// row 1
		{	 0,	6*i*l,2*i*l*l,		0, -6*i*l,4*i*l*l},	// row 0
	};
	
	for (int f = 0; f < 6; f++)
	{
		for (int g = 0; g < 6; g++)
		{
			d_ForLMtx[number*6+f] += el3*stiffarray[f][g]*d_DisLMtx[number*6+g];
		}
	}
}


void CirMain::Trans_Multip() // angle in rads
{
	for (int i = 0; i < i_nodenm-1; i++)
	{
		Trans_Single(i, d_PosAMtx[i]);
	}
}
void CirMain::Trans_Single(int i_node, double d_GloAgl) // angle in rads
{
	int number = i_node;
	
	double c = cos(d_GloAgl);
	double s = sin(d_GloAgl);

	double transarray[6][6] =
	{
		{		c,		s,		0,		0,		0,		0},	// row 0
		{	   -s,		c,		0,		0,		0,		0},	// row 1		
		{		0,		0,		1,		0,		0,		0},	// row 2
		{		0,		0,		0,		c,		s,		0},	// row 3
		{		0,		0,		0,	   -s,		c,		0},	// row 4
		{		0,		0,		0,		0,		0,		1},	// row 5
	};
	
	for (int f = 0; f < 6; f++)
	{
		for (int g = 0; g < 6; g++)
		{
			d_DisLMtx[number*6+f] += transarray[f][g]*d_DisMtx[number*3+g];
		}

	}
}
void CirMain::Trans_StsTMtx() // angle in rads
{
	for (int i = 0; i < i_nodenm-1; i++)
	{
		d_StsTMtx[i*2] =   -d_ForLMtx[i*6]/0.01	-d_ForLMtx[i*6+2]*0.01/0.00000001;
		d_StsTMtx[i*2+1] = -d_ForLMtx[i*6]/0.01	+d_ForLMtx[i*6+2]*0.01/0.00000001;
	}
}


void CirMain::Cin_PosIMtx()
{
	for (int i = 0; i < i_nodenm; i++)
	{
	//	d_PosIMtx[2*i+0] = 0;
	//	d_PosIMtx[2*i+1] = 2*i;
		d_PosIMtx[2*i+0] = 0.1*i;
		d_PosIMtx[2*i+1] = 0;
	}
}
void CirMain::Cin_DisMtx()
{
//	for (int i = 0; i < i_nodenm; i++)
//	{
//		d_DisMtx[3*i+0] = 0;
//		d_DisMtx[3*i+1] = i;
//	}
	///*
	for (int i = 0; i < i_nodenm; i++)
	{
		for (int j = 0; j < 3; j++)
		{	
			switch (j)
			{
				case 0:	
					d_DisMtx[i*3+j] = 0.001*i;  break;
				case 1:	
					d_DisMtx[i*3+j] = 0.001*i*i;  break;
				case 2:	
					d_DisMtx[i*3+j] = 0.02;  break;
			}
		}
	}
//	*/
}


void CirMain::Cout_PosIMtx()
{
	printf("Position Initial Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("section %d", i);
		printf("\n");

		for (int j = 0; j < 2; j++)
		{	
			printf("|	%.4f	|", d_PosIMtx[i*2+j]);
			switch (j)
			{
				case 0:	
					printf("Px");  break;
				case 1:	
					printf("Py");  break;
			}
			printf("\n");
		}
	}
}
void CirMain::Cout_PosFMtx()
{
	printf("Position Final Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("node %d", i);
		printf("\n");

		for (int j = 0; j < 2; j++)
		{	
			printf("|	%.4f	|", d_PosFMtx[i*2+j]);
			switch (j)
			{
				case 0:	
					printf("Px");  break;
				case 1:	
					printf("Py");  break;
			}
			printf("\n");
		}
	}
}
void CirMain::Cout_PosAMtx()
{
	printf("Angle Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("section %d", i);
		printf("|	%.4f	|", d_PosAMtx[i]);
		printf("\n");
	}
}
void CirMain::Cout_PosLMtx()
{
	printf("Length Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("element %d", i);
		printf("|	%.4f	|" "|	%.4f	|",  d_PosLIMtx[i], d_PosLFMtx[i]);
		printf("\n");
	}
}
void CirMain::Cout_DisMtx()
{
	printf("Global Displacement - Global Frame Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("node %d", i);
		printf("\n");

		for (int j = 0; j < 3; j++)
		{	
			printf("|	%.4f	|", d_DisMtx[i*3+j]);
			
			switch (j)
			{
				case 0:	
					printf("Ux");  break;
				case 1:	
					printf("Uy");  break;
				case 2:	
					printf("W");  break;
			}
			printf("\n");
		}
//		printf("\n");
	}
}

void CirMain::Cout_ForMtx()
{
	printf("Local Force - Global Frame Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("element %d", i);
		printf("\n");

		for (int j = 0; j < 6; j++)
		{	
			printf("|	%.0f	|", d_ForMtx[i*6+j]);
			
			switch (j)
			{
				case 0:	
					printf("Fix");  break;
				case 1:	
					printf("Fiy");  break;
				case 2:	
					printf("Mi");  break;
				case 3:	
					printf("Fjx");  break;
				case 4:	
					printf("Fjx");  break;
				case 5:	
					printf("Mj");  break;
			}
			printf("\n");
		}
//		printf("\n");
	}
}



void CirMain::Cout_DisLMtx()
{
	printf("Local Displacement - Local Frame Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("element %d", i);
		printf("\n");

		for (int j = 0; j < 6; j++)
		{	
			printf("|	%.5f	|", d_DisLMtx[i*6+j]);
			
			switch (j)
			{
				case 0:	
					printf("uix");  break;
				case 1:	
					printf("uiy");  break;
				case 2:	
					printf("oi");  break;
				case 3:	
					printf("ujx");  break;
				case 4:	
					printf("ujy");  break;
				case 5:	
					printf("oj");  break;
			}
			printf("\n");
		}
//		printf("\n");
	}
}
void CirMain::Cout_ForLMtx()
{
	printf("Local Force - Local Frame Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("element %d", i);
		printf("\n");

		for (int j = 0; j < 6; j++)
		{	
			printf("|	%.0f	|", d_ForLMtx[i*6+j]);
			
			switch (j)
			{
				case 0:	
					printf("uix");  break;
				case 1:	
					printf("uiy");  break;
				case 2:	
					printf("oi");  break;
				case 3:	
					printf("ujx");  break;
				case 4:	
					printf("ujy");  break;
				case 5:	
					printf("oj");  break;
			}
			printf("\n");
		}
//		printf("\n");
	}
}

void CirMain::Cout_StsAMtx()
{
	printf("Axial Stress Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("element %d", i);
		printf("|	%.0f	|" "|	%.0f	|",  d_ForLMtx[i*6]/0.01,   (d_PosLFMtx[i]-d_PosLIMtx[i])/d_PosLIMtx[i]*d_ElaMod);
		printf("\n");
	}
}
void CirMain::Cout_StsBMtx()
{
	printf("Bending Stress Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("element %d", i);
		printf("|	%.0f	|" "|	%.0f	|",  -d_ForLMtx[i*6+2]*0.01/0.00000001,  
			0.01*d_ElaMod*((6*(d_DisLMtx[i*6+4]-d_DisLMtx[i*6+1])/d_PosLFMtx[i]/d_PosLFMtx[i])-
			(2*(2*d_DisLMtx[i*6+5]+d_DisLMtx[i*6+2])/d_PosLFMtx[i]))
			);
		printf("\n");
	}
}

///*
void CirMain::Cout_StsTMtx()
{
	printf("Bending Stress Matrix \n");	
	for (int i = 0; i < i_nodenm; i++)
	{
		printf("element %d", i);
		printf("|	%.0f	|" "|	%.0f	|",  d_StsTMtx[i*2], d_StsTMtx[i*2+1]);
		printf("\n");
	}
}
//*/
/*
void CirMain::Cout_StsTMtx()
{

}
*/

//	Engine *p_Engine;
//	p_Engine = engOpen("null");
/*
  reduced_row_echelon_form(M);
  
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 4; ++j)
      std::cout << M[i][j] << '\t';
    std::cout << "\n";
  }
*/

//*/
