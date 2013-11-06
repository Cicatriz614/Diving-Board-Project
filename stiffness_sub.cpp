void CirMain::Stiff_Single(int i_node, double d_Length, double d_CrArea, double d_MomInt, double d_GloAgl) // angle in rads
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
			d_ForMtx[number*6+f] += el3*stiffarray[f][g]*d_DisMtx[number*6+g];
		}

	}
}
