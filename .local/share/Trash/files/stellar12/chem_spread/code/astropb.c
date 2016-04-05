#include <stdio.h>
#include <math.h>
double function(double x)
{
	double y,Z,a,p,b,c,e;
	Z=0.02;
	a=0.004;
	p=0.00019;
	b=(Z*pow(10,x)+a)/p;
	c=Z*pow(10,x)/0.43429444819;
	e=exp((-pow(Z*pow(((10,x)+a),2),2)-a*a)/2*p);
	y=b*c*e;
	return y;
}

void main()
{
	double nrbins;
	double i, initialpt, finalpt, bsize, result;
	
	printf ("Please specify the number of bins ");
	scanf ("%e", &nrbins);
	printf ("At what value does the interval begin? ");
	scanf ("%e", &initialpt);
	printf ("At what value does the interval end? ");
	scanf ("%e", &finalpt);
	printf ("\n x                   F(x)\n");
	printf (" ---                  ------ \n");
	
	i=initialpt;
	bsize=(finalpt-initialpt)/nrbins;
	printf ("%g bsize  \n",bsize);
	printf ("%g finalpt  \n",finalpt);
	printf ("%g initialpt  \n",initialpt);
	printf ("%g nrbins  \n",nrbins);

	while (i<finalpt)
		{
		  /* result=function(i);
		 printf ("%g                   %g\n", i, result); */
		 i=i+bsize;
		}

}
	
