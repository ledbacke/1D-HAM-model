#pragma once
#ifndef MATRIX_H
#define MATRIX_H

double *GaussSeidell(const double subd[], const double diag[], double superd[], double b[],int size)
{
    //1. LDU-ontbinding
    //--------------------------------------------------------------------------------------------------------
    const int n = size;   // get the length of the diagonal
 
    //initialize the output vectors
    double *u   = new double [n-1];     // superd   (from 0 to n-2);
    double *d   = new double [n];       // diag     (from 0 to n-1);     
    double *l   = new double [n-1];     // subd     (from 1 to n-1);
    for(int i=0;i<n-1;i++) u[i] = superd[i];
    for(int i=0;i<n;i++)   d[i] = diag[i];
    for(int i=0;i<n-1;i++) l[i] = subd[i+1];
        
    //Perform LU decomposition
    d[0] = diag[0];
    u[0] = superd[0];
    for (int i=1;i<n;i++)
    {
        l[i-1]  = superd[i]/d[i-1];                           // Update the lover triangle vector
        d[i]    = diag[i] - (subd[i]/d[i-1])*superd[i];       // Update the diaginal  
    }
    
    // 2 . SOLVE X
    double *x = new double [n];
    double *y = new double [n];

    //Solve Tridiagonal system LUx=b (Forward);

    // Step 1 : Solve Ly=b for y
    y[0] = b[0];
    for (int i=1; i<n; i++)        y[i] = b[i] - l[i-1]*y[i-1];
    
    // Step 2 : Solve Ux=y for x (Backward)
    x[n-1] = y[n-1]/d[n-1];  
    for (int i=(n-2);i>=0;i--)     x[i] = (y[i]-u[i]*x[i+1])/d[i];
    
    return (x);
}
    
    
double *Thomas_algorithm(const double subd[], const double diag[], double superd[], double q[],int size) 
{
          
        /*
         solves Ax = d where A is a tridiagonal matrix consisting of vectors a, b, c
         x[] - initially contains the input vector y, and returns the solution x.  
         n   - number of equations
         a[] - main diagonal        indexed from [1, ..., n - 1]
         b[] - subdiagonal          indexed from [0, ..., n - 1]
         c[] - superdiagonal        indexed from [0, ..., n - 2]
		 q[] - output				indexed from [0, ..., n - 1]
         */
        const int n = size;
		static double *l,*u,*d,*y, *xnew;
		l		= new double [n];
		u		= new double [n];
		d		= new double [n];
		y		= new double [n];
		xnew	= new double [n];

		// A = LU  
		//STAP 1: LU Decomposition
		d[0] = diag[0];
		u[0] = superd[0];
		l[0] = 0;  //bestaat niet
		for (int i = 1; i < n-1; i++) 
        {
			l[i]	= subd[i]/d[i-1];
			u[i]	= superd[i];
			d[i]	= diag[i]-l[i]*u[i-1];
			
        }
		l[n-1] = subd[n-1]/d[n-2];
		d[n-1] = diag[n-1]-l[n-1]*u[n-2];
		u[n-1] = 0;  //bestaat niet

		//STAP 2:forward subs: Ly=q => y
		y[0] = q[0];
		for(int i=1;i<n;i++)	y[i]=q[i]-l[i]*y[i-1];
		
		//Backward Ux = y =>x
		xnew[n-1] = y[n-1]/d[n-1];
        for (int i = (n-2); i>= 0;i--)	xnew[i] = (y[i]-u[i]*xnew[i + 1])/d[i];
		
		return &xnew[0];
}


double *GSOR(const double subd[], const double diag[], double superd[],double x_old[], double b[],double &omega,int size)
{
    //1. LDU-ontbinding
    //--------------------------------------------------------------------------------------------------------
    const int n = size;   // get the length of the diagonal
 
    //initialize the output vectors
    double *u   = new double [n];   // superd;
    double *d   = new double [n];   // diag
    double *l   = new double [n];   // superd
     
    //Perform LU decomposition
    d[0] = diag[0];
    for (int i=1;i<n;i++)
    {
        l[i-1]  = superd[i-1]/d[i-1];                           // Update the lover triangle vector
        d[i]    = diag[i] - (subd[i-1]/d[i-1])*superd[i-1];   // Update the diaginal  
    }
    
    // 2 . SOLVE X
    double *x = new double [n];
    double *y = new double [n];

    //Solve Tridiagonal system LUx=b (Forward);

    // Step 1 : Solve Ly=b for y
    y[0] = b[0];
    for (int i=1; i<n; i++)        y[i] = b[i] - l[i-1]*y[i-1];
    
    // Step 2 : Solve Ux=y for x (Backward)
    x[n-1] = y[n-1]/d[n-1];  
    for (int i=(n-2);i>=0;i--)     x[i] = (y[i]-u[i]*x[i+1])/d[i];

	//Relaxation
    double sum1	=	0.0;
	double sum2	=	0.0;
    double *x_new = new double[n];
	
	double **A = new double *[n];
	for(int i=0; i<n; i++)
    {
		A[i] = new double[n];
	}
	for(int i=0;i<n;i++)
    {
		for(int j=0; j<n; j++)			A[i][j]=0.0;
	}
	for(int rij = 0; rij<n; rij++)
    {
		if(rij == 0)
        {
			A[0][0] = diag[0];
			A[0][1] = superd[0];
		}
		else if(rij == n-1)
        {
			A[n-1][n-1] = diag[n-1];
			A[n-1][n-2] = subd[n-1];
		}
		else
        {
			for(int kolom=0 ; kolom<(rij-1);kolom++)  A[rij][kolom] = 0;
			
			A[rij][rij-1] = subd[rij];
			A[rij][rij]   = diag[rij];
			A[rij][rij+1] = superd[rij];
			for(int kolom=rij+2; kolom<=rij; kolom++) A[rij][kolom]=0;
		}
	}
	///*
	for(int i=0;i<n;i++)    sum1 = sum1 + (x[i]-x_new[i]);
	if(sum1>2)              omega = 2/sum1;
    else                    omega = 0.7;
	
	//*/
	for(int i=0;i<n;i++)
    {
		sum1 = 0.0;
		sum2 = 0.0;
		for(int j=0; j<i;j++)			sum1 = sum1+A[i][j]*x_new[j];
		
		for(int j=i+1;j<n;j++)			sum2 = sum2 + A[i][j]*x[j];
		
		x_new[i] = (1-omega)*(x[i]) + (omega)*(-sum1 -sum2 + b[i])/A[i][i];

	}
	return &x_new[0];
};

//Aitken Interface Underrelaxation uses a scalar value to underrelax the entire vector of the interface variables
	double *relaxation_Aitken(double mc[], double mc_p[], double mc_pp[], double & om_p, int aantal_midcellen)
    {
	//m_c				new estimated value
	//m_c_p				last value
	//m_c_pp			second last value
	//om_p				previous value of the relaxation factor
	

	//Aitken Interface Underrelaxation uses a scalar value to underrelax the entire vector of the interface variables
	double om;
	double *delta		= new double [aantal_midcellen];
	double *delta_prev	= new double [aantal_midcellen];
	double relax_f_maxspeed = 1.5;
	double relax_f_minslow	= 0.5;
	double Mc	 = 0.;
	double Mc_p	 = 0.;
	double Mc_pp = 0.;
	double Delta = 0.;
	double Delta_prev = 0.;
	
	double *x_new = new double[aantal_midcellen];
	x_new = mc;

	//Calculate Residual
	for(int i =0; i<aantal_midcellen; i++)
	{
		delta[i]		= mc[i] - mc_p[i] ;
		delta_prev[i]	= mc_p[i] - mc_pp[i];

		Delta       = fabs(delta[i]) + Delta;
		Delta_prev  = fabs(delta_prev[i]) + Delta_prev;
		Mc          = mc[i]    + Mc;
		Mc_p        = mc_p[i]  + Mc;
		Mc_pp       = mc_pp[i] + Mc;
	}
	Mc			= Mc/double(aantal_midcellen);
	Mc_p		= Mc_p/double(aantal_midcellen);
	Mc_pp		= Mc_pp/double(aantal_midcellen);
	Delta		= Delta/double(aantal_midcellen);
	Delta_prev	= Delta_prev/double(aantal_midcellen);

	if(delta-delta_prev > 0.01)
	{
		om = relax_f_maxspeed;
	}
	else
    {
        om = fabs((om_p)*(Delta)*((Delta_prev-Delta)/((Delta_prev-Delta)*(Delta_prev-Delta))));
        if (om > relax_f_maxspeed) om = relax_f_maxspeed;
        if (om < relax_f_minslow)  om = relax_f_minslow;
     
    }
	//Update
	for(int i =0;i<aantal_midcellen;i++)		x_new [i]= (mc[i]*om) + (mc_p[i]*(1-om));
	
	return &x_new[0];
}


#endif