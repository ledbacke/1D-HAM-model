//include guard
#pragma once
#ifndef HAM_H
#define HAM_H


//CONSTANTS
//****************************************************************************
double const Pi	    = 4*atan(1.0);	// Pi
const double Hev	= 2501000;      // latent heat of evaporation			[J/kg]
double Rv           = 461.5;        // specifieke gasconstante waterdamp (J/kgK)               
double Cpliq        = 4187;         // Heat capacity of the liquid phase [J/kgK]
double Cpvap        = 1926;         // Heat capacity of the vapour phase [J/kgK]   

//FUNCTIONS
//****************************************************************************
double max(double a, double b)
{
    double x = a > b ? a:b;
    return (x);
}
double pvsat(double T)
//T    = Temperature in °C
{ 
	double k;
	double pvsat;
	k=4.39553-6.2442*(double(T+273.16)/1000.)+9.953*pow((double(T+273.16)/1000.),2)-5.151*pow((double(T+273.16)/1000.),3);
	pvsat = 217.5e5*exp(2.3026*k*(1-(647.4/(T+273.16))));
	//pvsat = exp(65.8904-(7066.27/(T+273.15))-(5.976*log((T+273.15))));
    return (pvsat);
}

double pv(double T, double RH)
//T    = Temperature in °C
//Pvap = Vapour pressure in Pa
//RH   = % (base 100)
{
	return (pvsat(T)*(RH/100));
}

// Heat Capacity
double cp(double cp_mat, double w, double ro_dry)
{
	return (cp_mat+((w*Cpliq)/ro_dry));			//[J/kgK]
}

//Moisture capacity [kg/m³] - (derivative sorption isotherm to rel. humidity)
double MoistureCapacity_kgm3(double A,double B, double C, double RH)
{
    double rh = RH/100.;
    double Cap = ((1*(A+(B*rh)+(C*rh*rh)))-(rh*(B+2*C*rh)))/((A+(B*rh)+(C*rh*rh))*(A+(B*rh)+(C*rh*rh)));
    return Cap;
}

//Energy   
double energy(double ro, double cp, double T)
{
	return (ro*cp*T);
}

//Moisture Content [kg/m³]
double MoistureContent_kgm3(double A,double B, double C,double RH)
{
    double rh = RH/100.;
    double MC = rh/(A+(B*rh)+(C*rh*rh));
    return MC;
};

//Vapour Transfer Coeff
double Beta_Int(double beta, double alpha_int, double T_int, double Tsurf,double windspeed)
{
	double b = 0.;
    if(beta ==  1 || beta == 1)          b = (27+0.73*(T_int-Tsurf))*0.000000001;        //in
    else if(beta == 0)
    {        
        double alpha_int_conv = alpha_int -4.5;
        double alpha_min = 3.5;
        b = 7.e-9 * max(alpha_min,alpha_int_conv);
    }
    else b = beta;
	return b;
}
double Beta_Ext(double beta, double alpha_ext, double T_ext, double Tsurf,double windspeed)
{
	double b = 0.;
    if(beta == -1 || beta == 1)     b = 49.9*0.000000001*pow(windspeed,0.85);   //out
    else if(beta == 0)
    {
        double alpha_ext_conv = alpha_ext -6.5;
        double alpha_min = 3.5;
        b = 7.e-9 * max(alpha_min,alpha_ext_conv);
    }
    else b = beta;
	return b;
}

double Gemiddelde_mat(double mat1, double mat2, double lenght1, double lenght2)
{
           return (lenght1+lenght2)*(mat1*mat2)/(mat1*lenght2+mat2*lenght1);
}



#endif