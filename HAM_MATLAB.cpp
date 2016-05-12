///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//     Gateway function for HAM.cpp
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This code uses the HAM function developped originally in C++ 
// It is based on the code of M. Steeman and Wufi

// Input of geometry and material properties is given in in the screen "GUI_HAMMAT.fig"

# include <iostream>            // std::endl
# include <fstream>     
# include <iomanip>             // std::setw   std::setiosflags
# include <new>                 // ::operator new[]
# include <math.h>
# include <string>
//# include "Geometry.h"
# include "Ham.h"
# include "Discretize.h"
# include "Matrix.h"
//# include "Store.h"
# include "mex.h" /* Always include this */

int HAM(double* Geometry,                           //INPUT 1
        double* ThicknessLayer,                     //INPUT 2          
        double  GridExponent,                       //INPUT 3
        double* NumberMeshCells,                    //INPUT 4
        double* DryDensity,                         //INPUT 5
        double* ThermalChar,                        //INPUT 6
        double* HygricChar,                         //INPUT 7
        double* InitialConditions,                  //INPUT 8
        double* BoundaryConditionsLeft,             //INPUT 9
        double* BoundaryConditionsRight,            //INPUT 10
        double* ParametersTime,                     //INPUT 11
        double* TimeOfMeasurement,                  //INPUT 12
        double *MeasuredWeight,                     //INPUT 13
        double *TimeOfCalculation,                  //OUTPUT
        double *CalulatedWeight,                    //OUTPUT
        double  LayerDiffussionCoeffToEstimate,     //INPUT
        double *CalculatedDiffusionCoefficient,     //OUTPUT
        int     NumberOfMeasurements                //INPUT
        );

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//     MATLAB FUNCTION
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void mexFunction
        (
            int nlhs,               /* Number of output variables */
            mxArray *plhs[],        /* Array of mxArray pointers to the output variables */ 
            int nrhs,               /* Number of input variables  */
            const mxArray *prhs[]   /* Array of mxArray pointers to the input variables  */
        ) 
{

         
    /* MACROS   for output and input arguments */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //OUTPUTS
    #define CALCULATED_WEIGHT_GAIN              plhs[0]
    #define CALCULATED_DIFFUSION_COEFF          plhs[1]
        
    //INPUTS
    #define GEOMETRY                            prhs[0]
    #define THICKNESSLAYER                      prhs[1]
    #define GRIDEXPONENT                        prhs[2]
    #define MESH                                prhs[3]
    #define DRYDENSITY                          prhs[4]
    #define THERMALCHAR                         prhs[5]
    #define HYGRICCHAR                          prhs[6]
    #define INITIALCONDITIONS                   prhs[7]
    #define BOUNDARYCONDITIONSLEFT              prhs[8]
    #define BOUNDARYCONDITIONSRIGHT             prhs[9]
    #define PARAMETERSTIME                      prhs[10]
    #define ESTIMATEDIFFUSIONCOEFF              prhs[11]
    #define MEASURED_WEIGHT_GAIN                prhs[12]
         
        
    
    //CHECK NUMBER INPUTS AND OUTPUTS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if(nrhs != 13)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:nrhs","Eleven input arguments required.");
    if(nlhs != 2)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:nlhs","Two output arguments required.");
    
    //INPUTS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //GEOMETRY
    if(!mxIsDouble(GEOMETRY) || mxIsComplex(GEOMETRY))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix GEOMETRY must be type double.");
    if(mxGetN(GEOMETRY) != 3)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input GEOMETRY must contain three columns.");
    double* Geometry                = new double[3];            // 1 = Hoogte   2 = Breedte   3= AantalLayers
    Geometry                        = mxGetPr(GEOMETRY);
    int NumberMaterialLayers        = int(Geometry[2]); 
    if(NumberMaterialLayers<1 || fmod(NumberMaterialLayers,1)!= 0.)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport","NumberMaterialLayers is not Correct (<1 or not Int): %f", NumberMaterialLayers);
    else mexPrintf("\n \n Geometry \n");
    for (int i =0; i<3;i++) mexPrintf("%8.4f ",Geometry[i]); 
    
    //THICKNESSLAYER
    if(!mxIsDouble(THICKNESSLAYER) || mxIsComplex(THICKNESSLAYER))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix THICKNESSLAYER must be type double.");
    if(mxGetM(THICKNESSLAYER) != NumberMaterialLayers)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input THICKNESSLAYER must contain three columns.");
    double* ThicknessLayer      = new double(NumberMaterialLayers);
    ThicknessLayer              = mxGetPr(THICKNESSLAYER); 
    mexPrintf("\n \n ThicknessLayer \n");
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",ThicknessLayer[i]); 
        
    //GRIDEXPONENT
    if(!mxIsDouble(GRIDEXPONENT) || mxIsComplex(GRIDEXPONENT) || mxGetNumberOfElements(GRIDEXPONENT) != 1)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","Input GRIDEXPONENT must be a type double.");
    if (!(mxIsScalar(GRIDEXPONENT)))
        mexErrMsgIdAndTxt( "MyToolbox:VapourTransport:invalidInputType","Input GRIDEXPONENT must be a scalar.");
    double GridExponent    = mxGetScalar(GRIDEXPONENT);
    mexPrintf("\n \n GridExponent \n");
    mexPrintf("%8.4f ",GridExponent);
   
    
    //MESH (=Number of Cells for Each Material Layer)
    if(!mxIsDouble(MESH) || mxIsComplex(MESH))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix MESH must be type double.");
    if(mxGetM(MESH) != NumberMaterialLayers)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input MESH must contain three columns.");
    double* NumberMeshCells     = new double[NumberMaterialLayers];
    NumberMeshCells             = mxGetPr(MESH); 
    int TotalNumberMeshCells    = 0;
    for(int i=0;i<NumberMaterialLayers;i++)     TotalNumberMeshCells = TotalNumberMeshCells + int(NumberMeshCells[i]);    
    mexPrintf("\n \n NumberMeshCells \n");
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",NumberMeshCells[i]); 
    
    
    //DRYDENSITY
    if(!mxIsDouble(DRYDENSITY) || mxIsComplex(DRYDENSITY))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix DRYDENSITY must be type double.");
    if(mxGetM(DRYDENSITY) != NumberMaterialLayers)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input DRYDENSITY must contain # rows.");
    double *DryDensity      = new double [NumberMaterialLayers];
    DryDensity              = mxGetPr(DRYDENSITY);
    mexPrintf("\n \n DryDensity \n");
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",DryDensity[i]); 
    
    
    //THERMALCHAR
    if(!mxIsDouble(THERMALCHAR) || mxIsComplex(THERMALCHAR))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix DRYDENSITY must be type double.");
    if(mxGetM(THERMALCHAR) != NumberMaterialLayers)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input DRYDENSITY must contain # columns.");
    double *ThermalChar     = new double [NumberMaterialLayers*3];
    ThermalChar             = mxGetPr(THERMALCHAR);
    mexPrintf("\n \n ThermalChar \n");
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",ThermalChar[i]); 
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",ThermalChar[NumberMaterialLayers+i]); 
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",ThermalChar[2*NumberMaterialLayers+i]/100.); 
    
    
    //HYGRICCHAR
    if(!mxIsDouble(HYGRICCHAR) || mxIsComplex(HYGRICCHAR))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix HYGRICCHAR must be type double.");
    if(mxGetM(HYGRICCHAR) != NumberMaterialLayers)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input HYGRICCHAR must contain # columns.");
    double *HygricChar      = new double [NumberMaterialLayers*4];
    HygricChar              = mxGetPr(HYGRICCHAR);
    mexPrintf("\n \n HygricChar \n");
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",HygricChar[i]); 
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",HygricChar[NumberMaterialLayers+i]); 
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",HygricChar[2*NumberMaterialLayers+i]); 
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",HygricChar[3*NumberMaterialLayers+i]); 
   
    
    //INITIALCONDITIONS (inside the sample)
    if(!mxIsDouble(INITIALCONDITIONS) || mxIsComplex(INITIALCONDITIONS))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix INITIALCONDITIONS must be type double.");
    if(mxGetM(INITIALCONDITIONS) != NumberMaterialLayers)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input INITIALCONDITIONS must contain # rows.");
    double *InitialConditions       = new double [NumberMaterialLayers*2];
    InitialConditions               = mxGetPr(INITIALCONDITIONS);
    mexPrintf("\n \n InitialConditions \n");
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",InitialConditions[i]); 
    for (int i =0; i<NumberMaterialLayers;i++) mexPrintf("%8.4f ",InitialConditions[NumberMaterialLayers+i]); 
    
    //BOUNDARYCONDITIONS
    //LEFT
    if(!mxIsDouble(BOUNDARYCONDITIONSLEFT) || mxIsComplex(BOUNDARYCONDITIONSLEFT))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix BOUNDARYCONDITIONSLEFT must be type double.");
    double *BoundaryConditionsLeft      = new double [5];
    BoundaryConditionsLeft              = mxGetPr(BOUNDARYCONDITIONSLEFT);
    mexPrintf("\n \n BoundaryConditionsLeft \n");
    for (int i =0; i<5;i++) mexPrintf("%8.4f ",BoundaryConditionsLeft[i]); 
    
    //RIGHT
    if(!mxIsDouble(BOUNDARYCONDITIONSRIGHT) || mxIsComplex(BOUNDARYCONDITIONSRIGHT))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix BOUNDARYCONDITIONSRIGHT must be type double.");
    double *BoundaryConditionsRight     = new double [5];
    BoundaryConditionsRight             = mxGetPr(BOUNDARYCONDITIONSRIGHT);
    mexPrintf("\n \n BoundaryConditionsRight \n");
    for (int i =0; i<5;i++) mexPrintf("%8.4f ",BoundaryConditionsRight[i]); 
   
    
    //PARAMETERSTIME
    if(!mxIsDouble(PARAMETERSTIME) || mxIsComplex(PARAMETERSTIME))
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","input matrix PARAMETERSTIME must be type double.");
    if(mxGetN(PARAMETERSTIME) != 3)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notColumnVector","Input PARAMETERSTIME must contain 3 columns.");
    double *ParametersTime      = new double [3];
    ParametersTime              = mxGetPr(PARAMETERSTIME);
    mexPrintf("\n \n ParametersTime \n");
    for (int i =0; i<5;i++) mexPrintf("%8.4f ",ParametersTime[i]); 
   
    
    //ESTIMATEDIFFUSIONCOEFF (Layer for which to estimate diff coeff)
    //If this is not checked on in the GUI, this value is "-1"
    if(!mxIsDouble(ESTIMATEDIFFUSIONCOEFF) || mxIsComplex(ESTIMATEDIFFUSIONCOEFF) || mxGetNumberOfElements(GRIDEXPONENT) != 1)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notDouble","Input ESTIMATEDIFFUSIONCOEFF must be a type double.");
    if (!(mxIsScalar(ESTIMATEDIFFUSIONCOEFF)))
        mexErrMsgIdAndTxt( "MyToolbox:VapourTransport:invalidInputType","Input ESTIMATEDIFFUSIONCOEFF must be a scalar.");
    double LayerDiffussionCoeffToEstimate    = mxGetScalar(ESTIMATEDIFFUSIONCOEFF);
    mexPrintf("\n \n LayerDiffussionCoeffToEstimate \n");
    mexPrintf("%8.4f ",LayerDiffussionCoeffToEstimate); 
   
    
    //MEASURED_WEIGHT_GAIN
    if(mxGetM(MEASURED_WEIGHT_GAIN) != 2)
        mexErrMsgIdAndTxt("MyToolbox:VapourTransport:notRowVector","Input MEASURED_WEIGHT_GAIN must contain two rows.");
    int NumberOfMeasurements = int(mxGetN(MEASURED_WEIGHT_GAIN));               // Get number of columns = number of measurements
    mexPrintf("\nTotal number of measurements: %i \n", NumberOfMeasurements);
    double *TimeDataMeas    = new double [2*NumberOfMeasurements];              // Store in 1D-array
    TimeDataMeas            = mxGetPr(MEASURED_WEIGHT_GAIN);                    // Pointer to the beginning of the matrix
    double *TimeMeasured    = new double [NumberOfMeasurements];                // 1xN input matrix (N = number of measurements)
    double *MeasuredWeightG = new double [NumberOfMeasurements];                // 1xN input matrix (N = number of measurements)
    for(int i=0;i<NumberOfMeasurements;i++)
    {
        mexPrintf("\n %8.4f %8.4f",TimeDataMeas[2*i], TimeDataMeas[(2*i)+1]);
        TimeMeasured[i]      = TimeDataMeas[2*i];
        MeasuredWeightG[i]   = TimeDataMeas[(2*i)+1];
    }
    mexPrintf("\n \n TimeMeasured \n");
    for (int i =0; i<int(NumberOfMeasurements);i++) mexPrintf("%8.4f ",TimeMeasured[i]);
    mexPrintf("\n MeasuredWeight \n");
    for (int i =0; i<int(NumberOfMeasurements);i++) mexPrintf("%8.4f ",MeasuredWeightG[i]);
    
    
    
    
    
    // OUTPUTS (= Matrix of calculated weight gain) 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Make space for the output variables
    
    //CALCULATED_WEIGHT_GAIN
    double *CalculatedWeightGain;       //1xN output matrix (N = number of measurements)
   	CALCULATED_WEIGHT_GAIN = mxCreateDoubleMatrix(0,0,mxREAL);                             //empty array
    mxSetM(CALCULATED_WEIGHT_GAIN,2);                                                      //total number of rows = 2
    mxSetN(CALCULATED_WEIGHT_GAIN,NumberOfMeasurements);                                   //total number of columns
    mxSetData(CALCULATED_WEIGHT_GAIN, mxMalloc(sizeof(double)*2*NumberOfMeasurements));    //Allocate memory
    CalculatedWeightGain = mxGetPr(CALCULATED_WEIGHT_GAIN);
    
    //CALCULATED_DIFFUSION_COEFF (If needed --> Input Screen HAMMAT)
    double *CalculatedDiffCoef;         //1xL (L = number of material layers)
    CALCULATED_DIFFUSION_COEFF = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetM(CALCULATED_DIFFUSION_COEFF,1);
    mxSetN(CALCULATED_DIFFUSION_COEFF,1);
    mxSetData(CALCULATED_DIFFUSION_COEFF, mxMalloc(sizeof(double)*1*1));
    CalculatedDiffCoef = mxGetPr(CALCULATED_DIFFUSION_COEFF);
    
    
    //CALL THE COMPUTATIONAL ROUTINE (= HAM-code in C++)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    double *TimeCalculated      = new double [NumberOfMeasurements];
    double *CalculatedWeight    = new double [NumberOfMeasurements];
        
    int ResultWeight =  HAM(Geometry, 
                            ThicknessLayer, 
                            GridExponent,
                            NumberMeshCells,
                            DryDensity,
                            ThermalChar,
                            HygricChar,
                            InitialConditions,
                            BoundaryConditionsLeft,
                            BoundaryConditionsRight,
                            ParametersTime,
                            TimeMeasured, 
                            MeasuredWeightG, 
                            TimeCalculated, 
                            CalculatedWeightGain,
                            LayerDiffussionCoeffToEstimate,
                            CalculatedDiffCoef,
                            NumberOfMeasurements);
    if (ResultWeight != 0) mexErrMsgIdAndTxt("MyToolbox:VapourTransport","HAM calculation went wrong. Check Logfiles");
    else mexPrintf("\n Calculation Succeeded. \n");
            
    return;

}


//  END GATEWAY FUNCTION
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    HAM FUNCTION
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int HAM(double* Geometry,                           //INPUT 1
        double* ThicknessLayer,                     //INPUT 2          
        double  GridExponent,                       //INPUT 3
        double* NumberMeshCells,                    //INPUT 4
        double* DryDensity,                         //INPUT 5
        double* ThermalChar,                        //INPUT 6
        double* HygricChar,                         //INPUT 7
        double* InitialConditions,                  //INPUT 8
        double* BoundaryConditionsLeft,             //INPUT 9
        double* BoundaryConditionsRight,            //INPUT 10
        double* ParametersTime,                     //INPUT 11
        double* TimeOfMeasurement,                  //INPUT 12
        double *MeasuredWeight,                     //INPUT 13
        double *TimeOfCalculation,                  //OUTPUT
        double *CalulatedWeight,                    //OUTPUT
        double  LayerDiffussionCoeffToEstimate,     //INPUT
        double *CalculatedDiffusionCoefficient,     //OUTPUT
        int     NumberOfMeasurements                //INPUT
        )
{
    
    //////////////////////////////////////////////////////////////////////////////////////
    // I. DEFINE VARIABLES
    //////////////////////////////////////////////////////////////////////////////////////
    
    
    // *********************************************************************************************************************************
    // GEOMETRY
    // *********************************************************************************************************************************
        //GEOMETRY     
        //Dimensions Sample [m]
        double Height           = Geometry[0];
        double Width            = Geometry[1];

        //Number of Material Layers [-]
        const int NumberLayers	= int(Geometry[2]);

        //THICKNESSLAYER [m]
        //double *ThicknessLayer = new double[NumberLayers];
        double TotalThicknessSample    = 0.;
        for(int i=0;i<NumberLayers;i++) TotalThicknessSample = TotalThicknessSample + ThicknessLayer[i];
    
    
    // *********************************************************************************************************************************
    // MESH
    // *********************************************************************************************************************************

        //Number of Cells for each layer
        int* NumberCellsLayer = new int[NumberLayers];
        for(int i=0;i<NumberLayers;i++)    NumberCellsLayer[i]      = int(NumberMeshCells[i]);
        
        //Total number of mesh cells and nodes (2 boundary nodes + midpoint mesh cells)
        int Cells    = 0;       
        for(int i=0;i<NumberLayers;i++) Cells = Cells + int(NumberCellsLayer[i]);
        const int Cellen = Cells;     
        const int Nodes  = Cells + 2;
    
        // Grid concentration (0 = equidistant mesh)
        double Gridconcentration  = GridExponent;

    
    // *********************************************************************************************************************************
    // MATERIAL CONDITIONS
    // *********************************************************************************************************************************

        //INITIAL CONDITIONS
        //------------------------------------------------------------------------------------------------------------------------------
        //Temperature[°C] & Relative humidity[%] in the material
        double *Temp_Start_Mat  = new double [NumberLayers];
        double *RH_Start_Mat    = new double [NumberLayers];
        for(int i=0;i<NumberLayers;i++)
        {
            Temp_Start_Mat[i]   = InitialConditions[i];
            RH_Start_Mat[i]     = InitialConditions[NumberLayers+i];        
        }

        //MATERIAL PROPERTIES
        //------------------------------------------------------------------------------------------------------------------------------
        //BASIC
        //Dry density [kg/m³]
        double *Ro_dry = new double [NumberLayers];
        for(int i=0;i<NumberLayers;i++)     Ro_dry[i] = DryDensity[i];

        //THERMAL
        //Specific heat capacity, Dry [J/kgK]
        double *Cp_mat = new double [NumberLayers];
        for(int i=0;i<NumberLayers;i++)     Cp_mat[i] = ThermalChar[i];

        //Thermal conductivity [W/mK]
        //Lambda = lambda_dry + lambda_wet_correction*water_content
        //Dry
        double *Lambda_dry = new double[NumberLayers];
        for(int i=0;i<NumberLayers;i++)     Lambda_dry[i] = ThermalChar[NumberLayers+i];       //Oak
    
        //Correction for wet material [%/%M-%]
        double *Lambda_wet_correction = new double[NumberLayers];
        for(int i=0;i<NumberLayers;i++)     Lambda_wet_correction[i] = (ThermalChar[2*NumberLayers+i]/100.);
    

        //HYGRO-THERMAL
        // SORPTION ISOTHERM: Moisture content [kg/m³] & moisture capacity
        // Het verband tussen de relatieve vochtigheid (') en het vochtgehalte (w) wordt weergegeven
        // door de sorptie-isotherm
        // --> Calculated by function: double *MoistureContent_kgm3 = new double [NumberLayers];
        //     w(RH) = A.RH²+B.RH+C
        double *Amat = new double[NumberLayers];
        double *Bmat = new double[NumberLayers];
        double *Cmat = new double[NumberLayers];
        for(int i=0;i<NumberLayers;i++)
        {
            Amat[i] = HygricChar[i];
            Bmat[i] = HygricChar[NumberLayers+i];
            Cmat[i] = HygricChar[2*NumberLayers+i];     
        }
               
        //VAPOUR DIFFUSION:Vapour Diffusion Resistance [kg/(m.s.Pa)]
        double *DiffusionCoefficient_Pressure = new double [NumberLayers];
        for(int i=0;i<NumberLayers;i++)         DiffusionCoefficient_Pressure[i] = HygricChar[3*NumberLayers+i];                            //Layer 1 = Chipboard (Averaged Diff Coeff)
        
        int     Layer_DiffusionCoefficient = int(LayerDiffussionCoeffToEstimate);                                 //Layer for which to calculate diffusion coeffcient  
        if(Layer_DiffusionCoefficient==-1)
        {
            //No Calculation
        }


    // *********************************************************************************************************************************
    // CLIMATE
    // *********************************************************************************************************************************
        //Temperature [°C]
        double Tleft        = BoundaryConditionsLeft[0];
        double Tright       = BoundaryConditionsRight[0];        
    
        //Relative humidity [%]
        double RHleft       = BoundaryConditionsLeft[1];
        double RHright      = BoundaryConditionsRight[1];

        //Atmospheric pressure [Pa]
        double Pleft        = BoundaryConditionsLeft[2];								// luchtdruk 101325Pa
        double Pright       = BoundaryConditionsRight[2];								// luchtdruk van de omgeving [Pa]
    
        //Transfer coefficients
        //Heat transfer [W/m²K]
        //Wufi beta     = 8
        //Wufi Outdoor  = 17
        double Alpha_Left	= BoundaryConditionsLeft[3];
        double Alpha_Right	= BoundaryConditionsRight[3];

        //Water vapour transfer [kg/m²sPa]
        //Fixed Value:
            // Beta_indoor  = 18.5e-9 s/m
            // Beta_outdoor = 140 e-9 s/m 
        //calculated by formula:
            // Default waarde voor bi =  1 (calculated by formula)
            // Default waarde voor be = -1
        //Wufi [Kunzel (2011)]: 7e-9*alpha_conv  (stel beta = 0)
        double Beta_Left	= BoundaryConditionsLeft[4];
        double Beta_Right	= BoundaryConditionsRight[4];
        
        //Windspeed [m/s]
        double WindspeedLeft    = 0.;
        double WindspeedRight   = 0.;
        
   
    
    // *********************************************************************************************************************************
    // ITERATION PROCEDURE
    // *********************************************************************************************************************************
        //Time functions [s]
        double StartTime            = ParametersTime[0];
        double TimeStep             = ParametersTime[1];              //1   min
        double StopTime             = ParametersTime[2];             //2594701 = 100 dagen
    
        //parameters voor de berekening tot convergentie
        int    Iteratie             = 0;
        int    MaxNumberIteratie    = 200;
        double DTmax                = 0.0001;												//convergentie-eis temperatuur
        double DPmax                = 0.0001;												//convergentie-eis dampdruk
        double DWeightmax           = 0.00001;                                              //convergentie-eis gewicht vocht materiaal
        double RelaxFactorT         = 1.;
        double RelaxFactorPv        = 1.;   
    
   
     // *********************************************************************************************************************************
     // LOGFILES
     // *********************************************************************************************************************************
    
        std::ofstream F_Error;
        F_Error.open("Errors.log");
    
        std::ofstream F_Mesh;
        F_Mesh.open("Mesh.log");
    
        std::ofstream F_Vochtkar ;
        F_Vochtkar.open("MoistChar.log");
	
        std::ofstream F_Warmtekar ;
        F_Warmtekar.open("HeatChar.log");
    
        std::ofstream F_Weight ;
        F_Weight.open("WeightSample.log");
    
        std::ofstream Z_values;
        Z_values.open("Z_values.log");
        Z_values.close();
    
        std::ofstream R_values;
        R_values.open("R_values.txt");
        R_values.close();
    
        std::ofstream Te;
        Te.open("Temp.txt");
        Te.close();
    
        std::ofstream P;
        P.open("Pvap.txt");
        P.close();
    
        std::ofstream F_RH;
        F_RH.open("RelHum");
        F_RH.close();
    
        std::ofstream F_C;
        F_C.open("Moist.txt");
        F_C.close();
    
        std::ofstream a_values;
        a_values.open("a_values.txt");
        a_values.close();
		
        std::ofstream b_values;
        b_values.open("b_values.txt");
        b_values.close();
		
        std::ofstream c_values;
        c_values.open("c_values.txt");
        c_values.close();
		
        std::ofstream d_values;
        d_values.open("d_values.txt"); 
        d_values.close();
    
        std::ofstream aV_values;
        aV_values.open("aV_values.txt");
        aV_values.close();
		
        std::ofstream bV_values;
        bV_values.open("bV_values.txt");
        bV_values.close();
		
        std::ofstream cV_values;
        cV_values.open("cV_values.txt");
        cV_values.close();
            
        std::ofstream dV_values;
        dV_values.open("dV_values.txt");
        dV_values.close();
        
     // *********************************************************************************************************************************
     // TIME & ITERATION PROCESS
     // *********************************************************************************************************************************
    
        // TIME & ITERATION PROCESS
        double Time             = StartTime;
        double Time_Previous    = 0.;
        int t                   = 0; //teller van measured array
        double *Ctemp           = new (std::nothrow)  double[Cellen];
        double *Cpres           = new (std::nothrow)  double[Cellen];
        double *Cweight         = new (std::nothrow)  double[Cellen];
        double dtemp            = 0.;
        double dpress           = 0.;
        double dweight          = 0.;       //Difference between 
    
        //Ham-model
        int check1              = 0;
        int check2              = 0;
        int check3              = 0;
        int check               = 0;
    
        //Error in Calculated Weight
        double *Residuals       = new double [NumberOfMeasurements];
        double SSE              = 0.;
        double SSE_min          = 0.;  
        double CalculateDiffusionCoefficient = 0.;
        
     // *********************************************************************************************************************************
     // PROPERTIES OF EACH NODE
     // *********************************************************************************************************************************
    
        // Temperature Related
        double *T           = new (std::nothrow) double [Nodes];                    // Temperatuur voor iedere knoop
        double *TPI			= new (std::nothrow) double [Nodes]; 
        double *TPII        = new (std::nothrow) double [Nodes];
        double *TPT         = new (std::nothrow) double [Nodes];
        double *Tvapour     = new (std::nothrow) double[Cellen];
        double *Lambda      = new (std::nothrow) double [Nodes];                    // lambda-waarde per knoop
        double *LambdaPI    = new (std::nothrow) double [Nodes];                    // lambda-waarde per knoop
        double *R           = new (std::nothrow) double [Nodes-1];                  // warmteweerstand [tussenknoopwaarde]
        double  R_tot       = 0;
	    
        double *Ro          = new (std::nothrow) double[Nodes];
        double *Cp          = new (std::nothrow) double[Nodes];                     // warmtecapaciteit  - waarde per knoop 
        double *Energy      = new (std::nothrow) double[Nodes];                     // inwendige energie - waarde per knoop
        double *EnergyPI    = new (std::nothrow) double[Nodes];                     // inwendige energie - waarde per knoop vorige iteratiestap
        double *EnergyPT    = new (std::nothrow) double[Nodes];                     // inwendige energie - waarde per knoop vorige tijdstap
	
        double q_init       = 0.;
        double q            = 0.;                                                    // heat flux density [W/m²]
    
    
        //Matrix for Temperature
        double *apT	= new (std::nothrow) double[Cellen];
        double *aeT	= new (std::nothrow) double[Cellen];
        double *awT	= new (std::nothrow) double[Cellen];                            // aantal  = cellen
        double *dT	= new (std::nothrow) double[Cellen];                            // aantal  = cellen
    
        //Vapour Flow related
        double *rh              = new (std::nothrow) double [Cellen];               // hulpwaarde die eerst de rh in de knopen (dus niet randknopen) uitrekent
        double *RH              = new (std::nothrow) double [Nodes];				// RH [knoopwaarde]
        double *RHPI            = new (std::nothrow) double [Nodes];				// knoopwaarde vorige iteratiestap
        double *Psat            = new (std::nothrow) double [Nodes];				// Saturated Vapour Pressure	
        double *PsatPI          = new (std::nothrow) double [Nodes];				// knoopwaarde vorige iteratiestap
        double *Pv              = new (std::nothrow) double [Nodes];				// Vapour Pressure
        double *PvPI            = new (std::nothrow) double [Nodes];				// knoopwaarde vorige iteratiestap
        double *PvPII           = new (std::nothrow) double [Nodes];
        double *PvPT            = new (std::nothrow) double [Nodes];
        double *Pvapour         = new (std::nothrow) double [Cellen];
        double *Perm_air        = new (std::nothrow) double [Nodes];  
        double *Perm_mat        = new (std::nothrow) double [Nodes];
        double *VapRes          = new (std::nothrow) double [Nodes+1];              // tussenknoopwaarde
        double *Z               = new (std::nothrow) double [Nodes+1];              // Vapour resistance [tussenknoopwaarde]
        double  Z_tot           = 0;                                                // Total Vapour Resisitance of the material
    
        double *C               = new (std::nothrow) double [Nodes];				// Moisture capacity [kg/m³]
        double *CPI             = new (std::nothrow) double [Nodes];				// knoopwaarde vorige iteratiestap
        double *CPII            = new (std::nothrow) double [Nodes];
        double *CPT             = new (std::nothrow) double [Nodes];				// knoopwaarde vorige tijdstap
        double *RoKsi           = new (std::nothrow) double [Nodes];				// knoopwaarde
        double *RoKsiPI         = new (std::nothrow) double [Nodes];				// knoopwaarde vorige iteratiestap
    
        double *AmountMoisture  = new (std::nothrow) double [NumberLayers];
   
        double g1               = 0.;
        double g2               = 0.;
        double g1P              = 0.;
        double g2P              = 0.;
        double g                = 0.;                                               // Water Vapour Flux density [kg/m²s]
	
        //Matrix for Vapour Pressure
        double *apV	=new (std::nothrow) double[Cellen];
        double *aeV	=new (std::nothrow) double[Cellen];
        double *awV	=new (std::nothrow) double[Cellen];
        double *dV	=new (std::nothrow) double[Cellen];
    
        //Total Weight of Sample
        double DryWeight    = 0.;
        double WeightSample = 0.;
        
        //CHECK LOCAL VARIABLES
        if (T == nullptr)
        {
            F_Error << "Error 1: memory T could not be allocated"<<std::endl;
            return 1;
        }
        if (TPI == nullptr)
        {
            F_Error << "Error 1: memory TPI could not be allocated"<<std::endl;
        	return 1;
    	}
    	if (Lambda == nullptr)
        {
    		F_Error << "Error 1: memory lambda could not be allocated"<<std::endl;
    		return 1;
    	}
    	if (LambdaPI == nullptr)
        {
            F_Error << "Error: memory lambdaPI could not be allocated"<<std::endl;
            return 1;
        }
        if (R == nullptr)
        {
            F_Error << "Error: memory R could not be allocated"<<std::endl;
            return 1;
        }
        if (Ro == nullptr)
        {
        	F_Error << "Error 1: memory ro could not be allocated"<<std::endl;
            return 1;
        }
        if (Cp == nullptr)
        {
            F_Error << "Error: memory cp could not be allocated"<<std::endl;
            return 1;
        }
        if (Energy == nullptr)
        {
            F_Error << "Error: memory energy could not be allocated"<<std::endl;
            return 1;
        }
        if (EnergyPI == nullptr)
        {
            F_Error << "Error: memory EnergyPI could not be allocated"<<std::endl;
            return 1;
        }
        if (EnergyPT == nullptr)
        {
            F_Error << "Error: memory EnergyPT could not be allocated"<<std::endl;
            return 1;
        }
        F_Error<<"GEOMETRY"<<std::endl;
        F_Error<<"----------------------------------------------------"<<std::endl;
        F_Error<<Height<<"   "<<Width<<"   "<<std::endl<<NumberLayers<<std::endl;
        for(int i=0;i<NumberLayers;i++) F_Error<<ThicknessLayer[i]<<" ";
        F_Error<<std::endl;
        F_Error<<std::endl;
        
        F_Error<<"MESH"<<std::endl;
        F_Error<<"----------------------------------------------------"<<std::endl;
        F_Error<<std::endl;
        
        F_Error<<"INITIAL  CONDITIONS"<<std::endl;
        F_Error<<"----------------------------------------------------"<<std::endl;
        F_Error<<"TEMPERATURE [°C] & RELATIVE HUMIDITY [%]"<<std::endl;
        
        for (int i=0; i < 2*NumberLayers; i++)
        {
                F_Error<<"InitialConditions[i]: "<<InitialConditions[i]<<std::endl;
        }
        F_Error<<std::endl;
        for (int l=0; l < NumberLayers; l++)
        {
                F_Error<<"Temp_Start_Mat[l]: "<<Temp_Start_Mat[l]<<std::endl;
                F_Error<<"RH_Start_Mat[l]: "<<RH_Start_Mat[l]<<std::endl;
                F_Error<<"NumberCellsLayer[l]: "<<NumberCellsLayer[l]<<std::endl;
        }
        F_Error<<std::endl;
        
        F_Error<<"HYGRO-THERMAL"<<std::endl;
        F_Error<<"----------------------------------------------------"<<std::endl;
        for (int l=0; l < NumberLayers; l++)
        {
            F_Error<<Amat[l]<<"   "<<Bmat[l]<<"   "<<Cmat[l]<<"   "<<DiffusionCoefficient_Pressure[l]<<std::endl;
        }
        F_Error<<Layer_DiffusionCoefficient<<std::endl;
        F_Error<<std::endl;  
        
        F_Error<<"CLIMATE"<<std::endl;
        F_Error<<"----------------------------------------------------"<<std::endl;
        F_Error<<"Tleft: "<<Tleft<<"   ";
        F_Error<<"Tright: "<<Tright<<std::endl;
        F_Error<<"RHleft: "<<RHleft<<"    ";
        F_Error<<"RHright: "<<RHright<<std::endl;
        F_Error<<"Pleft: "<<Pleft<<"   ";
        F_Error<<"Pright: "<<Pright<<std::endl;
        F_Error<<"Alpha_Left: "<<Alpha_Left<<"   ";
        F_Error<<"Alpha_Right: "<<Alpha_Right<<std::endl;
        F_Error<<"Beta_Left: "<<Beta_Left<<"   ";
        F_Error<<"Beta_Right: "<<Beta_Right<<std::endl;
        
        F_Error<<"ITERATION PROCEDURE"<<std::endl;
        F_Error<<"----------------------------------------------------"<<std::endl;
        F_Error<<"StartTime: "<<StartTime<<"   ";
        F_Error<<"StopTime: "<<StopTime<<std::endl;
   
    //////////////////////////////////////////////////////////////////////////////////////
    // II. DISCRETISATION OF WALL
    //////////////////////////////////////////////////////////////////////////////////////
        
        // distance between nodes
        double *DeltaX  = new (std::nothrow) double [Nodes-1];    
        if (DeltaX == nullptr)
        {
            F_Error << "Error: memory DeltaX could not be allocated"<<std::endl;
            return 1;
        }
        // length of one meshcell
        double *Length  = new (std::nothrow) double [Nodes];
        if (Length == nullptr)
        {
            F_Error << "Error: memory Length could not be allocated"<<std::endl;
            return 1;
        }
        //Discretise
        int Mesh_Ok = discretize(NumberLayers,Cells,NumberCellsLayer,ThicknessLayer,Gridconcentration,Length,DeltaX);
        if (Mesh_Ok != 0)
        {
            F_Error <<"Mesh not ok!"<<std::endl;
            return 1;
        }
        else
        {
            F_Mesh<<"Length of all the meshcells"<<std::endl;
            for (int i=0; i < Nodes; i++)       
            F_Mesh <<Length[i]<<"  "<<std::endl;
        }
        F_Mesh.close();
    
      
    //////////////////////////////////////////////////////////////////////////////////////
    //III. SET INITIAL VALUES
    //////////////////////////////////////////////////////////////////////////////////////
   
    // *********************************************************************************************************************************
    // STATE VARIABLES NODES
    // *********************************************************************************************************************************
    
        //TEMPERATURE [°C] & RELATIVE HUMIDITY [%]
        T[0]        = Tleft;
        T[Nodes-1]  = Tright;
        
        RH[0]       = RHleft;
        RH[Nodes-1] = RHright;
        int x=0;
        
        for (int l=0; l < NumberLayers; l++)
        {
            for (int i=0; i < NumberCellsLayer[l]; i++)
            {
                T[i+x+1]	= Temp_Start_Mat[l];
                RH[i+x+1]	= RH_Start_Mat[l];
            }
            x = x + NumberCellsLayer[l];
        }
        for (int i=1; i<(Nodes-1); i++)
        {        
            if (RH[i] < 0)          F_Error<<"Wrong input RH. Smaller dan 0!"<<std::endl;
            else if (RH[i] < 1)     F_Error<<"RH is smaller than 1. Remember base 100"<<std::endl;
            else if (RH[i] > 100)   F_Error<<"RH is larger than 100%"<<std::endl;
        }
    
    
        // VAPOUR PRESSURE [Pa]
        Pv[0]          = pv(T[0],RH[0]);
        Pv[Nodes-1]    = pv(T[Nodes-1],RH[Nodes-1]);
        for (int i=1; i<(Nodes-1);i++)		Pv[i]     = pv(T[i],RH[i]);
	
        //WRITE OUT RESULTS	
        Te.open("Temp.txt",std::ios::app);
        Te<<"START"<<std::endl;
    	Te<<Time<<"     "<<TimeOfMeasurement[t]<<"     ";
        for(int i= 0;i<Nodes;i++)		Te<<std::setprecision(5)<<T[i]<<"      ";
        Te<<std::endl;
        Te.close();
	
        P.open("Pvap.txt",std::ios::app);
        P<<"START"<<std::endl;
        P<<Time<<"     ";
        for(int i= 0;i<Nodes;i++)		P<<std::setprecision(5)<<Pv[i]<<"      ";
        P<<std::endl;
        P.close();
        
        F_RH.open("RelHum.txt",std::ios::app);
        F_RH<<Time<<"     ";
        for(int i= 0;i<Nodes;i++)		F_RH<<std::setprecision(5)<<RH[i]<<"      ";
        F_RH<<std::endl;
        F_RH.close();
        
        
        //*
        // CALCULATE VAPOUR DIFUSSION RESISTANCE
        Perm_air[0]        = 0.;
        Perm_air[Nodes -1] = 0.;
        Perm_mat[0]		   = 0.;
        Perm_mat[Nodes -1] = 0.;
        x = 0;
        for(int l=0; l < NumberLayers; l++)
        {
            for(int i=0; i < NumberCellsLayer[l]; i++)			Perm_mat[i+x+1] = DiffusionCoefficient_Pressure[l];
            x = x + NumberCellsLayer[l];
        }
    
        // Moisture content C[kg/m³]
        // Moisture capacity RoKsi[kg/m³] (derivative sorption isotherm to rel. humidity)
        C[0]				= 0.;       
        C[Nodes-1]		    = 0.;
        RoKsi[0]			= 0.;     
        RoKsi[Nodes -1]     = 0.;
    
        x=0;
        for (int l=0;l < NumberLayers;l++)
        {
            AmountMoisture[l] = 0.;
            for (int i=0; i < NumberCellsLayer[l]; i++)
            {
                C[x+i+1]            = MoistureContent_kgm3 (Cmat[l],Bmat[l],Amat[l],RH[x+i+1]);
                RoKsi[x+i+1]        = MoistureCapacity_kgm3(Cmat[l],Bmat[l],Amat[l],RH[x+i+1]);
                AmountMoisture[l]   = AmountMoisture[l] + C[x+i+1]*Length[x+i+1];                       //(kg/m³)*m = kg/m²
            }
            x = x + NumberCellsLayer[l];
        }
        ///*
        
        F_C.open("Moist.txt",std::ios::app);    //initial Moisture Content
        F_C<<Time<<"     ";
        for(int i= 0;i<Nodes;i++)		F_C<<std::setprecision(5)<<C[i]<<"      ";
        F_C<<std::endl;
        F_C.close();
        
        F_Vochtkar <<"Time: " << Time << std::endl<<std::endl;
        F_Vochtkar <<"Node ";
        F_Vochtkar <<std::setw(15)<<"Perm [kg/(m?Pa)]";
        F_Vochtkar <<std::setw(15)<<"Roksi[kg/m³]";
        F_Vochtkar <<std::endl;
        for (int i=0; i< Nodes;i++)
        {
            //F_Vochtkar <<i+1;
            //F_Vochtkar <<std::setprecision(3)<<std::setw(18)<<Perm_mat[i];
            F_Vochtkar <<std::setprecision(3)<<std::setw(18)<<RoKsi[i];
            //F_Vochtkar <<std::endl;
        }
        F_Vochtkar <<std::endl;
        //F_Vochtkar.close();
        
        
    // *********************************************************************************************************************************
    // THERMAL CHAR
    // *********************************************************************************************************************************
        //CALCULATE HEAT RESISTANCE
        Lambda[0]        = 0.;
        Lambda[Nodes -1] = 0.;
        Ro[0]			 = 0.;
        Ro[Nodes -1]	 = 0.;
        Cp[0]            = 0.;
        Cp[Nodes -1]     = 0.;
        Energy[0]        = 0.;
        Energy[Nodes -1] = 0.;
	
        x=0;
        for(int l=0; l<NumberLayers;l++)
        {
            for(int i=0; i < NumberCellsLayer[l]; i++)
            {
                Lambda[x+i+1]   = Lambda_dry[l] + C[x+i+1]*Lambda_wet_correction[l];
                Ro[x+i+1]       = Ro_dry[l];										//[kg/m³]
                Cp[x+i+1]       = Cp_mat[l]+((C[x+i+1]*Cpliq)/Ro[x+i+1]);
                Energy[x+i+1]   = ((Ro[x+i+1]*Cp_mat[l])+(C[x+i+1]*Cpliq))*T[x+i+1];				//[J/m³]
            }
            x = x + NumberCellsLayer[l];
        }
    
        F_Warmtekar <<"Time: " << Time << std::endl<<std::endl;
        F_Warmtekar <<"Node ";
        F_Warmtekar <<std::setw(12)<<"Lambda []";
        F_Warmtekar <<std::setw(12)<<"Cp[]";
        F_Warmtekar <<std::setw(12)<<"Energy [J]" ;
        F_Warmtekar <<std::endl;
        for (int i=0; i< Nodes;i++)
        {
            F_Warmtekar <<i+1;
            F_Warmtekar <<std::setprecision(3)<<std::setw(12)<<Lambda[i];
            F_Warmtekar <<std::setprecision(3)<<std::setw(12)<<Cp[i];
            F_Warmtekar <<std::setiosflags(std::ios::scientific);
            F_Warmtekar <<std::setprecision(3)<<std::setw(12)<<Energy[i];
            F_Warmtekar <<std::endl;
        }
        F_Warmtekar.close();
    
        //CALCULATE (DRY) WEIGHT
        for(int l=0; l<NumberLayers;l++)		DryWeight = DryWeight + (Ro_dry[l]*Height*Width*ThicknessLayer[l]);										//[kg/m³]
        WeightSample = 0.;
        for (int l=0; l< NumberLayers; l++)		WeightSample = WeightSample + AmountMoisture[l];    //kg/m?
        WeightSample = WeightSample*Height*Width;
        t = 0;
        TimeOfCalculation[t]    = Time;
        CalulatedWeight[t]      = DryWeight + WeightSample;
        t++;
        F_Weight << Time;
        F_Weight << std::setw(12)<<DryWeight;
        F_Weight << std::setw(12)<<WeightSample;
        F_Weight << std::setw(12)<<DryWeight + WeightSample;
        F_Weight << std::endl;
    
   ///*   
    //////////////////////////////////////////////////////////////////////////////////////
    // I. CALCULATE MOISTURE AND HEAT TRANSPORT: LOOP FROM STARTTIME TO STOPTIME
    //////////////////////////////////////////////////////////////////////////////////////
        do
        {
            //////////////////////////////////////////////////////////////////////////////////////
            //A. SET VALUES FIRST ITERATION STEP
            //////////////////////////////////////////////////////////////////////////////////////
        
            // IF THIS IS THE STARTTIME (Time = 0)
            // ---------------------------------------------------------------------------------------------------------------------------------
       
            //SET FIRST TIMESTEP
            if(Time == StartTime) Time = Time+TimeStep;
       
       
       
            //IF THIS IS THE FIRST TIMESTEP (Time = Timestep):
            // ---------------------------------------------------------------------------------------------------------------------------------
            //		values are similar to initial values; 
            //		first timestep/ first iterationstep: there are no values stored yet in the storage array
            //		inital values will funtion as 'values of previous timestep'
        
            if(Time == StartTime +  TimeStep)
            {
                //Set values Previous iteration step and previous time
                for (int i=0; i<Nodes; i++)
                {
                    TPI[i]      = T[i];
                    TPT[i]      = T[i];
                    RHPI[i]     = RH[i];
                    PvPI[i]     = Pv[i];
                    PvPT[i]     = Pv[i];
                    PsatPI[i]   = Psat[i];
                    CPI[i]      = C[i];
                    CPT[i]      = C[i];
                    RoKsiPI[i]  = RoKsi[i];
                    EnergyPI[i] = Energy[i];
                    EnergyPT[i] = Energy[i];
                }
            
                //store initial values in the storage array
                //Store::setStorageVars(Time,+11*Nodes+2);
                double Time_Previous = 0;
            }//end first timestep
            
            //////////////////////////////////////////////////////////////////////////////////////
            //B. UPDATE VALUES NEXT TIMESTEPS
            //////////////////////////////////////////////////////////////////////////////////////
            else if(Time > StartTime + TimeStep)
            {
                for (int i = 0; i < Nodes; i++)
                {
                    TPT[i]          = TPI[i];
                    CPT[i]          = CPI[i];
                }
                //g1P					= Store::getStorageVars( +11*Nodes);
                //g2P					= Store::getStorageVars( +11*Nodes+1);
                //Time_Previous		= Store::getStorageVars( +11*Nodes+2);
                
                //TEMPERATURE
    
                //BOUNDARY VALUES NEW TIMESTEP (FROM CLIMATE FILE)
                T[0]            = Tleft;
                T[Nodes-1]      = Tright;
            
                RH[0]           = RHleft;
                RH[Nodes-1]     = RHright;
            
                Psat[0]        = pvsat(T[0]);
                Psat[Nodes-1]  = pvsat(T[Nodes-1]);
	
                Pv[0]          = pv(T[0],RH[0]);
                Pv[Nodes-1]    = pv(T[Nodes-1],RH[Nodes-1]);  
				
                
            }// einde setting variables in timestep/iterationstep 
    
            //////////////////////////////////////////////////////////////////////////////////////
            //C. CALCULATE HEAT AND VAPOUR FLOW
            //////////////////////////////////////////////////////////////////////////////////////
            Iteratie = 0;
            dtemp    = 10.;
            dpress   = 1000.;
            check    = 0;
            RelaxFactorT         = 1.;
            RelaxFactorPv        = 1.;
        
            //ONE TIMESTEP!!!!!!!!!!!! (ALL ITERATIONS IN HERE)
            do
            {
                Iteratie++;
                
                // *********************************************************************************************************************************
                // 1) HEAT AND VAPOUR RESISTANCE
                // *********************************************************************************************************************************
            
                //Z-values
                //Berekenen van Beta indien deze waarden niet aangepast zijn in de Input-file
                Beta_Left   = Beta_Int(Beta_Left,Alpha_Left,Tleft,T[1],WindspeedLeft);
                Beta_Right  = Beta_Int(Beta_Right,Alpha_Right,Tright,T[Nodes-2],WindspeedRight);
                Z[0]        = (1/Beta_Left)   + ((Length[1]/2.)/Perm_mat[1]);
                Z[Cellen]   = (1/Beta_Right)  + ((Length[Cellen-1]/2.)/Perm_mat[Cellen]);   
                for(int i=1; i<Cellen ; i++)       	Z[i] = ((Length[i]/2.)/Perm_mat[i])+((Length[i+1]/2.)/Perm_mat[i+1]);
            
            
                //R-values
                R[0]        = (1/Alpha_Left)  + ((Length[1]/2.)/Lambda[1]);
                R[Cellen]   = (1/Alpha_Right) + ((Length[Cellen-1]/2.)/Lambda[Cellen]);
                for(int i=1; i<Cellen; i++)			R[i]  = ((Length[i]/2.)/Lambda[i])+((Length[i+1]/2.)/Lambda[i+1]);
			
        
                // *********************************************************************************************************************************
                // 2) Temperature : X x T = D
                // *********************************************************************************************************************************
            	
                //Coefficienten
                for(int i=1;i<Cellen+1;i++)
                {
                    awT[i-1] =  1/R[i-1];
                    apT[i-1] = -1/R[i-1]-1/R[i]-(Ro[i]*Cp[i]*Length[i]/TimeStep);
                    aeT[i-1] =  1/R[i];
                }
                for(int i=2;i<Cellen;i++)
                {
                    //dT[i-1]  = (-Ro[i]*Cp[i]*Length[i]*TPT[i]) -  Hev*(((Pv[i-1]-Pv[i])/Z[i])-((Pv[i]-Pv[i+1])/Z[i+1]));
                    dT[i-1]  = (-Ro[i]*Cp[i]*Length[i]*TPT[i]/TimeStep);
                    //dT[i-1]  = (-Ro[i]*Cp[i]*Length[i]*TPI[i]) -  Hev*(((Pv[i-1]-Pv[i])/Z[i])-((Pv[i]-Pv[i+1])/Z[i+1]))+ ((EnergyPI[i]-EnergyPT[i])*Length[i]);
                    //dT[i-1]  = (-Ro[i]*Cp[i]*Length[i]*TPI[i]/TimeStep) + ((EnergyPI[i]-EnergyPT[i])*Length[i]/TimeStep);  
                }
                //dT[0]		   = (-Ro[1]*Cp[1]*Length[1]*TPT[1]/TimeStep) -  Hev*(((Pv[0]-Pv[1])/Z[0])-((Pv[1]-Pv[2])/Z[1])) - awT[0]*T[0];
                //dT[Cellen-1] = (-Ro[Cellen]*Cp[Cellen]*Length[Cellen]*TPT[Cellen]/TimeStep) - Hev*(((Pv[Cellen-1]-Pv[Cellen])/Z[Cellen-1])-(Pv[Cellen+1]-Pv[Cellen])/Z[Cellen])- aeT[Cellen-1]*T[Cellen+1];
                dT[0]		   = (-Ro[1]*Cp[1]*Length[1]*TPT[1]/TimeStep)  - awT[0]*T[0];
                dT[Cellen-1] = (-Ro[Cellen]*Cp[Cellen]*Length[Cellen]*TPT[Cellen]/TimeStep) - aeT[Cellen-1]*T[Cellen+1];
                //dT[0]		   = (-Ro[1]*Cp[1]*Length[1]*TPI[1]) -  Hev*(((Pv[0]-Pv[1])*TimeStep/Z[0])-((Pv[1]-Pv[2])/Z[1]))+ (((EnergyPI[1]-EnergyPT[1])/TimeStep)*Length[1])-awT[0]*T[0];
                //dT[Cellen-1] = (-Ro[Cellen]*Cp[Cellen]*Length[Cellen]*TPI[Cellen]/TimeStep) - Hev*(((Pv[Cellen-1]-Pv[Cellen])/Z[Cellen-1])-((Pv[Cellen]-Pv[Cellen+1])/Z[Cellen])) + ((EnergyPI[Cellen]-EnergyPT[Cellen])*Length[Cellen]/TimeStep)-aeT[Cellen-1]*T[Cellen+1];
                //dT[0]		   = (-Ro[1]*Cp[1]*Length[1]*TPI[1]/TimeStep)+ ((EnergyPI[1]-EnergyPT[1])*Length[1])-awT[0]*T[0];
                //dT[Cellen-1]   = (-Ro[Cellen]*Cp[Cellen]*Length[Cellen]*TPI[Cellen]/TimeStep)+ ((EnergyPI[Cellen]-EnergyPT[Cellen])*Length[Cellen])-aeT[Cellen-1]*T[Cellen+1];
		
                //Gauss-Seidel om de waarden voor T te vinden
                //for(int i=0;i<Cellen;i++)			Tvapour[i]  = TPI[i+1];
                //Tvapour = GSOR(awT,apT,aeT,TPI,dT,RelaxFactor,Cellen);
                Tvapour = Thomas_algorithm(awT,apT,aeT,dT,Cellen);
                //T = relaxation_Aitken(T, TPI, TPII, RelaxFactorT, Nodes);
                for(int i=0;i<Cellen;i++)			T[i+1]      = Tvapour[i];
                RelaxFactorT = 0.7;
                for(int i=1;i<Cellen+1;i++) T[i] = (T[i]*RelaxFactorT) + (TPI[i]*(1-RelaxFactorT));
                T[0]        = Tleft;
                T[Nodes-1]  = Tright;
     
                
                // *********************************************************************************************************************************
                // 3) Calculation of Vapour Pressure : Y x Pv = D
                // *********************************************************************************************************************************
            	for(int i=1; i<(Cellen+1); i++)
                {
                    awV[i-1] =   1/Z[i-1];
                    apV[i-1] =  -1/Z[i-1]-1/Z[i]-(RoKsi[i]*Length[i]/(TimeStep*pvsat(TPI[i])));
                    aeV[i-1] =   1/Z[i];
                }
                for(int i=2; i<(Cellen);i++)
                {
                    //dV[i-1]   = -(RoKsi[i]*Length[i]*(PvPT[i]/pvsat(TPT[i]))); // in de eerste tijdstap, de eerste it is PI=PT
                    dV[i-1]   = -(RoKsi[i]*Length[i]*(PvPI[i]/pvsat(TPI[i]))/TimeStep) + ((CPI[i]-CPT[i])*Length[i])/TimeStep; // in de eerste tijdstap, de eerste it is PI=PT
                }
		
                //dV[0]		    = -(RoKsi[1]*Length[1]*(PvPT[1]/pvsat(TPT[1]))/TimeStep) - awV[0]*Pv[0];
                //dV[Cellen-1]  = -(RoKsi[Cellen]*Length[Cellen]*(PvPT[Cellen]/pvsat(TPT[Cellen]))/TimeStep) -aeV[Cellen-1]*Pv[Cellen+1];
                dV[0]		    = -(RoKsi[1]*Length[1]*(PvPI[1]/pvsat(TPI[1]))/TimeStep) + ((CPI[1]-CPT[1])*Length[1]/TimeStep)- awV[0]*Pv[0];
                dV[Cellen-1]    = -(RoKsi[Cellen]*Length[Cellen]*(PvPI[Cellen]/pvsat(TPI[Cellen]))/TimeStep) + ((CPI[Cellen]-CPT[Cellen])*Length[Cellen]/TimeStep)-aeV[Cellen-1]*Pv[Cellen+1];
		    
                
                std::ofstream bV_values;
                bV_values<<std::setprecision(10);
                bV_values.open("bV_values.txt",std::ios::app);
                for(int i=0;i<Cellen;i++)			bV_values<<apV[i]<<std::endl;
                bV_values.close();
            
                std::ofstream dV_values;
                dV_values.open("dV_values.txt",std::ios::app);
                dV_values<<std::setprecision(10);
                dV_values<<Time<<"   "<<Iteratie<<"  "<<dV[0]<<"  "<<-(RoKsi[1]*Length[1]*(PvPI[1]/pvsat(TPI[1]))/TimeStep)<<"  "<<((CPI[1]-CPT[1])*Length[1]/TimeStep)<<"               "<<PvPI[1]<<"   "<<Pv[1]<<"   "<<-awV[0]<<" "<<Pv[0]<<std::endl;
                for(int i=2; i<(Cellen);i++)
                {
                    dV_values<<Time<<"   "<<Iteratie<<"  "<<dV[i-1]<<"  "<<-(RoKsi[i]*Length[i]*(PvPI[i]/pvsat(TPI[i]))/TimeStep)<<"  "<<((CPI[i]-CPT[i])*Length[i]/TimeStep)<<"               "<<PvPI[i]<<"   "<<Pv[i]<<std::endl;
                }
                dV_values<<Time<<"   "<<Iteratie<<"  "<<dV[Cellen-1]<<"  "<<-(RoKsi[Cellen]*Length[Cellen]*(PvPI[Cellen]/pvsat(TPI[Cellen]))/TimeStep)<<"  "<<((CPI[Cellen]-CPT[Cellen])*Length[Cellen]/TimeStep)<<"               "<<PvPI[Cellen]<<"   "<<Pv[Cellen]<<"   "<<-aeV[Cellen-1]<<" "<<Pv[Cellen+1]<<std::endl;
                dV_values.close();
                
            
                //Gauss-Seidel om de waarden voor Pv te vinden
                //for(int i=0;i<Cellen;i++)	Pvapour[i]  = PvPI[i+1];
                //GSOR(awV,apV,aeV,Pvapour,dV,RelaxFactor,Cellen);
                Pvapour = Thomas_algorithm(awV,apV,aeV,dV,Cellen);
                for(int i=0 ; i<Cellen ; i++) 	Pv[i+1]     = Pvapour[i];
                RelaxFactorPv = 1;
                for(int i=1;i<Cellen+1;i++) Pv[i] = (Pv[i]*RelaxFactorPv) + (PvPI[i]*(1-RelaxFactorPv));
                       
                RH[0]           = RHleft;
                RH[Nodes-1]     = RHright;
                Pv[0]           = pv(T[0],RH[0]);
                Pv[Nodes-1]     = pv(T[Nodes-1],RH[Nodes-1]);
            
                for (int i=1; i<Cellen+1; i++)
                {
                    if(Pv[i] < 0)                   Pv[i] = 1.0e-5;
                    else if(Pv[i] > pvsat(T[i]))	Pv[i] = pvsat(T[i]);
                }
            
            
                // *********************************************************************************************************************************
                // 4)  UPDATE the MATERIAL PROPERTIES  for the next iteration step in the same timestep
                // *********************************************************************************************************************************
                // Moisture content C[kg/m³]
                // Moisture capacity RoKsi[kg/m³] (derivative sorption isotherm to rel. humidity)
                C[0]				= 0.;       
                C[Nodes-1]		    = 0.;
                RoKsi[0]			= 0.;     
                RoKsi[Nodes -1]     = 0.;
	
                //C[i] vochtgehalte in de knopen van de cellen
                x=0;
                for (int l=0;l < NumberLayers;l++)
                {
                    AmountMoisture[l] = 0.;
                    for (int i=0; i < NumberCellsLayer[l]; i++)
                    {
                        RH[x+i+1]           = 100*(Pv[x+i+1]/pvsat(T[i]));
                        C[x+i+1]            = MoistureContent_kgm3 (Cmat[l],Bmat[l],Amat[l],RH[x+i+1]);
                        RoKsi[x+i+1]        = MoistureCapacity_kgm3(Cmat[l],Bmat[l],Amat[l],RH[x+i+1]);
                        AmountMoisture[l]   = AmountMoisture[l] + C[x+i+1]*Length[x+i+1];  //(kg/m?*m = kg/m?
                    }
                    x = x + NumberCellsLayer[l];
                }
			
            
                for (int i=0; i< Nodes;i++) F_Vochtkar <<std::setprecision(3)<<std::setw(18)<<RoKsi[i];
                F_Vochtkar <<std::endl;
                        
            
                //Heat resistance
                //Boundary values ->0
                Lambda[0]        = 0.;
                Lambda[Nodes -1] = 0.;
                Ro[0]			 = 0.;
                Ro[Nodes -1]	 = 0.;
                Cp[0]            = 0.;
                Cp[Nodes -1]     = 0.;
                Energy[0]        = 0.;
                Energy[Nodes -1] = 0.;
            
                x=0;
                for(int l=0; l<NumberLayers;l++)
                {
                    for(int i=0; i < NumberCellsLayer[l]; i++)
                    {
                        Lambda[x+i+1]   = Lambda_dry[l] + C[x+i+1]*Lambda_wet_correction[l];
                        Ro[x+i+1]       = Ro_dry[l];                                                        //[kg/m³]
                        Cp[x+i+1]       = Cp_mat[l]+((C[x+i+1]*Cpliq)/Ro[x+i+1]);
                        Energy[x+i+1]   = ((Ro[x+i+1]*Cp_mat[l])+(C[x+i+1]*Cpliq))*T[x+i+1];				//[J/m³]
                    }
                    x = x + NumberCellsLayer[l];
                }

                //Difussion resistance
                Perm_air[0]        = 0.;
                Perm_air[Nodes -1] = 0.;
                Perm_mat[0]		   = 0.;
                Perm_mat[Nodes -1] = 0.;
                x = 0;
                for(int l=0; l < NumberLayers;l++)
                {
                    for(int i=0; i < NumberCellsLayer[l]; i++)  Perm_mat[i+x+1] = DiffusionCoefficient_Pressure[l];
                    x = x + NumberCellsLayer[l];
                }
            
            
                // *********************************************************************************************************************************
                // 5)  CONTROLEREN OF OPLOSSING GECONVERGEERD IS
                // *********************************************************************************************************************************
                
                for(int i=1; i<Cellen+1; i++)
                {
                    Ctemp[i-1] = abs(T[i]  - TPI[i]);
                    Cpres[i-1] = abs(Pv[i] - PvPI[i]);
                }
                //opsporen grootste waarde voor ctemp en cpres:
                dtemp   = Ctemp[0];
                dpress  = Cpres[0];
                for(int i=1; i<Cellen; i++)
                {
                    if (Ctemp[i]>dtemp)			dtemp   = Ctemp[i];
                    if (Cpres[i]>dpress)       	dpress  = Cpres[i];
                } 

                // *********************************************************************************************************************************
                // 6)  UPDATE STATE VARIABLES (VALUES FOUND IN PREVIOUS ITERATION STEP)
                // *********************************************************************************************************************************
                // Nieuw berekende waarden opslaan in the "value previous  iteration step"
            
                for (int i=0; i<Nodes; i++)
                {
                    TPII[i]		= TPI[i];
                    TPI[i]		= T[i];
                    RHPI[i]		= RH[i];
                    PvPII[i]	= PvPI[i];
                    PvPI[i]		= Pv[i];
                    CPI[i]		= C[i];
                    RoKsiPI[i]	= RoKsi[i];
                    EnergyPI[i]	= Energy[i];
                }
            
                if (dtemp > DTmax) check1 = 0;
                else check1 = 1;
                if (dpress > DPmax) check2 = 0;
                else check2 = 1;
                if (Iteratie > MaxNumberIteratie) check3 = 1;
                else check3 = 0;
                check = (check1*check2)+check3;
            
            
            
                //uitrekenen Psat en RH op basis van de nieuw bekomen temperatuur [enkel knopen in materiaal, niet de 'boundary values']
                
                std::ofstream F_Test;
                F_Test.open("test.log",std::ios::app);
                F_Test<<Time<<"    "<<Iteratie<<"  "<<check1<<"  "<<check2<<"  "<<check3<<"  "<<check;
                for (int i=1; i<Cellen+1; i++)
                {
                    F_Test<<std::setprecision(3)<<std::setw(10)<<Pv[i];
                    //F_Test<<std::setprecision(3)<<std::setw(10)<<PvPI[i];
                }
            
                for (int i=1; i<Cellen+1; i++)
                {
                    F_Test<<std::setprecision(3)<<std::setw(10)<<RH[i];
                }
                F_Test<<std::endl;
                F_Test.close();
                
            
            
            }while (check == 0);    //einde while-lus: indien resultaten voldoende convergent naar volgende stap
        
            //////////////////////////////////////////////////////////////////////////////////////
            //D. TIMESTEP DONE: RESULTS WRITE OUT 
            //////////////////////////////////////////////////////////////////////////////////////
            F_Error<<Time<<"   "<<Iteratie<<"   "<<dtemp<<"   "<<dpress<<std::endl;
    
            x=0;
            for (int l=0; l<NumberLayers; l++)
            {
                AmountMoisture[l] = 0.;
                for(int i=0; i<NumberCellsLayer[l]; i++)       AmountMoisture[l]   = AmountMoisture[l] + C[x+i+1]*Length[x+i+1];  //(kg/m?*m = kg/m?
            }
            
            WeightSample = 0.;
            for (int l=0; l< NumberLayers; l++)		WeightSample = WeightSample + AmountMoisture[l];    //kg/m?
            WeightSample = WeightSample*Height*Width;
        
            if(Time == TimeOfMeasurement[t])
            {
                TimeOfCalculation[t]    = Time;
                CalulatedWeight[t]      = WeightSample;
                t++;
            }
        
            //Write out Weight
            std::ofstream F_Weight ;
            F_Weight.open("WeightSample.log",std::ios::app);
            F_Weight << Time;
            F_Weight << std::setw(12)<<DryWeight;
            F_Weight << std::setw(12)<<WeightSample;
            F_Weight << std::setw(12)<<DryWeight + WeightSample;
            F_Weight << std::endl;
            F_Weight.close();
        
            //Write out Temp and Vap Pressure
            std::ofstream Te;
            Te.open("Temp.txt",std::ios::app);
            Te<<Time<<"     ";
            for(int i= 0;i<Nodes;i++)		Te<<std::setprecision(5)<<T[i]<<"      ";
            Te<<std::endl;
            Te.close();
	
            std::ofstream P;
            P.open("Pvap.txt",std::ios::app);
            P<<Time<<"     ";
            for(int i= 0;i<Nodes;i++)		P<<std::setprecision(5)<<Pv[i]<<"      ";
            P<<std::endl;
            P.close();
        
            std::ofstream F_RH;
            F_RH.open("RelHum.txt",std::ios::app);
            F_RH<<Time<<"     ";
            for(int i= 0;i<Nodes;i++)		F_RH<<std::setprecision(5)<<RH[i]<<"      ";
            F_RH<<std::endl;
            F_RH.close();
        
            //std::ofstream C_RH;
            F_C.open("Moist.txt",std::ios::app);
            F_C<<Time<<"     ";
            for(int i= 0;i<Nodes;i++)		F_C<<std::setprecision(5)<<C[i]<<"      ";
            F_C<<std::endl;
            F_C.close();
        
            /*
            std::ofstream Test_Tvalues;
            Test_Tvalues.open("Test_Tvalues.txt",std::ios::app);
            Test_Tvalues<<Time<<"  "<<Iteratie<<"   ";
            for(int i=0;i<Cellen+2;i++)			Test_Tvalues<<T[i]<<std::endl;
            Test_Tvalues.close();
            
            std::ofstream R_values;
            R_values.open("R_values.txt",std::ios::app);
            for(int i=0;i<Cellen+1;i++)			R_values<<R[i]<<std::endl;
            R_values.close();
        
            std::ofstream a_values;
            a_values.open("a_values.txt",std::ios::app);
            a_values<<std::setprecision(10);
            for(int i=0;i<Cellen;i++)			a_values<<awT[i]<<std::endl;
            a_values.close();
		
            std::ofstream b_values;
            b_values.open("b_values.txt",std::ios::app);
            b_values<<std::setprecision(10);
            for(int i=0;i<Cellen;i++)			b_values<<apT[i]<<std::endl;
            b_values.close();
		
            std::ofstream c_values;
            c_values<<std::setprecision(10);
            c_values.open("c_values.txt",std::ios::app);
            for(int i=0;i<Cellen;i++)			c_values<<aeT[i]<<std::endl;
            c_values.close();
		
            std::ofstream d_values;
            d_values<<std::setprecision(10);
            d_values.open("d_values.txt",std::ios::app); 
            for(int i=0;i<Cellen;i++)			d_values<<dT[i]<<std::endl;
            d_values<<std::endl;
            d_values.close();
            */
        
            // txt van de matrix
            /*
            std::ofstream Z_values;
            Z_values.open("Z_values.log",std::ios::app);
            Z_values<<Iteratie<<"   "<<Beta_Left<<"  "<<Beta_Right<<"  ";
            for(int i=0;i<Cellen+1;i++)			Z_values<<Z[i]<<"  ";
            Z_values<<std::endl;
            Z_values.close();
            
            std::ofstream aV_values;
            aV_values.open("aV_values.txt",std::ios::app);
            aV_values<<std::setprecision(10);
            for(int i=0;i<Cellen;i++)			aV_values<<awV[i]<<std::endl;
            aV_values.close();
		
            std::ofstream bV_values;
            bV_values<<std::setprecision(10);
            bV_values.open("bV_values.txt",std::ios::app);
            for(int i=0;i<Cellen;i++)			bV_values<<apV[i]<<std::endl;
            bV_values.close();
		
            std::ofstream cV_values;
            cV_values.open("cV_values.txt",std::ios::app);
            cV_values<<Iteratie<<"   ";
            cV_values<<std::setprecision(10);
            for(int i=0;i<Cellen;i++)			cV_values<<aeV[i];
            cV_values<<std::endl;
            cV_values.close();
            
            std::ofstream dV_values;
            dV_values.open("dV_values.txt",std::ios::app);
            dV_values<<std::setprecision(10);
            for(int i=0;i<Cellen;i++)			dV_values<<dV[i]<<"  "<<-(RoKsi[i]*Length[i]*(PvPI[i]/PsatPI[i]))<<"  "<<((CPI[i]-CPT[i])*Length[i])<<std::endl;
            dV_values.close();
            */
        
       
        
            Time = Time + TimeStep;

        }while(Time<StopTime);
    
    	
            
    //////////////////////////////////////////////////////////////////////////////////////
    // II. COMPARE THE CALCULATED WEIGHT WITH THE MEASURED WEIGHT
    //////////////////////////////////////////////////////////////////////////////////////
    /*
        if(LayerDiffussionCoeffToEstimate >=0)
        {
            //Residuals
            for(int i=0;i<NumberOfMeasurements;i++)       Residuals[i] = (CalulatedWeight[i] - MeasuredWeight[i]);
    
            //Sum of Squares
            SSE = 0;        
            for(int i=0;i<NumberOfMeasurements;i++) SSE = SSE + Residuals[i]*Residuals[i];
            if(SSE < SSE_min)
            {
                SSE_min = SSE;  
                //CalculateDiffusionCoefficient = DiffusionCoefficient_Pressure[l];
            }
    
            //Propose new diffusion coefficient
        }
    */
    
    
    /*    
        F_Error<<"tot eind. Joepi  ";
        
        
        //CLEAR MEMORY
        //CLEAR MEMORY
        delete [] apT   ;
        delete [] aeT   ;
        delete [] awT   ;
        delete [] dT    ;
        delete [] apV   ;
        delete [] aeV   ;
        delete [] awV   ;
        delete [] dV;
        delete [] AmountMoisture;

        delete [] C;
        delete [] CPI;
        delete [] CPT;
        delete [] Cp_mat;
        delete [] Ctemp;
        delete [] Cpres;
        delete [] Cweight;

        delete [] DiffusionCoefficient_Pressure;

        delete [] Lambda_dry;
        delete [] Lambda_wet_correction;

        delete [] Lambda;
        delete [] LambdaPI;

        delete [] Cp;
        
        delete [] Energy;
        delete [] EnergyPI;
        delete [] EnergyPT;
               
        delete [] NumberCellsLayer;

        delete [] Psat;
        delete [] PsatPI;
        delete [] Pv;
        delete [] PvPI;
        delete [] PvPII;
        delete [] PvPT;
        delete [] Pvapour;
        delete [] Perm_air;
        delete [] Perm_mat;

        delete [] R;
        delete [] rh;
        delete [] RH;
        delete [] RHPI;
        delete [] RH_Start_Mat;
        delete [] Ro_dry;
        delete [] Ro;

        delete [] RoKsi;
        delete [] RoKsiPI;


        delete [] ThicknessLayer;

        delete [] T;
        delete [] TPI;
        delete [] TPII;
        delete [] TPT;
        delete [] Tvapour;

        delete [] VapRes;
        delete [] Z;
        
        
        delete [] DeltaX;
        delete [] Length;


        F_Error.close();
        F_Vochtkar.close();
        F_Warmtekar.close();
        F_Weight.close();
    */
           
    
    
    
    //Als alles OK is
    return 0;

}

