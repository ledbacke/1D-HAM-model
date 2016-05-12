#pragma once
#ifndef DISCRETIZE_H
#define DISCRETIZE_H

#include <math.h>       // ceil 
#include <new>          // std::nothrow

int discretize( int      layers,
                int      cellen,
                int    * numbercellslayer,
                double * thicknesslayer,
                double   gridconcentration,
                double * length,
                double * deltax
        );


int discretize( int      layers,
                int      cellen,
                int    * numbercellslayer,
                double * thicknesslayer,
                double   gridconcentration,
                double * length,
                double * deltax
        )
{
    //Total number of mesh cells have to be Odd
    double even_cells = 0 ;
	for (int i=0; i<layers; i++)
    {
		even_cells = numbercellslayer[i]/2.;
		numbercellslayer[i]= int(2*ceil(even_cells));
	}
    
    //CALCULATE GRID
    double noemer           = 0.;
    int CellsPerLayer       = 0;            //Number of cells in layer i
    int TotalNumberCells    = 0;            //Total number of cells in the construction element
    double* LengthCell      = new(std::nothrow) double[cellen];           //Length of the cell
    
    ///*
    //Divide layers into meshcells
    for (int l=0;l<layers;l++)
    {
		noemer = 0;
		CellsPerLayer = int(numbercellslayer[l]);
				
		for (int i=0; i<(CellsPerLayer/2);i++)
		{
			noemer = noemer + (2*exp(double(-gridconcentration)*(1.-(2.*double((1+i))/double(numbercellslayer[l])))));
		}
	    
		for (int i=0; i<(CellsPerLayer/2); i++)
            LengthCell[TotalNumberCells+i] = (exp(double(-gridconcentration)*(1.-(2.*double((1+i))/double(numbercellslayer[l])))))*thicknesslayer[l]/noemer;
        for (int i=0; i<(CellsPerLayer/2); i++)				  
            LengthCell[TotalNumberCells+(CellsPerLayer-1)-i]=LengthCell[TotalNumberCells+i];
		TotalNumberCells = TotalNumberCells+CellsPerLayer;
	}
    
    ///*
    //Distance between nodes
    // node = center of a mesh cell
    deltax[0]                   = LengthCell[0]/2.;
    deltax[TotalNumberCells]    = LengthCell[TotalNumberCells-1]/2.;  
    for (int i=1; i<TotalNumberCells;i++)  deltax[i] = ((LengthCell[i-1])/2.)+((LengthCell[i])/2.);
    
    //Distance for a cell coupled to a node
    const int nodes = TotalNumberCells +2;
    length[0] = 0;
	length[nodes-1] = 0;
    for (int i=1; i<(nodes-1); i++) 		length[i]=LengthCell[(i-1)];
	
    delete [] LengthCell; 
	//*/
    
    return 0;        
};

#endif