#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h"
#include "./immune_submodels.h"
#include "./receptor_dynamics.h"
#include "./internal_mRNA_dynamics.h"


#ifndef __DC_LNP_mRNA__
#define __DC_LNP_mRNA__

extern Submodel_Information DC_LNP_mRNA_info; 
	
void DC_LNP_mRNA_model_setup( void );

void introduce_DC_cells( void );

void inject_mRNA_loaded_LNP( void );

// don't put into individual cell models
void DC_LNP_mRNA_model( Cell* pCell, Phenotype& phenotype, double dt );

void DC_LNP_mRNA_model_phenotype_func( Cell* pCell, Phenotype& phenotype, double dt ); 

// this needs to be done on faster time scale; 
void DC_LNP_mRNA_main_model( double dt ); 


#endif 