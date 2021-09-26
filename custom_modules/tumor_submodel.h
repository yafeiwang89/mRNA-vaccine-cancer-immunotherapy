#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

#include "./submodel_data_structures.h"

#include "./immune_submodels.h"

#ifndef __tumor_submodel__
#define __tumor_submodel__

extern Submodel_Information tumor_submodel_info;

void tumor_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt );
void tumor_phenotype( Cell* pCell, Phenotype& phenotype, double dt );
void tumor_mechanics( Cell* pCell, Phenotype& phenotype, double dt );

// this damage response will need to be added to the "infected cell response" model
void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt );

void tumor_submodel_setup( void );

#endif
