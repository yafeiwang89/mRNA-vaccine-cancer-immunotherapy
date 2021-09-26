#include "./tumor_submodel.h"

using namespace PhysiCell;

std::string tumor_submodel_version = "0.1.0";

Submodel_Information tumor_submodel_info;

void tumor_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
	// elastic adhesions
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );

	return;
}

void tumor_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{

	static int start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	static int end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	static int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );  
	static int necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 

	static int debris_index = microenvironment.find_density_index( "debris");
	static int ISF_index = microenvironment.find_density_index( "ISF" );

	phenotype.secretion.secretion_rates[ISF_index] = parameters.doubles( "ISF_secretion_rate" ); // 0.1

	phenotype.motility.is_motile = false;

	// T-cell based death
	TCell_induced_apoptosis(pCell, phenotype, dt );

	// if I am dead, remove all adhesions
	if( phenotype.death.dead == true )
	{
		// detach all attached cells
		// remove_all_adhesions( pCell );
		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
		phenotype.secretion.secretion_rates[ISF_index] = 0.0; 

		pCell->functions.update_phenotype = NULL;
	}

	// oxygen based proliferation and death model
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	// std::cout << phenotype.death.rates[necrosis_index] << std::endl;

	return;
}

void tumor_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I'm dead, don't bother
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions,
		// since those are part of mechanics.
		// remove_all_adhesions( pCell );

		// Let's just fully disable now.
		pCell->functions.custom_cell_rule = NULL;
		pCell->functions.contact_function = NULL;	
		return;
	}

	return;
}

void tumor_submodel_setup( void )
{
	Cell_Definition* pCD;

	// set up any submodels you need

	// receptor trafficking
	//receptor_dynamics_model_setup(); // done
	// pathogen replication
	//internal_pathogen_model_setup();
	// single-cell response
	// internal_pathogen_response_model_setup();

	tumor_submodel_info.name = "tumor model";
	tumor_submodel_info.version = tumor_submodel_version;
	// set functions
	tumor_submodel_info.main_function = NULL;
	tumor_submodel_info.phenotype_function = tumor_phenotype;
	tumor_submodel_info.mechanics_function = tumor_mechanics;

	// what microenvironment variables do you expect?
	// tumor_submodel_info.microenvironment_variables.push_back( "TNF" );

	// what custom data do I need?
	//tumor_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	tumor_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "tumor cell" );
	pCD->functions.update_phenotype = tumor_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = tumor_submodel_info.mechanics_function;
	pCD->functions.contact_function = tumor_contact_function;

	// necrosis rate for tumor cells
	// static int necrosis_index = pCD->phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
	// pCD->parameters.max_necrosis_rate = 10; 

	return;
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	static int debris_index = microenvironment.find_density_index( "debris" );
	static int nAb = microenvironment.find_density_index( "Ig" ); 

	static bool CD8_T_enabled = parameters.bools( "cytotoxic_T_response" );  // false;
	static bool antibody_enabled = parameters.bools( "antibody_response" );  // false;

	static double Ig_50 = 0.6;
	static double apoptosis_rate = 5.32e-5; 
	static double Ig_induce_apoptosis = 0.001; 

	double Ig_conc = pCell->nearest_density_vector()[nAb];  
	double E_Ig = Ig_conc /( Ig_conc + Ig_50); 


	if ( antibody_enabled )
	{	
		// Hill function: (1- E)*r_0 + E*r_max
		phenotype.death.rates[apoptosis_index] = (1- E_Ig)*apoptosis_rate + E_Ig* Ig_induce_apoptosis;

	}


	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] && CD8_T_enabled )
	{
		// make sure to get rid of all adhesions!
		// detach all attached cells
		// remove_all_adhesions( pCell );

		#pragma omp critical
		{
			std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from T cell contact" << std::endl;
		}

		// induce death
		pCell->start_death( apoptosis_index );
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];

		pCell->functions.update_phenotype = NULL;
	}

	return;
}
