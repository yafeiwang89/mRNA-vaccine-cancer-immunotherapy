#include "./receptor_dynamics.h" 

using namespace PhysiCell; 

std::string receptor_model_version = "0.2.0"; 

Submodel_Information receptor_dynamics_info; 

void receptor_dynamics_model_setup( void )
{
		// set version 
	receptor_dynamics_info.name = "receptor dynamics"; 
	receptor_dynamics_info.version = receptor_model_version; 
		// set functions 
	receptor_dynamics_info.main_function = receptor_dynamics_main_model; 
	receptor_dynamics_info.phenotype_function = NULL; // pushed into the "main" model  
	receptor_dynamics_info.mechanics_function = NULL; 	
		// what microenvironment variables do you need 

		// what cell variables and parameters do you need? 

	receptor_dynamics_info.cell_variables.push_back( "uptaken_LNP" ); 

	receptor_dynamics_info.cell_variables.push_back( "unbound_external_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound_external_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "unbound_internal_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound_internal_ACE2" ); 
	
	receptor_dynamics_info.cell_variables.push_back( "ACE2_binding_rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2_endocytosis_rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2_cargo_release_rate" ); 	
	receptor_dynamics_info.cell_variables.push_back( "ACE2_recycling_rate" ); 
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	receptor_dynamics_info.register_model(); 
	
	return; 
}

void receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int DC_type = get_cell_definition( "DC" ).type; 
	
	// bookkeeping -- find microenvironment variables we need
	static int LNP_external = microenvironment.find_density_index( "Lipid-nanoparticles" ); 
	static int LNP_endosome = pCell->custom_data.find_variable_index( "uptaken_LNP" ); 

	// bookkeeping -- find custom data we need 
	
	static int nR_EU = pCell->custom_data.find_variable_index( "unbound_external_ACE2" ); 
	static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_ACE2" ); 
	static int nR_IU = pCell->custom_data.find_variable_index( "unbound_internal_ACE2" ); 
	static int nR_IB = pCell->custom_data.find_variable_index( "bound_internal_ACE2" ); 
	
	static int nR_bind = pCell->custom_data.find_variable_index( "ACE2_binding_rate" ); 
	static int nR_endo = pCell->custom_data.find_variable_index( "ACE2_endocytosis_rate" ); 
	static int nR_release = pCell->custom_data.find_variable_index( "ACE2_cargo_release_rate" ); 	
	static int nR_recycle = pCell->custom_data.find_variable_index( "ACE2_recycling_rate" ); 
	
	// do nothing if dead 
	if( phenotype.death.dead == true )
	{ return; } 

	// if not Dendritic cell, do nothing         
	if( pCell->type != DC_type )
	{ return; } 
	
	// actual model goes here 
	
	// internalized LNP tells us how many have recently bound to receptors 
	double newly_bound = phenotype.molecular.internalized_total_substrates[LNP_external]; 
	// if it tried to bind to more LNP than there are receptors, compensate 
	double excess_binding = newly_bound - pCell->custom_data[nR_EU]; 
	if( excess_binding > 0.0 )
	{
		// don't bring in more LNP than there are receptors 
		newly_bound = pCell->custom_data[nR_EU]; 
		// dump any excess back into the microenvironment
		static double one_LNP_to_density = 1.0 / microenvironment.mesh.dV; 
		// this needs omp critical because 2 cells writing to 1 voxel is not thread safe 
		#pragma omp critical
		{
			pCell->nearest_density_vector()[LNP_external] += excess_binding * one_LNP_to_density; 
		}
	}
	phenotype.molecular.internalized_total_substrates[LNP_external] = 0.0; 
	
	// add newly bound receptor to R_EB	
	pCell->custom_data[nR_EB] += newly_bound; 
	
	// remove newly bound receptor from R_EU 
	pCell->custom_data[nR_EU] -= newly_bound; 

	
	// endocytosis 
	
	double dR_IB = dt*pCell->custom_data[nR_endo]*pCell->custom_data[nR_EB];
	if( dR_IB > pCell->custom_data[nR_EB] )
	{ dR_IB = pCell->custom_data[nR_EB]; }
	pCell->custom_data[nR_EB] -= dR_IB; // move from external bound
	pCell->custom_data[nR_IB] += dR_IB; // move to internal bound
	
	// viral release from endosomes 
	
	double dR_IU = dt*pCell->custom_data[nR_release]*pCell->custom_data[nR_IB];
	if( dR_IU > pCell->custom_data[nR_IB] )
	{ dR_IU = pCell->custom_data[nR_IB]; }
	pCell->custom_data[nR_IB] -= dR_IU; // move from internal bound 
	pCell->custom_data[nR_IU] += dR_IU; // move to internal unbound 
	pCell->custom_data[LNP_endosome] += dR_IU; // release LNP in the endosome 
	
	// receptor recycling 
	
	double dR_EU = dt*pCell->custom_data[nR_recycle]*pCell->custom_data[nR_IU];
	if( dR_EU > pCell->custom_data[nR_IU] )
	{ dR_EU = pCell->custom_data[nR_IU]; }
	pCell->custom_data[nR_IU] -= dR_EU; // move from internal unbound 
	pCell->custom_data[nR_EU] += dR_EU; // move to external unbound 
	
	// update the LNP uptake rate, this needs to be modified based on nanobio project !!!!!!!!!!!!!!!!
	
	phenotype.secretion.uptake_rates[LNP_external] = 
		pCell->custom_data[nR_bind] * pCell->custom_data[nR_EU]; 
	
	return; 
}

void receptor_dynamics_main_model( double dt )
{
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{ receptor_dynamics_model( pC, pC->phenotype , dt ); }
	}
	
	return; 
}
	
