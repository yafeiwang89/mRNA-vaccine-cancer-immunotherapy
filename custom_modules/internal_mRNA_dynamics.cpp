#include "./internal_mRNA_dynamics.h" 

using namespace PhysiCell; 

std::string internal_mRNA_dynamics_version = "0.4.0"; 

Submodel_Information internal_mRNA_dynamics_info; 


void internal_mRNA_model_setup( void )
{
		// set version
	internal_mRNA_dynamics_info.name = "internal mRNA dynamics"; 
	internal_mRNA_dynamics_info.version = internal_mRNA_dynamics_version; 
		// set functions 
	internal_mRNA_dynamics_info.main_function = NULL; 
	internal_mRNA_dynamics_info.phenotype_function = internal_mRNA_model; 
	internal_mRNA_dynamics_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_mRNA_dynamics_info.register_model();
	
	return; 
}

void internal_mRNA_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int epitope_external = microenvironment.find_density_index( "epitope" ); 

	static int LNP_endosome = pCell->custom_data.find_variable_index( "uptaken_LNP" ); 
	static int LNP_cytoplasm = pCell->custom_data.find_variable_index( "free_LNP" ); 
	static int LNP_uncoated = pCell->custom_data.find_variable_index( "uncoated_LNP" ); 
	static int mRNA_index = pCell->custom_data.find_variable_index( "mRNA" ); 
	static int epitope_index = pCell->custom_data.find_variable_index( "antigenic_epitope" ); 



/*	
	static bool done = false; 
	extern Cell* pInfected; 
	if( pCell == pInfected && 1 == 0 )
	{
		std::cout << std::endl << "viral dynamics : " << __LINE__ << " " 
			<< phenotype.molecular.internalized_total_substrates[ nV_external ] << " " 
			<< phenotype.molecular.internalized_total_substrates[ nA_external ] << " " 
			<< pCell->custom_data[nV_internal] << " " 
			<< pCell->custom_data[nUV] << " " 
			<< pCell->custom_data[nR] << " " 
			<< pCell->custom_data[nP] << " " 	
			<< pCell->custom_data[nA_internal] << " " 
			<< std::endl; 		
	}
*/	
	
	// copy virions from "internalized variables" to "custom variables"
/*	
	pCell->custom_data[nV_internal] = 
		phenotype.molecular.internalized_total_substrates[nV_external]; 
	// this transfer is now handled in receptor dynamics 
*/		

  // This resets the internal assembled virion count 
  // so we are commenting it out 
	// pCell->custom_data[nA_internal] = 
	// 	phenotype.molecular.internalized_total_substrates[nA_external]; 
		
	// actual model goes here 

	// LNP escape from endosome into cytoplasm
  double dEndo = dt * parameters.doubles("LNP_escape_rate") * pCell->custom_data[LNP_endosome];
	if( dEndo > pCell->custom_data[LNP_endosome] )
	{ dEndo = pCell->custom_data[LNP_endosome]; } 
	pCell->custom_data[LNP_endosome] -= dEndo; 
	pCell->custom_data[LNP_cytoplasm] += dEndo; 


	// uncoat LNP within the cytoplasm
	double dCyto = dt * parameters.doubles("uncoated_LNP_rate") * pCell->custom_data[LNP_cytoplasm] ;
	if( dCyto > pCell->custom_data[LNP_cytoplasm] )
	{ dCyto = pCell->custom_data[LNP_cytoplasm]; } 
	pCell->custom_data[LNP_cytoplasm] -= dCyto; 
	pCell->custom_data[LNP_uncoated] += dCyto; 

	// convert uncoated LNP to usable mRNA (assume one LNP has one mRNA)
	pCell->custom_data[mRNA_index] = pCell->custom_data[LNP_uncoated]; 

	
	// use mRNA to create epitope protein 
	double dP = dt * parameters.doubles("epitope_translation_rate") * pCell->custom_data[mRNA_index];
	pCell->custom_data[epitope_index] += dP; 


	// degrade mRNA and protein 
	
	
/*
	// set export rate 
	phenotype.secretion.net_export_rates[epitope_external] = 
		parameters.doubles("epitope_export_rate") * pCell->custom_data[epitope_index]; 
*/

	
	return; 
}






