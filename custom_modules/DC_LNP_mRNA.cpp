#include "./DC_LNP_mRNA.h" 
#include <algorithm> 

using namespace PhysiCell; 

std::string DC_LNP_mRNA_version = "0.1.0"; 

Submodel_Information DC_LNP_mRNA_info; 

void DC_LNP_mRNA_model_setup( void )
{

	// set up any submodels you need 
	receptor_dynamics_model_setup();

	internal_mRNA_model_setup(); 


	// set version 
	DC_LNP_mRNA_info.name = "DC LNP mRNA"; 
	DC_LNP_mRNA_info.version = DC_LNP_mRNA_version; 
		// set functions 
	DC_LNP_mRNA_info.main_function = DC_LNP_mRNA_main_model; 
	DC_LNP_mRNA_info.phenotype_function = DC_LNP_mRNA_model_phenotype_func; 
	DC_LNP_mRNA_info.mechanics_function = NULL; 	
	
	// what microenvironment variables do you need 	
	// what cell variables and parameters do you need? 
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	DC_LNP_mRNA_info.register_model(); 

	// set functions for the corresponding cell definition 
	Cell_Definition* pCD;
	pCD = find_cell_definition( "DC" ); 
	pCD->functions.update_phenotype = DC_LNP_mRNA_info.phenotype_function;
	
	// temporary ignore them now !!!!!!!!!!!!!!!!!!
	// pCD->functions.custom_cell_rule = DC_LNP_mRNA_info.mechanics_function; 
	// pCD->functions.contact_function = epithelium_contact_function; 
	
	return; 
}

void DC_LNP_mRNA_model_phenotype_func( Cell* pCell, Phenotype& phenotype, double dt )
{
	// receptor dynamics 
	// requires faster time scale - done in main function 
	
	internal_mRNA_dynamics_info.phenotype_function(pCell,phenotype,dt); 
	// internal_mRNA_model(pCell,phenotype,dt);

	return; 
}


void introduce_DC_cells( void )
{	

	static bool immune_cells_introduced = false; 
	static double immune_activation_time = 3*60-1; 

	if( PhysiCell_globals.current_time > immune_activation_time - 0.01*diffusion_dt && immune_cells_introduced == false )
	{
		std::cout << "Therapy activated!" << std::endl << std::endl; 
		immune_cells_introduced = true; 
	}
	else
	{
		return;
	}

	double tumor_radius = parameters.doubles( "tumor_radius" ) ;  
	
	// if this goes wackadoodle, choose 250 
	if( tumor_radius < 250.0 )
	{ tumor_radius = 250.0; }
	
	std::cout << "current tumor radius: " << tumor_radius << std::endl; 
	
	// now seed immune cells 
	
	int number_of_immune_cells = 100; 
		// parameters.ints("number_of_immune_cells"); // 7500; // 100; // 40; 
	double radius_inner = tumor_radius + 50; 
		// parameters.doubles("initial_min_immune_distance_from_tumor"); // 30.0; // 75 // 50; 
	double radius_outer = radius_inner + 100; 
		// parameters.doubles("thickness_of_immune_seeding_region"); // 75.0; // 100; // 1000 - 50.0; 
	
	double mean_radius = 0.5*(radius_inner + radius_outer); 
	double std_radius = 0.33*( radius_outer-radius_inner)/2.0; 

	static int LNP_index = microenvironment.find_density_index( "Lipid-nanoparticles" ); 
	
	for( int i=0 ;i < number_of_immune_cells ; i++ )
	{
		double theta = UniformRandom() * 6.283185307179586476925286766559; 
		// double phi = acos( 2.0*UniformRandom() - 1.0 );  
		double radius = NormalRandom( mean_radius, std_radius ); 
		
		Cell* pCell = create_cell( get_cell_definition("DC") ); 
		static int LNP_endosome = pCell->custom_data.find_variable_index( "uptaken_LNP" ); 

		pCell->assign_position( radius*cos(theta), radius*sin(theta), 0.0); // if 3D: radius*cos(theta)*sin(phi)

		// may simplify the internalization dynamics by remove receptor dynamics, and insert LNP directly in the beginning
		pCell->phenotype.molecular.internalized_total_substrates[ LNP_index ] = 5; // actually don't use this func now
		pCell->custom_data[LNP_endosome] = 5; 
	}
	
	return; 
}



void inject_mRNA_loaded_LNP( void )  
{

	static double clearance_rate = parameters.doubles( "clearance_rate" );   //  0.0005; 
	static double doses = parameters.doubles( "mRNA_doses" );  // 2    inject doses

	static double therapy_dt = 3; // update therapy every 3 minutes
	static double next_therapy_time = 0.0; 
	
	static double tolerance = 0.01 * diffusion_dt; 	
	static int therapy_update_time_multiplier = 0; 

	// DON'T CHECK FOR THERAPY EVERY 0.01 MINUTES! 
	// Let's check it every 3 minutes instead
	
	if( PhysiCell_globals.current_time <= next_therapy_time - tolerance )
	{ return; }
	else
	{ 
		therapy_update_time_multiplier++;
		next_therapy_time = therapy_update_time_multiplier*therapy_dt;     
	}
 
    static double start_inject = 60*parameters.doubles( "start_inject_time" ); // when starting to inject mRNA     
	static int dose_dt = 24*60*parameters.doubles( "injection_frequency" );   // inject LNP every n days 
	static int next_dose_time = 0;
	static int dose_update = 0;

	if ( PhysiCell_globals.current_time <= start_inject - tolerance )
	{ return; }

	if ( PhysiCell_globals.current_time > next_dose_time - tolerance )
	{
		dose_update++;
		next_dose_time = start_inject + dose_update*dose_dt; 
	}

	double elapsed_therapy_time = PhysiCell_globals.current_time - (next_dose_time - dose_dt); 
    
	static int LNP_index = microenvironment.find_density_index( "Lipid-nanoparticles" ); 
	// If therapy is started, this time is non-negative. 
	if( elapsed_therapy_time >= -tolerance )
	{
		// compute the current boundary / blood value, 
		// based on the model that the drug is globally cleared (e.g., in the renal system)
		double dose = doses * exp( -clearance_rate * elapsed_therapy_time ); 
		
		// make sure the kth condition is set to be Dirichlet ( turn on boundry condition !!!!!! )
		microenvironment.set_substrate_dirichlet_activation(LNP_index, true); 
		
		#pragma omp parallel for 
		for( int i=0 ; i < microenvironment.mesh.voxels.size() ;i++ )
		{
			if( microenvironment.mesh.voxels[i].is_Dirichlet == true )
			{ microenvironment.update_dirichlet_node( i, LNP_index, dose ); }
		}		
	}
		
	
	return; 

}


void DC_LNP_mRNA_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int DC_type = get_cell_definition( "DC" ).type; 
	
	// bookkeeping -- find microenvironment variables we need

	// bookkeeping -- find custom data we need 
	static double DCprob = parameters.doubles( "DC_leave_prob" ); 
	extern double DCAMOUNT; //declare existance of counter
	// do nothing if dead 
	if( phenotype.death.dead == true )
	{ return; } 

	// if not DC, do nothing 
	if( pCell->type != DC_type )
	{ return; } 
	
	// (Adrianne) if DC is already activated, then check whether it leaves the tissue
	if( pCell->custom_data["activated_immune_cell"] >  0.5 && UniformRandom() < DCprob)
	{
		// (Adrianne) DC leaves the tissue and so we lyse that DC
		std::cout<<"DC leaves tissue"<<std::endl;
		pCell->lyse_cell(); 
		#pragma omp critical 
		{ DCAMOUNT++; } // add one	
		return;
		
	}
	
	return; 
}


void DC_LNP_mRNA_main_model( double dt )
{
	extern double DCAMOUNT;
	extern std::vector<int>history;
	DCAMOUNT=0;
	
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{ DC_LNP_mRNA_model( pC, pC->phenotype , dt ); }
	}
	std::rotate(history.rbegin(),history.rbegin()+1,history.rend());
	history.front() = DCAMOUNT;
	
	/* std::copy(history.begin(), history.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl; 
	 */
	return; 
}
