


all_virion = sum( sum(MCDS.continuum_variables(1).data(:,:,1)) )*8000;

index = find(MCDS.discrete_cells.metadata.type==6);  % DC cell

b = MCDS.discrete_cells.custom.uptaken_LNP(index);

c = MCDS.discrete_cells.custom.assembled_virion(index);