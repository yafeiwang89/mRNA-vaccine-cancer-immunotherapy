


all_virion = sum( sum(MCDS.continuum_variables(1).data(:,:,1)) )*8000;

index = find(MCDS.discrete_cells.metadata.type==6);  % DC cell

b = MCDS.discrete_cells.custom.uptaken_LNP(index);

c = MCDS.discrete_cells.custom.free_LNP(index);

d = MCDS.discrete_cells.custom.mRNA(index);

e = MCDS.discrete_cells.custom.antigenic_epitope(index);

index11 = find(MCDS.discrete_cells.custom.TCell_contact_time > 0);  % DC cell
