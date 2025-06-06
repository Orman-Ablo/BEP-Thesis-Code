function [selected_indices] = select_spots(roixy,params)
%select_global_spots Select spots for global fitting uniformly over the
%FOV. This is done by dividing the FOV in a 10x10 grid and selecting uniformly
%from each grid cell.

grid_size = 10;
nr_of_cells = grid_size*grid_size;
Ncfg_total = params.Ncfg_total;
Ncfg = params.Ncfg;
nr_of_spots_per_cell = floor(Ncfg/nr_of_cells);
selected_indices = zeros(Ncfg,1);
Ncfg_max = 1e5; % To prevent the selection process to become slow

imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
cellSizeX = imgSizeX/grid_size;
cellSizeY = imgSizeY/grid_size;

% For each spot, determine in which cell it is located.
indices_per_cell = cell(grid_size,grid_size);
for i=1:min(Ncfg_total,Ncfg_max)
    xi = ceil(roixy(1,i)/cellSizeX);
    yi = ceil(roixy(2,i)/cellSizeY);
    indices_per_cell{xi,yi} = [indices_per_cell{xi,yi}, i];
end

% Loop over the cells and select spots randomly.
jcfg_selected = 0;
for xi=1:grid_size
    for yi=1:grid_size
        
        % Handle when cell does not contain enough elements
        spot_indices = indices_per_cell{xi,yi};
        Nxy = min(numel(spot_indices),nr_of_spots_per_cell);
        jcfg_selected = jcfg_selected + Nxy;
        
        index = randperm(numel(spot_indices),Nxy);
        spots = spot_indices(index);
        selected_indices(jcfg_selected-Nxy+1:jcfg_selected) = spots;

    end
end

% Select the remaining spots randomly from the not yet chosen indices.
to_be_selected = Ncfg - jcfg_selected;
if to_be_selected>0
    not_selected = setdiff(1:Ncfg_total,selected_indices);
    index = randperm(numel(not_selected),to_be_selected);
    selected_indices(jcfg_selected+1:end) = not_selected(index);
end

end
