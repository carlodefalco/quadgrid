%% Create input files in JSON format for use
%% with the quadgrid library

%% properties of the grid

hx = hy = 0.025;
nx = ny = 2/hx;

in.grid_properties = struct ('hx', hx, 'hy', hy, 'nx', nx, 'ny', ny);
num_global_nodes = (in.grid_properties.nx+1)*(in.grid_properties.ny+1);

%% compute arrays to be used to define fields on the grid

x = 0 : hx : nx*hx; %% as we want to compute the velocity field
y = 0 : hy : ny*hy; %% as a function of the coordinates, we need
%%                  %% to compute the coordinates here, they need
%%                  %% not be saved in the JSON file though.

[x, y] = meshgrid (x, y); %% the ordering in quadgrid is explicitely
%%                        %% intended to be compatible to that in
%%                        %% Matlab/Octave

vx = 1e-3*(2*pi*(1-y))(:); %% grid arrays must be flattened
vy = 1e-3*(2*pi*(x-1))(:);

in.grid_vars = struct ('betax', vx, 'betay', vy, 'rho', zeros(num_global_nodes, 1), 'Jdrift_x', zeros(num_global_nodes,1), 
                       'Jdrift_y', zeros(num_global_nodes,1), 'Jdiff_x', zeros(num_global_nodes,1), 'Jdiff_y', zeros(num_global_nodes,1));

%% generate particles and assign position and properties

in.num_particles = 2.e5;
in.x             = randn(in.num_particles, 1) * .2/6 + 1.5;
in.y             = randn(in.num_particles, 1) * .2/6 + 1.5;
in.iprops.label  = [0:numel(in.x)](:);

%% 1. Parameters for the distribution
sig = 0.2/6; 

%% 2. Calculate squared distance from center for every particle
%% dist_sq will be a vector of size (200000, 1)
dist_sq = (in.x - 1.5).^2 + (in.y - 1.5).^2;

%% 3. Calculate rho for each particle (Peak density = 1.0)
%% This creates the (Np, 1) vector you need.
rho_vector = exp(-dist_sq / (2 * sig^2));

in.dprops = struct ('M',  2*pi*(.2/6)^2*ones(in.num_particles, 1)/in.num_particles, 'BETAx', zeros(in.num_particles, 1),
                    'BETAy', zeros(in.num_particles, 1), 'VX', zeros(in.num_particles, 1), 'VY', zeros(in.num_particles, 1),
                    'zero', zeros(in.num_particles, 1) );


%% grid and particles are read separately but can be stored in the same file
str = jsonencode (in);
fid = fopen ('driftdiffusion.json', 'w');
fwrite (fid, str, 'char');
fclose (fid);
