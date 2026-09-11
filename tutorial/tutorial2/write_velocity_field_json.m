
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

in.grid_vars = struct ('vx', vx, 'vy', vy, 'rho', zeros(num_global_nodes, 1));

%% generate particles and assign position and properties

in.num_particles = 2.e5;
in.x             = randn(in.num_particles, 1) * .2/6 + 1.5;
in.y             = randn(in.num_particles, 1) * .2/6 + 1.5;
in.iprops.label  = [0:numel(in.x)](:);

in.dprops = struct ('M',  2*pi*(.2/6)^2*ones(in.num_particles, 1)/in.num_particles,
		    'VX', zeros(in.num_particles, 1), 'VY', zeros(in.num_particles, 1));


%% grid and particles are read separately but can be stored in the same file
str = jsonencode (in);
fid = fopen ('velocity.json', 'w');
fwrite (fid, str, 'char');
fclose (fid);
