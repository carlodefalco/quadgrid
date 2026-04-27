
%% Create input files in JSON format for use
%% with the quadgrid library

%% properties of the grid
hy = 0.025;
hx = 0.025;
nx = 2/hx;
ny = 2/hx;
px=3;
py=3;
rx=2;
ry=2;

in.grid_properties = struct ('hx', hx, 'hy', hy, 'nx', nx, 'ny', ny,'rx',rx,'ry',ry,'px',px,'py',py);
num_global_nodes = ((nx - 1) * (px - rx) + (px + 1))*((ny - 1) * (py - ry) + (py + 1));

%% compute arrays to be used to define fields on the grid

x = 0 : hx : nx*hx; %% as we want to compute the velocity field
y = 0 : hy : ny*hy; %% as a function of the coordinates, we need
%%                  %% to compute the coordinates here, they need
%%                  %% not be saved in the JSON file though.




fy = @(x) 1e-3*(2*pi*(x-1));

z1= fy(x);
pp=spline(x,z1);
sconds=3*ones(size(pp.breaks,2)-2,1);% uguale al numero di interior breakpts, è il numero di condizioni di regolarità in ogni nodo.
sp_y=fn2fm(pp,'B-',sconds);

% We should multiply the coeffs of Bsplines in x with the ones in y, but since the function is constant in y they are 
% alwais 1 (because of partition of unity),so we need only the coeffs in x.


fx= @(y) 1e-3*(2*pi*(1-y));
z2=fx(y);
pp=spline(y,z2);
sconds=3*ones(size(pp.breaks,2)-2,1);
sp_x=fn2fm(pp,'B-',sconds);

% NB If the function were not separable in x and y the global coeffs would not have been the multiplication of the mono-dimensional coeffs
[vy,vx]=meshgrid(sp_y.coefs,sp_x.coefs);
vy=vy(:);
vx=vx(:);

%[x, y] = meshgrid (x, y); %% the ordering in quadgrid is explicitely
%%                        %% intended to be compatible to that in
%%                        %% Matlab/Octave

%tmp = 1e-3*(2*pi*(1-y));
%vx=tmp(:); %% grid arrays must be flattened
%tmp = 1e-3*(2*pi*(x-1));
%vy=tmp(:);

in.grid_vars = struct ('vx', vx, 'vy', vy);% 'rho', zeros(num_global_nodes, 1));

%% generate particles and assign position and properties

in.num_particles = 2.e5;
in.x             = randn(in.num_particles, 1) * .2/6 + 1.5;
in.y             = randn(in.num_particles, 1) * .2/6 + 1.5;
tmp  = [0:numel(in.x)-1];
in.iprops.label  = tmp(:);

in.dprops = struct ('VX', zeros(in.num_particles, 1), 'VY', zeros(in.num_particles, 1));
% 'M',  2*pi*(.2/6)^2*ones(in.num_particles, 1)/in.num_particles,

%% grid and particles are read separately but can be stored in the same file
str = jsonencode (in);
fid = fopen ('velocity.json', 'w');
fwrite (fid, str, 'char');
fclose (fid);



% capire cosa fa per input in codice mpm
