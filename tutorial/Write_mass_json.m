
clear;
close all;

hy = 0.025;
hx = 0.025;
nx = 2/hx;
ny = 2/hx;
% Degree
px=1;
py=1;
% Regularity
rx=0;
ry=0;


in.grid_properties = struct ('hx', hx, 'hy', hy, 'nx', nx, 'ny', ny,'rx',rx,'ry',ry,'px',px,'py',py);
num_global_nodes = ((nx - 1) * (px - rx) + (px + 1))*((ny - 1) * (py - ry) + (py + 1));


in.grid_vars = struct ('Mn', zeros(num_global_nodes, 1));% credo basti il gridname essendo default init


in.num_particles = 2.e5;
in.x             = randn(in.num_particles, 1) * .2/6 + 1.5;
in.y             = randn(in.num_particles, 1) * .2/6 + 1.5;
tmp  = [0:numel(in.x)-1];
in.iprops.label  = tmp(:);



masses=randn(in.num_particles, 1)*0.2+1;
masses(masses<0)=0.5;
in.dprops = struct ('mp', masses);


str = jsonencode (in);
fid = fopen ('mass.json', 'w');
fwrite (fid, str, 'char');
fclose (fid);



