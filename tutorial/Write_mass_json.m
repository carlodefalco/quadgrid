
clear;
close all;

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


in.grid_vars = struct ('Mn', zeros(num_global_nodes, 1));% credo basti il gridname essendo default init


in.num_particles = 2.e5;
in.x             = randn(in.num_particles, 1) * .2/6 + 1.5;
in.y             = randn(in.num_particles, 1) * .2/6 + 1.5;
tmp  = [0:numel(in.x)-1];
in.iprops.label  = tmp(:);

in.dprops = struct ('mp',  2*pi*(.2/6)^2*ones(in.num_particles, 1)/in.num_particles);


str = jsonencode (in);
fid = fopen ('mass.json', 'w');
fwrite (fid, str, 'char');
fclose (fid);



