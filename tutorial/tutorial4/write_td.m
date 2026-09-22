hx = 0.01;
hy = 0.002;

Lx = 45; 
H = 0.1;

nx = 45/hx;
ny = 2*H/hy;

in.grid_properties = struct ('hx', hx, 'hy', hy, 'nx', nx, 'ny', ny);
num_global_nodes = (in.grid_properties.nx+1)*(in.grid_properties.ny+1);

x = 0 : hx : 45;
y = 0 : hy : 2*H;

[x, y] = meshgrid (x, y);

v_avg = 1000;
vx = 1.5 * v_avg * (1 - ((y - H).^2) / (H^2)) (:) ;

in.grid_vars = struct ('betax', vx, 'betay', zeros(num_global_nodes,1), 'rho', zeros(num_global_nodes, 1), 'Jdrift_x', zeros(num_global_nodes,1), 
                       'Jdrift_y', zeros(num_global_nodes,1), 'Jdiff_x', zeros(num_global_nodes,1), 'Jdiff_y', zeros(num_global_nodes,1));
l=0.01;

in.x             = [];
in.y             = [];
x_centers = zeros(nx, 1);
for i = 1:nx
    x_centers(i) = hx*(i-0.5); 
end

mu = 5.0;            % Media della gaussiana
sigma = 1.0;         % Deviazione standard (varianza = 1)

phi_Lx = 0.5 * (1 + erf((Lx - mu) / (sigma * sqrt(2))));
phi_0  = 0.5 * (1 + erf((0 - mu) / (sigma * sqrt(2))));

NA_mtot = 2.506e5 / (phi_Lx - phi_0);

N0 = zeros(nx, 1);



for i = 1:nx
    xi = x_centers(i);
    
    % Calcolo della funzione di ripartizione negli estremi dell'intervallo [xi - 0.5*l, xi + 0.5*l]
    phi_plus  = 0.5 * (1 + erf(((xi + 0.5*l) - mu) / (sigma * sqrt(2))));
    phi_minus = 0.5 * (1 + erf(((xi - 0.5*l) - mu) / (sigma * sqrt(2))));
    
    N0(i) = round( NA_mtot * (phi_plus - phi_minus) );
    
    if N0(i) > 0
        % Genera le coordinate X (tutte identiche per questa colonna)
        x_p_colonna = ones(N0(i), 1) * xi;
        
        % Genera le coordinate Y scelte randomicamente tra 0 e +2H 
        if N0(i) == 1
            y_p_colonna = 0.1; % Se c'è una sola particella, va al centro
        else
            y_p_colonna = rand(N0(i), 1) * (2 * H);
        end
        
        in.x = [in.x; x_p_colonna];
        in.y = [in.y; y_p_colonna];
    end
end

in.num_particles = length(in.x);
in.iprops.label  = [0:numel(in.x)](:);

in.dprops = struct ('M', ones(in.num_particles, 1), 'BETAx', zeros(in.num_particles, 1),
                    'BETAy', zeros(in.num_particles, 1), 'VX', zeros(in.num_particles, 1), 'VY', zeros(in.num_particles, 1),
                    'zero', zeros(in.num_particles, 1) );

str = jsonencode (in);
fid = fopen ('td.json', 'w');
fwrite (fid, str, 'char');
fclose (fid);

