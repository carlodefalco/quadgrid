clear;
close all;


file = '/home/leonsalvu/Tesi/quadgrid/tutorial/particle.450.csv';

T = readtable(file);   % legge il CSV (anche con header)
xp = T{:,1};           
yp = T{:,2};           
% mp = T{:,3};
mp=ones(size(xp));

% figure;
% axis equal
% scatter3(xp,yp,hp);

figure;
scatter3(xp, yp, mp,20,mp, 'filled');
xlabel('xp'); ylabel('yp'); zlabel('mp');
colorbar; view(3);% axis equal;

xlim([0 2]);
ylim([0 2]);

% 
% tri = delaunay(xp, yp);
% figure;
% trisurf(tri, xp, yp, hp, hp, 'EdgeColor', 'none');
% xlabel('x'); ylabel('y'); zlabel('z');
% colorbar; axis tight; view(3); camlight; lighting gouraud;
