%% Compare mSim with analytical solution
%% Overview
% In this tutorial we will build a simple 3D cylindrical model with one
% well at the center and we will compare it with the Thiem, 1906 analytical
% solution.
%% Analytical solution 
% The analytical solution for radial flow to a well where the head h0 is known
% at some point at r0 distance close to well is given by :
%%
% $h(r)=\frac{Q}{2 \pi T}ln\frac{r}{r_0} + h_0$ 
%%%
% where Q is the pumping rate, T is the transmissivity.
%%% 
% To visualize the above equation we will generate random points within a
% disk of radious r0, triangulate them and plot
%%
% Lets define the inputs first
r0 = 2000; % [m]
h0 = 30; % [m]
Q = 500; % [m^3/day]
T = 50; % [m^2/day]
%%%
% Generate points on disk (This is actually a rectangle)
[disk_pntsX, disk_pntsY] = meshgrid(-r0:10:r0);
%%%
% Evaluate the Thiem equation
Thiem = @(r)(Q/(2*pi*T))*log(r/r0)+h0;
disk_pntsR = sqrt(disk_pntsX.^2 + disk_pntsY.^2);
%%%
% and plot them
mesh(disk_pntsX, disk_pntsY, Thiem(disk_pntsR),...
    'edgecolor','none','FaceColor','interp','FaceLighting','phong')
camlight right
view(-76, 13)
drawnow
%% Numerical solution with mSim
%%
% *Mesh generation*
%%%
% For the numerical solution we have first to generate a mesh.
%%%
% First we are going to generate a cylindrical polygon with the help of the
% parametric equation of a circle x = r*cos(t), y = r*sin(t) where t spans
% from 0 to 2*pi for a full circle
circle_bnd = [r0*cos(0:0.2:2*pi)' r0*sin(0:0.2:2*pi)'];
%%%
% Prepare the shapefile structure that is needed for the mesh
% generation for the domain
circle_domain.Geometry = 'Polygon';
circle_domain.X = [circle_bnd(:,1)' circle_bnd(1,1) nan]; 
circle_domain.Y = [circle_bnd(:,2)' circle_bnd(1,2) nan]; 
%%%
% and the well
well.Geometry = 'Point';
well.X = 0;
well.Y = 0;
well.DistMin = 1;
well.DistMax = r0;
well.LcMin = 1;
well.LcMax = 400;
%%%
% Construct the CSG object and generate the mesh
circle_dom = CSGobj_v2(2,1,10,10,1);
circle_dom = circle_dom.readshapefile(circle_domain);
circle_dom = circle_dom.readshapefile(well);
%%%
% define the mesh options. Here the 400 m is roughly the distance between
% the boundary points when r0 = 2000 and 0.2 discretization of the
% parametric variable that generated the circular points.
meshopt = msim_mesh_options;
meshopt.lc_gen = 400;
meshopt.embed_points = 1;
%%%
% run gmsh
circle_dom.writegeo('thiem_example', meshopt);
gmsh_path = '/home/giorgk/Downloads/gmsh-4.4.1-Linux64/bin/gmsh';
circle_dom.runGmsh('thiem_example',gmsh_path,[]);
%%
% read mesh into matlab
[p, MSH]=read_2D_Gmsh('thiem_example', 0, 0);
%% 
% *Boundary conditions*
%%%
% First we will identify the boundary the mesh nodes and assign a constant head of h0
id_bnd = find(sqrt(sum(p.^2,2)) > 1999);
CH = [id_bnd h0*ones(length(id_bnd),1)];
%%
% *Well flux*
%%%
% For the well we have to find first which node corresponds to the well
id_well = find(sqrt(sum(p.^2,2)) < 0.1);
FLUX_well = [id_well -Q];

%%
% *Assemble*
%%%
% To assemble we need the transmissivity, which we defined as constant on
% the mesh nodes
Tnd = T*ones(size(p,1),1);

%%%
% Set up the standard simulation options and assemble
simopt.dim=2;
simopt.el_type='triangle';
simopt.el_order='linear';
[Kglo, H]= Assemble_LHS(p, MSH(3,1).elem(1,1).id, Tnd , CH, [], simopt);
%%%
% Convert Point flux matrix to sparse vector
F=sparse(FLUX_well(:,1),1,FLUX_well(:,2),length(H),1);
%%%
% *Solve*
Hnumerical=solve_system(Kglo,H,F);
%% Compare the solutions
% To compare the the numerical solution of mSim with the analytical we have
% to evaluate the analytical over the points of the numerical.
Hanalytical = Thiem(sqrt(sum(p.^2,2)));
%%%
% The numerical solution contains the center of well where the
% analytical solution cannot be defined. Therefore we will make them nan to
% both solutions. 
Hanalytical(id_well,1) = nan;
Hnumerical(id_well,1) = nan;
%%%
% now we plot the relative error with respect the analytical solution
trisurf(MSH(3,1).elem(1,1).id, p(:,1), p(:,2), 100*(Hnumerical - Hanalytical)./Hanalytical)
title('Relative error wrt analytical solution')
zlabel('%')
view(0,0)
colormap cool
drawnow
%%
% We can confirm that the numerical solution follows the analytical very
% close as the error is less that 1 percent.




