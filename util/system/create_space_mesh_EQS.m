function [space, msh, geometry, boundaries] = create_space_mesh_EQS(geo_name, deg, nsubs)
%% Discretization parameters
degree      = deg*[1 1 1];
regularity  = degree - 1;
nsub        = nsubs*[1 1 1];
nquad       = degree + 1;

[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

% 'geo_2cubesa.txt'
% 'geo_cube.txt'

%% create spaces and meshes for subdomains
msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch

% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
end

%% EQS: create mp space and mesh
msh = msh_multipatch (msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
clear sp


end