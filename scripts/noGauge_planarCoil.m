% consts
eps0 = 8.854187e-12;
sig0 = 6e7;
nu0 = 1/(4*pi)*1e7;
nuR = 1;

% geo
dimCond = 0.003; % m
numLoops = 3;
type = 1; 
[geo, condVec, patchIndexHigh, patchIndexLow,epsVec] = cnstrct_paper_example_config(dimCond, numLoops, type);

% refinement
%deg = 2;
%subs = 2;

% export
exportToParaview = true;

% frequency
%f = 150e0;
omega = 2*pi*f;                             

sigma_coeffs = condVec * sig0;
epsilon_coeffs = epsVec * eps0;
nu_coeffs = ones(size(condVec)) * nu0;

% RHS functions
f = @(x, y, z) cat(1, ...
    zeros ([1, size(x)]), ...
    zeros ([1, size(x)]), ...
    zeros ([1, size(x)]));

rho = @(x,y,z) zeros(size(x));

%% EQS: create geometry, space and mesh
[space_eqs, msh_eqs, geometry, boundaries] = create_space_mesh_EQS(geo, deg, subs);

%% Boundaries
% right bottom conductor patch, positive x direction -> 2
highBndSide = intersect( find([boundaries.patches] == patchIndexHigh), find([boundaries.faces] == 2));
% middle conductor patch, z direction -> 5/6
lowBndSide = intersect( find([boundaries.patches] == patchIndexLow), union(find([boundaries.faces] == 5),find([boundaries.faces] == 6)));

% eqs sides
drchlt_sides_eqs    = [highBndSide, lowBndSide];
nmn_sides_eqs       = setdiff(1:numel([boundaries.patches]), drchlt_sides_eqs);
% maxwell sides
nmn_sides_maxwell_global    = [];
drchlt_sides_maxwell_global  = setdiff(1:6, nmn_sides_maxwell_global);

% Neumann functions
g_maxwell = @(x, y, z, ind) cat(1, ...
    zeros ([1, size(x)]), ...
    zeros ([1, size(x)]), ...
    zeros ([1, size(x)]));

g_eqs = @(x,y,z,ind) zeros(size(x));

% Dirichlet functions
h_maxwell = @(x, y, z, ind) cat(1, ...            
                    zeros([1, size(x)]), ...
                    zeros([1, size(x)]), ...
                    zeros([1, size(x)]));

% proper function for boundary
h_eqs = @(x,y,z,ind) eqsDirichletFun(x,y,z,ind,highBndSide,lowBndSide);

%% EQS: assemble matrices
[eqs_mat, eqs_rhs] = create_system_matrices_EQS(space_eqs, msh_eqs, epsilon_coeffs, sigma_coeffs, rho, omega);

%disp('Created EQS System')

%% EQS: Apply Neumann boundary conditions
eqs_rhs = apply_NeumannBC_eqs(eqs_rhs, boundaries, g_eqs, space_eqs, msh_eqs, nmn_sides_eqs);

%% EQS: Apply Dirichlet boundary conditions (L2-projection)
phi = zeros (space_eqs.ndof, 1);
[phi_drchlt, drchlt_dofs_eqs] = sp_drchlt_l2_proj (space_eqs, msh_eqs, h_eqs, drchlt_sides_eqs);
phi(drchlt_dofs_eqs) = phi_drchlt;
int_dofs_eqs = setdiff (1:space_eqs.ndof, drchlt_dofs_eqs);

%% EQS: Solve the linear system
eqs_rhs(int_dofs_eqs) = eqs_rhs(int_dofs_eqs) - eqs_mat(int_dofs_eqs, drchlt_dofs_eqs)*phi_drchlt;
phi(int_dofs_eqs) = eqs_mat(int_dofs_eqs, int_dofs_eqs) \ eqs_rhs(int_dofs_eqs);

%disp('Solved EQS')

%% MAXWELL: create geometry, space and mesh
[space, msh, ~, boundaries] = create_space_mesh_maxwell(geo, deg, subs);

% boundary sides
nmn_sides_maxwell = [];
for ni = nmn_sides_maxwell_global
    nmn_sides_maxwell   = union(nmn_sides_maxwell, find ([boundaries.faces] == ni));
end
drchlt_sides_maxwell  = setdiff(1:numel(boundaries), nmn_sides_maxwell);

%% MAXWELL: create system matrices
[mass_mat, stiff_mat, f_maxwell] = create_system_matrices_maxwell(space, msh, epsilon_coeffs, nu_coeffs, sigma_coeffs, f, omega);
% EQS modification
G = calculateGaugingMatrix(space_eqs, space, msh, sigma_coeffs + 1i * omega * epsilon_coeffs);
maxwell_rhs = f_maxwell - G' * phi;

%disp('Created Maxwell System')

%% MAXWELL: Apply Neumann boundary conditions
maxwell_rhs = apply_NeumannBC_maxwell(maxwell_rhs, boundaries, g_maxwell, space, msh, nmn_sides_maxwell);

%% MAXWELL: Apply Dirichlet boundary conditions (L2-projection)
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h_maxwell, drchlt_sides_maxwell);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
boundary_dofs = space.boundary.dofs;
nmn_dofs = setdiff(boundary_dofs, drchlt_dofs);

%% MAXWELL: Solve system
K = stiff_mat + mass_mat;
u(int_dofs) = K(int_dofs, int_dofs) \ (maxwell_rhs(int_dofs)); % No dirichlet here because 0

% plotting
if exportToParaview
    output_file1 = 'vts/planarCoilExampleNoGaugeB';
    output_file2 = 'vts/planarCoilExampleNoGaugeE';
    vtk_pts = [20 20 20]; {linspace(0, 1, 15), linspace(0, 1, 15), linspace(0, 1, 15)};
    sp_to_vtk(u, space, geometry, vtk_pts,output_file1,{'B'},{'curl'})
    scaling_coeffs = epsilon_coeffs;
    sp_to_vtk_u_gradphi_mp_epsScale_material(1j*omega*u, space, phi, space_eqs, geometry, vtk_pts, output_file2, scaling_coeffs, epsilon_coeffs, sigma_coeffs, nu_coeffs, {'J'}, {'value'})
    fprintf ('The result is saved in the files %s and %s \n \n', output_file1, output_file2);
end

csvwrite([output_file1,'_u.csv'],u);
csvwrite([output_file1,'_phi.csv'],phi);

%keyboard
function vals = eqsDirichletFun(x,y,z,ind,patchesHighBnd,patchesLowBnd)
    if ismember(ind,patchesHighBnd)
        vals = 1*ones(size(x));
        return
    end
    vals = zeros(size(x));
end
