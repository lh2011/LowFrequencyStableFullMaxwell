%% Problem specification
% consts
eps0 = 8.854187e-12;
sig0 = 6e7;
nu0 = 1/(4*pi)*1e7;
nuR = 1;

% geo
dimCond = 0.03; % m
numLoops = 3;
[geo, condVec, patchIndexHigh, patchIndexLow,epsVec,regionVec] = cnstrct_coil_with_ring(dimCond);

% refinement
%deg = 2;
%subs = 2;

% export
exportToParaview = true;

% frequency
%f = 150;
omega = 2*pi*f;                             

% metrial distribution
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
highBndSide = intersect( find([boundaries.patches] == patchIndexHigh), find([boundaries.faces] == 1));
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

%disp("Created EQS System")

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

%disp("Created Maxwell System")

%% MAXWELL: Apply Neumann boundary conditions
maxwell_rhs = apply_NeumannBC_maxwell(maxwell_rhs, boundaries, g_maxwell, space, msh, nmn_sides_maxwell);

%% MAXWELL: Apply Dirichlet boundary conditions (L2-projection)
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h_maxwell, drchlt_sides_maxwell);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
boundary_dofs = space.boundary.dofs;
nmn_dofs = setdiff(boundary_dofs, drchlt_dofs);

%% GAUGING: find conductor nodes, edges and boundaries
conductor_elems = find(sigma_coeffs > 0);

conductor_int_idx = [];
air_int_idx = [];

conductor_int_nodes = [];
air_int_nodes = [];

for i = 1:numel(sigma_coeffs)
    if any(conductor_elems(:) == i)
        conductor_int_idx = union(conductor_int_idx, space.gnum{i}(1:space.sp_patch{i}.ndof));
        conductor_int_nodes = union(conductor_int_nodes, space_eqs.gnum{i}(1:space_eqs.sp_patch{i}.ndof));
    else
        air_int_idx = union(air_int_idx, space.gnum{i}(1:space.sp_patch{i}.ndof));
        air_int_nodes = union(air_int_nodes, space_eqs.gnum{i}(1:space_eqs.sp_patch{i}.ndof));
    end
end

% calculate interface dofs through intersection between other domains
interface_idx = intersect(air_int_idx, conductor_int_idx);
interface_nodes = intersect(air_int_nodes, conductor_int_nodes);

% remove interface from both regions
air_nodes = setdiff(air_int_nodes, interface_nodes);
conductor_nodes = setdiff(conductor_int_nodes, interface_nodes);
air_idx = setdiff(air_int_idx, interface_idx);
conductor_idx = setdiff(conductor_int_idx, interface_idx);

% calculate dir node dofs for maxwell Bnd
[~, drchlt_dofs_eqs_maxwellBnd] = sp_drchlt_l2_proj (space_eqs, msh_eqs, @(x,y,z,ind) zeros(size(x)), drchlt_sides_maxwell);

% remove dir nodes
conductor_nodes = setdiff(conductor_nodes, drchlt_dofs_eqs_maxwellBnd);
air_nodes = setdiff(air_nodes, drchlt_dofs_eqs_maxwellBnd);
interface_nodes = setdiff(interface_nodes, drchlt_dofs_eqs_maxwellBnd);

%disp("calculated dofs for different regions")

%% GAUGING: construct tree with priorised conductor interfaces
[tree_dofs, cotree_dofs] = cotree_multipatch(space_eqs, space, boundaries, drchlt_sides_maxwell, interface_idx, geometry, 2);

%% GAUGING: Modify tree/cotree to split into both regions
% helpers
conductor_real_dofs = setdiff(conductor_idx, drchlt_dofs);
air_real_dofs = setdiff(air_idx, drchlt_dofs);
int_real_dofs = setdiff(interface_idx, drchlt_dofs);
% TCG splitting in all regions
cotree_dofs_air = intersect(cotree_dofs, air_real_dofs);
tree_dofs_air = intersect(tree_dofs, air_real_dofs);
cotree_dofs_cond = intersect(cotree_dofs, conductor_real_dofs);
tree_dofs_cond = intersect(tree_dofs, conductor_real_dofs);
cotree_dofs_int = intersect(cotree_dofs, int_real_dofs);
tree_dofs_int = intersect(tree_dofs, int_real_dofs);

% add interface dofs into cond dofs
cotree_dofs_cond = union(cotree_dofs_cond, cotree_dofs_int);
tree_dofs_cond = union(tree_dofs_cond, tree_dofs_int);
conductor_nodes = union(conductor_nodes, interface_nodes);
 
%disp('Constructed and split tree')

%% GAUGING: Calculate Gauging matrix
% calculate Coulomb Gauging matrix 
B_air = calculateGaugingMatrix(space_eqs, space, msh, epsilon_coeffs);
B_cond = G;

% split into tree and cotree contributions
B_air_c = B_air(air_nodes,cotree_dofs_air);
B_air_t = B_air(air_nodes,tree_dofs_air);

B_cond_c = B_cond(conductor_nodes,cotree_dofs_cond);
B_cond_t = B_cond(conductor_nodes,tree_dofs_cond);

%disp('Calculated Gauging matrix')

%% GAUGING: Build the linear system
% original system matrix
K = stiff_mat + mass_mat;

% simpler name conventions for the assembly of the bigger systems
Cc = cotree_dofs_cond;
Tc = tree_dofs_cond;
Ta = tree_dofs_air;
Ca = cotree_dofs_air;
Nc = conductor_nodes;
Na = air_nodes;
N = union(Na, Nc);

T = union(tree_dofs_air, tree_dofs_cond);
C = union(cotree_dofs_air, cotree_dofs_cond);

% scaling factors
gamma = (1 + omega) * (max(sigma_coeffs) + artificSig)/max(epsilon_coeffs);
beta = (1 + omega);

% matrix assembly
S_c = [ beta*B_cond_c,                              beta*B_cond(conductor_nodes, cotree_dofs_air); 
        gamma*B_air(air_nodes, cotree_dofs_cond)    gamma*B_air_c                               ];
S_t = [ beta*B_cond_t,                              beta*B_cond(conductor_nodes, tree_dofs_air);
        gamma*B_air(air_nodes, tree_dofs_cond),    gamma*B_air_t                               ];

K_cc = [K(Cc,Cc), K(Cc, Ca); K(Ca, Cc), K(Ca,Ca)];
K_ct = [K(Cc,Tc), K(Cc, Ta); K(Ca, Tc), K(Ca,Ta)];

mat = [K_cc, K_ct;
        S_c, S_t];

%disp('Calculated system matrices')

%% MAXWELL: solve gauged system (build Dir. contr.)

% naming convention for strong dirichlet boundaries
Cd = intersect(cotree_dofs, drchlt_dofs);
Td = intersect(tree_dofs, drchlt_dofs);

Cdc = intersect(conductor_int_idx, Cd);
Cda = intersect(air_idx, Cd);
Tdc = intersect(conductor_int_idx, Td);
Tda = intersect(air_idx, Td);

%% NOT REALLY CORRECT BUT NOT MATTER ATM %%
S_c_intDir = blkdiag(beta*B_cond(Nc,Cdc), gamma*B_air(Na,Cda));
S_t_intDir = blkdiag(beta*B_cond(Nc,Tdc), gamma*B_air(Na,Tda));
%% --
rhs = [ maxwell_rhs(Cc) - K(Cc,Cdc) * u(Cdc) - K(Cc,Cda) * u(Cda) - K(Cc, Tdc) * u(Tdc) - K(Cc, Tda) * u(Tda);
        maxwell_rhs(Ca) - K(Ca,Cdc) * u(Cdc) - K(Ca,Cda) * u(Cda) - K(Ca, Tdc) * u(Tdc) - K(Ca, Tda) * u(Tda);
        zeros(numel(T),1) - S_c_intDir * u(Cd) - S_t_intDir * u(Td)];

% solve sys
sol = mat \ rhs;

% reorder
u(Cc) = sol(1:numel(Cc));
off = numel(Cc);
u(Ca) = sol(off+(1:numel(Ca)));
off = off + numel(Ca);
u(Tc) = sol(off+(1:numel(Tc)));
off = off + numel(Tc);
u(Ta) = sol(off+(1:numel(Ta)));

%disp('Solved Maxwell step')

%% post
% plotting

if exportToParaview
    output_file2 = 'vts/inductiveCouplingExampleEQS';
    output_file1 = 'vts/inductiveCouplingExampleBOTH';
    vtk_pts = {linspace(0, 1, 15), linspace(0, 1, 15), linspace(0, 1, 15)};
    scaling_coeffs = sigma_coeffs;
    sp_to_vtk_u_gradphi_mp_epsScale_material(1j*omega*u, space, phi, space_eqs, geometry, vtk_pts, output_file1, scaling_coeffs, epsilon_coeffs, sigma_coeffs, nu_coeffs, {'J'}, {'value'})
    sp_to_vtk_u_gradphi_mp_epsScale_material(zeros(size(u)), space, phi, space_eqs, geometry, vtk_pts, output_file2, scaling_coeffs, epsilon_coeffs, sigma_coeffs, nu_coeffs, {'J'}, {'value'}) 
    fprintf ('The result is saved in the files %s and %s \n \n', output_file1, output_file2);
end
% csvwrite([output_file1,'_u.csv'],u);
% csvwrite([output_file1,'_phi.csv'],phi);

function vals = eqsDirichletFun(x,y,z,ind,patchesHighBnd,patchesLowBnd)
    if ismember(ind,patchesHighBnd)
        vals = 1*ones(size(x));
        return
    end
    vals = zeros(size(x));
end



% %%%%%%%%%  mass mat contribution for every patch
% Meps = sparse(space.ndof, space.ndof);
% for i=1:numel(epsilon_coeffs)
%     Meps = Meps + epsilon_coeffs(i) * op_u_v_mp (space, space, msh, @(x,y,z) ones(size(x)), i);
% end
% %%%%%% Hdiv space %%%%%%%
% %% Discretization parameters
% degree      = deg*[1 1 1];
% regularity  = degree - 1;
% nsub        = subs*[1 1 1];
% nquad       = degree + 1;
% [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo);
% 
% npatch = numel (geometry);
% 
% % 'geo_2cubesa.txt'
% % 'geo_cube.txt'
% 
% %% create spaces and meshes for subdomains
% msh = cell (1, npatch);
% sp = cell (1, npatch);
% for iptc = 1:npatch
%     [knots, zeta] = ...
%         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
%     [knots_hdiv, degree_hdiv] = knt_derham (knots, degree, 'Hdiv');
% 
%     % Construct msh structure
%     rule      = msh_gauss_nodes (nquad);
%     [qn, qw]  = msh_set_quad_nodes (zeta, rule);
%     msh{iptc} = msh_cartesian (zeta, qn, qw, geometry(iptc));
% 
%     % Construct space structure
%     scalar_spaces = cell (msh{iptc}.ndim, 1);
%     for idim = 1:msh{iptc}.ndim
%         scalar_spaces{idim} = sp_bspline (knots_hdiv{idim}, degree_hdiv{idim}, msh{iptc});
%     end
%     sp{iptc} = sp_vector (scalar_spaces, msh{iptc}, 'div-preserving');
%     clear scalar_spaces
% end
% 
% %% create mp space and mesh
% msh_Hdiv = msh_multipatch (msh, boundaries);
% space_Hdiv = sp_multipatch (sp, msh_Hdiv, interfaces, boundary_interfaces);
% clear sp
% 
% M_Hdiv = op_u_v_mp(space_Hdiv,space_Hdiv,msh_Hdiv);
% 
% d=M_Hdiv\Meps*u;
