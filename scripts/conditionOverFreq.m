%% consts
eps0 = 8.854187e-12;
sig0 = 6e7;
nu0 = 1/(4*pi)*1e7;
nuR = 1;

%% geo
%dimCond = 0.003; % m
%numLoops = 3;
%type = 1;
%[geo, condVec, patchIndexHigh, patchIndexLow,epsVec] = cnstrct_paper_example_config(dimCond, numLoops, type);
l = 0.22;
geo = ostrowski_cube(l,0.02);

% refinement
%deg = 2;
%subs = 2;

%% metrial distribution
%sigma_coeffs = condVec * sig0;
%epsilon_coeffs = epsVec * eps0;
%nu_coeffs = ones(size(condVec)) * nu0;
sigma_coeffs =      [0 0 0 0 5 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 5 0 0 0 0] * sig0;
epsilon_coeffs =    [5 5 5 5 5 5 5 5 5 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 5 5 5] * eps0;
nu_coeffs =         ones(size(sigma_coeffs))*nu0 * nuR;

%% RHS functions
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
drchlt_sides_eqs    = [3 5 8 10 11 13 16 18 21 36 38 41 43 44 46 49 51 54]; %[highBndSide, lowBndSide];
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
h_eqs = @(x,y,z,ind) z./l; %eqsDirichletFun(x,y,z,ind,highBndSide,lowBndSide); % z./l;

%% loop over different frequencies

condOriginal = nan(size(freqs));
condMethod1 = nan(size(freqs));
condMethod2 = nan(size(freqs));

resOrig = nan(size(freqs));
res1 = nan(size(freqs));
res2 = nan(size(freqs));
res2OLD = nan(size(freqs));

divOrig = nan(size(freqs));
div1 =  nan(size(freqs));
div2 = nan(size(freqs));
div2OLD = nan(size(freqs));

ranks = nan([3 max(size(freqs))]);

condEIGSOrig = nan(size(freqs));
condEIGS1 = nan(size(freqs));
condEIGS2 = nan(size(freqs));

eigValsOrig = [];
eigVals1 = [];
eigVals2 = [];

for iter = 1:numel(freqs)
    %% print
    fprintf("Starting condition number computation for f = %d\n", freqs(iter));
    %% omega
    omega = 2*pi*freqs(iter);

    %% EQS: create system
    [eqs_mat, eqs_rhs] = create_system_matrices_EQS(space_eqs, msh_eqs, epsilon_coeffs, sigma_coeffs, rho, omega);

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

    %% MAXWELL: create geometry, space and mesh
    [space, msh, ~, boundaries] = create_space_mesh_maxwell(geo, deg, subs);

    % boundary sides
    nmn_sides_maxwell = [];
    for ni = nmn_sides_maxwell_global
        nmn_sides_maxwell   = union(nmn_sides_maxwell, find ([boundaries.faces] == ni));
    end
    drchlt_sides_maxwell  = setdiff(1:numel(boundaries), nmn_sides_maxwell);

    %% MAXWELL: create system matrices
    [mass_mat, stiff_mat, rhs] = create_system_matrices_maxwell(space, msh, epsilon_coeffs, nu_coeffs, sigma_coeffs, f, omega);
    % EQS modification
    G = calculateGaugingMatrix(space_eqs, space, msh, sigma_coeffs + 1i * omega * epsilon_coeffs);
    rhs = rhs - G' * phi;

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

    %% GAUGING: Calculate Gauging matrix
    % scaling factors
    gamma = (1 + omega) * (max(sigma_coeffs) + artificSig)/max(epsilon_coeffs);
    beta = (1 + omega);

    % calculate Coulomb Gauging matrix
    B_air = calculateGaugingMatrix(space_eqs, space, msh, epsilon_coeffs);
    B_cond = G;

    % split into tree and cotree contributions
    B_air_c = gamma * B_air(air_nodes,cotree_dofs_air);
    B_air_t = gamma * B_air(air_nodes,tree_dofs_air);

    B_cond_c = beta * B_cond(conductor_nodes,cotree_dofs_cond);
    B_cond_t = beta * B_cond(conductor_nodes,tree_dofs_cond);

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

    % matrix assembly
    S_cOLD = blkdiag(beta*B_cond_c, gamma*B_air_c);
    S_tOLD = blkdiag(beta*B_cond_t, gamma*B_air_t);

    S_c = [ beta*B_cond_c,                              beta*B_cond(conductor_nodes, cotree_dofs_air); 
            gamma*B_air(air_nodes, cotree_dofs_cond)    gamma*B_air_c                               ];
    S_t = [ beta*B_cond_t,                              beta*B_cond(conductor_nodes, tree_dofs_air);
            gamma*B_air(air_nodes, tree_dofs_cond),     gamma*B_air_t                               ];

    K_cc = [K(Cc,Cc), K(Cc, Ca); K(Ca, Cc), K(Ca,Ca)];
    K_ct = [K(Cc,Tc), K(Cc, Ta); K(Ca, Tc), K(Ca,Ta)];

    matSys2 = [ K_cc, K_ct;
                S_c, S_t];

    matOrig = K(int_dofs, int_dofs);

    %% compute condition number of original system and method II
    condOriginal(iter) = cond(matOrig);
    condMethod2(iter) = cond(matSys2);

    % naming convention for strong dirichlet boundaries
    Cd = intersect(cotree_dofs, drchlt_dofs);
    Td = intersect(tree_dofs, drchlt_dofs);

    Cdc = intersect(conductor_int_idx, Cd);
    Cda = intersect(air_idx, Cd);
    Tdc = intersect(conductor_int_idx, Td);
    Tda = intersect(air_idx, Td);

    % matrix and rhs assembly with dirichlet contribution
    S_c_intDir = [beta*B_cond(Nc,Cdc), beta*B_cond(Nc,Cda); gamma*B_air(Na,Cdc), gamma*B_air(Na,Cda)];
    S_t_intDir = [beta*B_cond(Nc,Tdc), beta*B_cond(Nc,Tda); gamma*B_air(Na,Tdc), gamma*B_air(Na,Tda)];

    rhsSys2 = [ rhs(Cc) - K(Cc,Cdc) * u(Cdc) - K(Cc,Cda) * u(Cda) - K(Cc, Tdc) * u(Tdc) - K(Cc, Tda) * u(Tda);
        rhs(Ca) - K(Ca,Cdc) * u(Cdc) - K(Ca,Cda) * u(Cda) - K(Ca, Tdc) * u(Tdc) - K(Ca, Tda) * u(Tda);
        zeros(numel(T),1) - S_c_intDir * u(Cd) - S_t_intDir * u(Td)];

    %% build system for method I
    % precompute contribution
    B_contr_air = (B_air_t\B_air_c);
    B_contr_cond = (B_cond_t\B_cond_c);

    % gauging contribution from conductor region
    gauge_mat_cond_ct = K(cotree_dofs_cond, tree_dofs_cond) * B_contr_cond;
    gauge_mat_ac_ct = K(cotree_dofs_air, tree_dofs_cond) * B_contr_cond;
    % gauging contribution form air region
    gauge_mat_ca_ct = K(cotree_dofs_cond, tree_dofs_air) * B_contr_air;
    gauge_mat_air_ct = K(cotree_dofs_air, tree_dofs_air) * B_contr_air;
    % assembly
    K(cotree_dofs_cond, cotree_dofs_cond) = K(cotree_dofs_cond, cotree_dofs_cond) - gauge_mat_cond_ct;
    K(cotree_dofs_air, cotree_dofs_cond) = K(cotree_dofs_air, cotree_dofs_cond) - gauge_mat_ac_ct;

    K(cotree_dofs_cond, cotree_dofs_air) = K(cotree_dofs_cond, cotree_dofs_air) - gauge_mat_ca_ct;
    K(cotree_dofs_air, cotree_dofs_air) = K(cotree_dofs_air, cotree_dofs_air) - gauge_mat_air_ct;
    cotree_dofs = union(cotree_dofs_air, cotree_dofs_cond);

    %% compute cond for method I
    matSys1 = K(cotree_dofs, cotree_dofs);
    condMethod1(iter) = cond(matSys1);

    %% compute solutions
    xOrig = zeros (space.ndof, 1);
    xSys2 = zeros (space.ndof, 1);
    xSys1 = zeros (space.ndof, 1);
    xSys2OLD = zeros (space.ndof, 1);

    xOrig(int_dofs) = matOrig \ rhs(int_dofs);

    xSys2tmp = matSys2 \ rhsSys2;
    xSys2OLD(Cc) = xSys2tmp(1:numel(Cc));
    off = numel(Cc);
    xSys2OLD(Ca) = xSys2tmp(off+(1:numel(Ca)));
    off = off + numel(Ca);
    xSys2OLD(Tc) = xSys2tmp(off+(1:numel(Tc)));
    off = off + numel(Tc);
    xSys2OLD(Ta) = xSys2tmp(off+(1:numel(Ta)));

    matSys2OLD = [ K_cc,        K_ct;
                    S_cOLD,     S_tOLD];

    xSys2tmpOLD = matSys2OLD \ rhsSys2;
    xSys2(Cc) = xSys2tmpOLD(1:numel(Cc));
    off = numel(Cc);
    xSys2(Ca) = xSys2tmpOLD(off+(1:numel(Ca)));
    off = off + numel(Ca);
    xSys2(Tc) = xSys2tmpOLD(off+(1:numel(Tc)));
    off = off + numel(Tc);
    xSys2(Ta) = xSys2tmpOLD(off+(1:numel(Ta)));

    xSys1(cotree_dofs) = matSys1 \ rhs(cotree_dofs);
    xSys1(tree_dofs_air) = - B_contr_air * xSys1(cotree_dofs_air);
    xSys1(tree_dofs_cond) = - B_contr_cond * xSys1(cotree_dofs_cond);

    %% compute residuals
    K = stiff_mat + mass_mat;
    K = K(int_dofs, int_dofs);
    rhs = rhs(int_dofs);
    nb = norm(rhs);
    resOrig(iter) = norm(K * xOrig(int_dofs) - rhs) / nb;
    res1(iter) = norm(K * xSys1(int_dofs) - rhs) / nb;
    res2(iter) = norm(K * xSys2(int_dofs) - rhs) / nb;
    res2OLD(iter) = norm(K * xSys2OLD(int_dofs) - rhs) / nb;

    S = calculateGaugingMatrix(space_eqs, space, msh, sigma_coeffs + 1j * omega * epsilon_coeffs);
    S = S(int_dofs_eqs, int_dofs);
    
    % M1 = op_gradu_gradv_mp (space_eqs, space_eqs, msh_eqs, @(x,y,z) ones(size(x)));
    % M1 = M1(int_dofs_eqs, int_dofs_eqs);

    divorigvec = S*xOrig(int_dofs);
    divOrig(iter) = norm(divorigvec);
    div1vec = S*xSys1(int_dofs);
    div1(iter) = norm(div1vec);
    div2vec = S*xSys2(int_dofs);
    div2(iter) = norm(div2vec);

    %ranks(1, iter) = rank(full(matOrig), 1e2);
    %ranks(2, iter) = rank(full(matSys1), 1e2);
    %ranks(3, iter) = rank(full(matSys2), 1e2);

end

filename = "vts/planarCoilExample.csv";
fprintf("Saved results in " + filename + "\n")
data = [freqs', condOriginal', condMethod1', condMethod2', resOrig', res1', res2', divOrig', div1', div2'];
writematrix(data, filename)

dat = readmatrix(filename);
loglog(dat(:,1), dat(:,5), 'LineWidth',2)
hold on
loglog(dat(:,1), dat(:,6), 'LineWidth',2)
loglog(dat(:,1), dat(:,7),'--', 'LineWidth',2)
legend("orig", "1", "2")
xlabel("f in hz")
title("||K*x-b||/||b||")

figure
loglog(dat(:,1), dat(:,2), 'LineWidth',2)
hold on
loglog(dat(:,1), dat(:,4), 'LineWidth',2)
loglog(dat(:,1), dat(:,3),'--', 'LineWidth',2)
legend("orig", "1", "2")
xlabel("f in hz")
title("\kappa(K)")

figure
loglog(dat(:,1), dat(:,8), 'LineWidth',2)
hold on
loglog(dat(:,1), dat(:,10), 'LineWidth',2)
loglog(dat(:,1), dat(:,9),'--', 'LineWidth',2)
legend("orig", "1", "2")
xlabel("f in hz")
title("||(S_\sigma+jwS_\epsilon)x||_{L2}")

% figure
% loglog(dat(:,1), condEIGSOrig, '-*', 'LineWidth',2)
% hold on
% loglog(dat(:,1), condEIGS1, '-*', 'LineWidth',2)
% loglog(dat(:,1), condEIGS2, '-*', 'LineWidth',2)
% legend("orig", "1", "2")
% xlabel("f in hz")
% title("cond (eigs)")

% figure
% scatter(freqs,abs((max(eigValsOrig))))
% hold on
% scatter(freqs,abs((min(eigValsOrig))))
% scatter(freqs,abs((max(eigVals2))), '*')
% hold on
% scatter(freqs,abs((min(eigVals2))), '*')
% title("|(lambda)|")
% legend("max orig", "min orig", "max 2", "min 2")
% xlabel("f in hz")
% ylabel("eigenvalue")
% set(gca,'xscale','log')
% set(gca,'yscale','log')

function vals = eqsDirichletFun(x,y,z,ind,patchesHighBnd,patchesLowBnd)
if ismember(ind,patchesHighBnd)
    vals = 1*ones(size(x));
    return
end
vals = zeros(size(x));
end

