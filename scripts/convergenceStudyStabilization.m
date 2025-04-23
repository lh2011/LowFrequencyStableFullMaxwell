%% Problem specification
clear all;
close all;

eps0 = 8.854187e-12;
sig0 = 6e7;
artificSig = 6e7;
nu0 = 1/(4*pi)*1e7;

srf = nrbsquare([0 0],pi,pi);
geo = nrbtform(nrbextrude(srf,[0 0 pi]),vectrans(pi/2*ones(3,1)));
%subs = 4;

exportToParaview = false;
for f = [10 1000000]
    omega = 2*pi*f;                             
    
    sigma_coeffs =      [0] * sig0; %[0 0 0 0 5 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 5 0 0 0 0] * sig0;                      % electric conductivity in every patch -> has to be suited to the geometry
    epsilon_coeffs =    [1] * eps0; %[5 5 5 5 5 5 5 5 5 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 5 5 5] * eps0;                      % electric permittivity
    nu_coeffs = nu0*ones(size(sigma_coeffs));
    
    nuFun  = @(x, y, z) nu0*ones(size(x));
    sigmaFun = @(x, y, z) sigma_coeffs(1)*ones(size(x));
    epsiFun = @(x, y, z) eps0*ones(size(x));
    kappaFun = @(x, y, z) sigmaFun(x,y,z) + 1i*omega*epsiFun(x,y,z);
    
    nu_mat = @(x, y, z) cat(1,reshape(nuFun(x,y,z),[1,size(x)]),...
                              reshape(nuFun(x,y,z),[1,size(x)]),...
                              reshape(nuFun(x,y,z),[1,size(x)]));
    kappa_mat = @(x, y, z) cat(1,reshape(kappaFun(x,y,z),[1,size(x)]),...
                                 reshape(kappaFun(x,y,z),[1,size(x)]),...
                                 reshape(kappaFun(x,y,z),[1,size(x)]));
    
    Afun = @(x, y, z) cat(1, ...
        reshape(    sin(x).*cos(y).*cos(z),[1,size(x)]),...
        reshape(-2* cos(x).*sin(y).*cos(z),[1,size(x)]),...
        reshape(    cos(x).*cos(y).*sin(z),[1,size(x)]));
    
    Bfun = @(x, y, z) cat(1, ...
        reshape(-3*cos(x).*sin(y).*sin(z),[1,size(x)]),...
        zeros([1,size(x)]),...
        reshape(3*sin(x).*sin(y).*cos(z),[1,size(x)]));
    
    Ufun = @(x, y, z) cos(x).*cos(y).*cos(z);
    
    % negative E-field - just grad(u)
    nEfun = @(x, y, z) -cat(1,reshape(sin(x).*cos(y).*cos(z),[1,size(x)]),...
                              reshape(cos(x).*sin(y).*cos(z),[1,size(x)]),...
                              reshape(cos(x).*cos(y).*sin(z),[1,size(x)]));
    
    
    %% Construct right hand sides
    ccAfun = @(x, y, z) cat(1,reshape(3 *   sin(x).*cos(y).*cos(z),[1,size(x)]),...
                              reshape(-6*   cos(x).*sin(y).*cos(z),[1,size(x)]),...
                              reshape(3 *   cos(x).*cos(y).*sin(z),[1,size(x)]));
    
    rhsVecFun = @(x, y, z) nu_mat(x,y,z).*ccAfun(x,y,z) + 1i*omega*kappa_mat(x,y,z).*Afun(x,y,z);
    rhsScaFun = @(x, y, z) 3*kappaFun(x,y,z).*Ufun(x,y,z);
    
    % rhsVec = Js - kappa*grad(u) <=> Js = rhsVec + kappa*grad(u)
    Js = @(x, y, z) rhsVecFun(x, y, z) + kappa_mat(x,y,z) .* nEfun(x,y,z);
    
    %rho = @(x,y,z) 3*sin(x).*sin(y).*sin(z)*(sig0 + eps0*omega*1i);
    
    % Neumann functions
    g_maxwell = @(x, y, z, ind) cat(1, ...
        zeros ([1, size(x)]), ...
        zeros ([1, size(x)]), ...
        zeros ([1, size(x)]));
    
    g_eqs = @(x,y,z,ind) zeros(size(x));
    
    % Dirichlet functions
    h_maxwell = @(x, y, z, ind) Afun(x,y,z);
    
    h_eqs = @(x,y,z,ind) Ufun(x,y,z);
    
    % boundary sides
    drchlt_sides_eqs    = 1:6; %[11 14 17 20]; %[3 5 8 10 11 13 16 18 21 36 38 41 43 44 46 49 51 54];
    nmn_sides_eqs       = setdiff(1:6, drchlt_sides_eqs);
    
    nmn_sides_maxwell_global    = [];
    drchlt_sides_maxwell_global  = setdiff(1:6, nmn_sides_maxwell_global);
    
    subsArr = [2 4 8 16];
    degrees = [1 2 3];
    
    H1 = zeros(numel(degrees),numel(subsArr));
    Hcurl = zeros(numel(degrees),numel(subsArr));
    L2 = zeros(numel(degrees),numel(subsArr));
    Hcurlsemi = zeros(numel(degrees),numel(subsArr));
    
    for d = 1:numel(degrees)
        for s = 1:numel(subsArr)
        
            deg = degrees(d);
            subs = subsArr(s);
    
            %% EQS: create geometry, space and mesh
            [space_eqs, msh_eqs, geometry, boundaries] = create_space_mesh_EQS(geo, deg, subs);
            
            %% EQS: assemble matrices
            [eqs_mat, eqs_rhs] = create_system_matrices_EQS(space_eqs, msh_eqs, epsilon_coeffs, sigma_coeffs, rhsScaFun, omega);
            
            disp("Created EQS System")
            
            %% EQS: Apply Neumann boundary conditions
            %eqs_rhs = apply_NeumannBC_eqs(eqs_rhs, boundaries, g_eqs, space_eqs, msh_eqs, nmn_sides_eqs);
            
            %% EQS: Apply Dirichlet boundary conditions (L2-projection)
            phi = zeros (space_eqs.ndof, 1);
            [phi_drchlt, drchlt_dofs_eqs] = sp_drchlt_l2_proj (space_eqs, msh_eqs, h_eqs, drchlt_sides_eqs);
            %phi(drchlt_dofs_eqs) = phi_drchlt;
            int_dofs_eqs = setdiff (1:space_eqs.ndof, drchlt_dofs_eqs);
            
            %% EQS: Solve the linear system
            %eqs_rhs(int_dofs_eqs) = eqs_rhs(int_dofs_eqs) - eqs_mat(int_dofs_eqs, drchlt_dofs_eqs)*phi_drchlt;
            phi(int_dofs_eqs) = eqs_mat(int_dofs_eqs, int_dofs_eqs) \ eqs_rhs(int_dofs_eqs);
            
            disp('Solved EQS')
            
            %% MAXWELL: create geometry, space and mesh
            [space, msh, ~, boundaries] = create_space_mesh_maxwell(geo, deg, subs);
            
            % boundary sides
            nmn_sides_maxwell = [];
            for ni = nmn_sides_maxwell_global
                nmn_sides_maxwell   = union(nmn_sides_maxwell, find ([boundaries.faces] == ni));
            end
            drchlt_sides_maxwell  = setdiff(1:numel(boundaries), nmn_sides_maxwell);
            
            %% MAXWELL: create system matrices
            %[mass_mat, stiff_mat, f_maxwell] = create_system_matrices_maxwell(space, msh, epsilon_coeffs, nu0*ones(size(sigma_coeffs)), sigma_coeffs, f, omega, @(x,y,z) ones(size(x)));
            [mass_mat, stiff_mat, f_maxwell] = create_system_matrices_maxwell(space, msh, epsilon_coeffs, nu0, sigma_coeffs, Js, omega);
            % EQS modification
            G = calculateGaugingMatrix(space_eqs, space, msh, sigma_coeffs + 1i * omega * epsilon_coeffs);
            maxwell_rhs = f_maxwell - G' * phi;
            
            disp("Created Maxwell System")
            
            %% MAXWELL: Apply Neumann boundary conditions
            %maxwell_rhs = apply_NeumannBC_maxwell(maxwell_rhs, boundaries, g_maxwell, space, msh, nmn_sides_maxwell);
            
            %% MAXWELL: Apply Dirichlet boundary conditions (L2-projection)
            u = zeros (space.ndof, 1);
            [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h_maxwell, drchlt_sides_maxwell);
            %u(drchlt_dofs) = u_drchlt;
            
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
            
            disp("calculated dofs for different regions")
            
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
            
            % tree_dofs_dirAndint = union(intersect(tree_dofs, interface_idx), intersect(tree_dofs, drchlt_dofs));
            % nodes_dirAndint = union(drchlt_dofs_eqs_maxwellBnd, interface_nodes);
            
            % add interface dofs into cond dofs
            cotree_dofs_cond = union(cotree_dofs_cond, cotree_dofs_int);
            tree_dofs_cond = union(tree_dofs_cond, tree_dofs_int);
            conductor_nodes = union(conductor_nodes, interface_nodes);
             
            disp('Constructed and split tree')
            
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
                    gamma*B_air(air_nodes, tree_dofs_cond),     gamma*B_air_t                               ];
            
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
                   
            %% --
            rhs = [ maxwell_rhs(Cc) - K(Cc,Cdc) * u(Cdc) - K(Cc,Cda) * u(Cda) - K(Cc, Tdc) * u(Tdc) - K(Cc, Tda) * u(Tda);
                    maxwell_rhs(Ca) - K(Ca,Cdc) * u(Cdc) - K(Ca,Cda) * u(Cda) - K(Ca, Tdc) * u(Tdc) - K(Ca, Tda) * u(Tda);
                    zeros(numel(T),1)];

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
            
            disp('Solved Maxwell step')
            
            %% post
            [H1(d,s)] = sp_h1_error(space_eqs, msh_eqs, phi, Ufun, nEfun);         
            [Hcurl(d,s),L2(d,s),Hcurlsemi(d,s)] = sp_hcurl_error(space, msh, u, Afun, Bfun);

            % if exportToParaview
            % output_file = 'vts/TESTCUBE';
            % vtk_pts = {linspace(0, 1, 15), linspace(0, 1, 15), linspace(0, 1, 15)};
            % scaling_coeffs = sigma_coeffs; %ones(size(epsilon_coeffs));
            % %sp_to_vtk(phi, space_eqs, geometry, vtk_pts, output_file,{'U', 'E'},{'value', 'gradient'})
            % sp_to_vtk(u, space, geometry, vtk_pts,output_file,{'A', 'B'},{'value', 'curl'})
            % %sp_to_vtk_u_gradphi_mp_epsScale_material(1j*omega*u, space, phi, space_eqs, geometry, vtk_pts, output_file, scaling_coeffs, epsilon_coeffs, sigma_coeffs, nu_coeffs, {'J'}, {'value'})
            % fprintf ('The result is saved in the file %s \n \n', output_file);
            % end
        
        end
    end
    
    %%
    fprintf('Convergence H1 of phi:\n')
    H1 = abs(H1);
    rates_H1 = diff(log10(H1),1,2)./repmat(diff(log10(subsArr)),size(H1,1),1);
    disp(rates_H1);
    
    fprintf('Convergence Hcurl of A:\n')
    Hcurl = abs(Hcurl);
    Hcurlsemi = abs(Hcurlsemi);
    rates_Hcurl = diff(log10(abs(Hcurl)),1,2)./repmat(diff(log10(subsArr)),size(Hcurl,1),1);
    disp(rates_Hcurl);
    
    
    if 1
        figure(1)
        clf()
        loglog(subsArr, Hcurl(1,:),'Marker','+')
        hold on
        % loglog(subsArr, Hcurl(2,:),'Marker','+')
        % loglog(subsArr, Hcurl(3,:),'Marker','+')
        loglog(subsArr, Hcurlsemi(1,:),'Marker','o')
        % loglog(subsArr, Hcurlsemi(2,:),'Marker','o')
        % loglog(subsArr, Hcurlsemi(3,:),'Marker','o')
        ylabel("value of error")
        xlabel("nsubs")
        legend("p=1", "p=2", "p=3", 'Location','southwest');
        grid on
    end
    
    writematrix([subsArr',H1',Hcurl',Hcurlsemi'],strcat('errors_stabilization_f',num2str(f),'_mat.csv'))
end