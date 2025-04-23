clear all;
close all;

eps0 = 8.854187e-12;
sig0 = 6e7;
nu0 = 1/(4*pi)*1e7;

srf = nrbsquare([0 0],pi,pi);
geo = nrbtform(nrbextrude(srf,[0 0 pi]),vectrans(pi/2*ones(3,1)));
%subs = 4;

%exportToParaview = false;
output_file = 'testtesttest';

for f = [10 1000000]
    omega = 2*pi*f;                             
    
    sigma_coeffs = [0]*sig0; %[0 0 0 0 5 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 5 0 0 0 0] * sig0;     
    epsilon_coeffs = [1]*eps0; %[5 5 5 5 5 5 5 5 5 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 5 5 5] * eps0;    
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
    
    % rho = @(x,y,z) 3*sin(x).*sin(y).*sin(z)*(sig0 + eps0*omega*1i);
    
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
    drchlt_sides_eqs    = 1:6;
    nmn_sides_eqs       = setdiff( 1:6, drchlt_sides_eqs);
    
    nmn_sides_maxwell_global    = [];
    drchlt_sides_maxwell_global  = setdiff(1:6, nmn_sides_maxwell_global);
    
    degrees = [1 2 3];
    subsArr = [2 4 8 16];
    
    H1 = zeros(numel(degrees),numel(subsArr));
    Hcurl = zeros(numel(degrees),numel(subsArr));
    L2 = zeros(numel(degrees),numel(subsArr));
    Hcurlsemi = zeros(numel(degrees),numel(subsArr));
    for d = 1:numel(degrees)
        for s = 1:numel(subsArr)
        
            subs = subsArr(s);
            deg = degrees(d);
            
            fprintf('Freq.: %0.2d, Deg.: %i, Sub.: %i\n',f,deg,subs);

            %% EQS: create geometry, space and mesh
            [space_eqs, msh_eqs, geometry, boundaries] = create_space_mesh_EQS(geo, deg, subs);
            
            %% EQS: assemble matrices
            [eqs_mat, eqs_rhs] = create_system_matrices_EQS(space_eqs, msh_eqs, epsilon_coeffs, sigma_coeffs, rhsScaFun, omega);
            
            %% EQS: Apply Neumann boundary conditions
            % eqs_rhs = apply_NeumannBC_eqs(eqs_rhs, boundaries, g_eqs, space_eqs, msh_eqs, nmn_sides_eqs);
            
            %% EQS: Apply Dirichlet boundary conditions (L2-projection)
            phi = zeros (space_eqs.ndof, 1);
            [phi_drchlt, drchlt_dofs_eqs] = sp_drchlt_l2_proj (space_eqs, msh_eqs, h_eqs, drchlt_sides_eqs);
            % phi(drchlt_dofs_eqs) = phi_drchlt;
            int_dofs_eqs = setdiff (1:space_eqs.ndof, drchlt_dofs_eqs);
            
            %% EQS: Solve the linear system
            % eqs_rhs(int_dofs_eqs) = eqs_rhs(int_dofs_eqs) - eqs_mat(int_dofs_eqs, drchlt_dofs_eqs)*phi_drchlt;
            phi(int_dofs_eqs) = eqs_mat(int_dofs_eqs, int_dofs_eqs) \ eqs_rhs(int_dofs_eqs);
            
            %disp('Solved EQS step')
            
            %% MAXWELL: create geometry, space and mesh
            [space, msh, ~, boundaries] = create_space_mesh_maxwell(geo, deg, subs);
            
            % boundary sides
            nmn_sides_maxwell = [];
            for ni = nmn_sides_maxwell_global
                nmn_sides_maxwell   = union(nmn_sides_maxwell, find ([boundaries.faces] == ni));
            end
            drchlt_sides_maxwell  = setdiff(1:numel(boundaries), nmn_sides_maxwell);
            
            %% MAXWELL: create system matrices
            [mass_mat, stiff_mat, maxwell_rhs] = create_system_matrices_maxwell(space, msh, epsilon_coeffs, nu0, sigma_coeffs, Js, omega);
            % EQS modification
            G = calculateGaugingMatrix(space_eqs, space, msh, sigma_coeffs + 1i*omega*epsilon_coeffs);
            maxwell_rhs = maxwell_rhs - G' * phi;
            
            %% MAXWELL: Apply Neumann boundary conditions
            % maxwell_rhs = apply_NeumannBC_maxwell(maxwell_rhs, boundaries, g_maxwell, space, msh, nmn_sides_maxwell);
            
            %% MAXWELL: Apply Dirichlet boundary conditions (L2-projection)
            u = zeros (space.ndof, 1);
            [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h_maxwell, drchlt_sides_maxwell);
            % u(drchlt_dofs) = u_drchlt;
            int_dofs = setdiff (1:space.ndof, drchlt_dofs);
            
            %% MAXWELL: Solve the linear system
            % maxwell_rhs(int_dofs) = maxwell_rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt - mass_mat(int_dofs, drchlt_dofs)*u_drchlt;
            u(int_dofs) = (stiff_mat(int_dofs, int_dofs) +  mass_mat(int_dofs, int_dofs)) \ maxwell_rhs(int_dofs);
            
            %disp('Solved Maxwell step')
            
            %% post
            [H1(d,s)] = sp_h1_error(space_eqs, msh_eqs, phi, Ufun, nEfun);    
            [Hcurl(d,s),L2(d,s),Hcurlsemi(d,s)] = sp_hcurl_error(space, msh, u, Afun, Bfun);
            
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

    fprintf('Convergence Hcurlsemi of A:\n')
    rates_Hcurlsemi = diff(log10(abs(Hcurlsemi)),1,2)./repmat(diff(log10(subsArr)),size(Hcurl,1),1);
    disp(rates_Hcurlsemi);
    
    
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
    
    writematrix([subsArr',H1',Hcurl',Hcurlsemi'],strcat('errors_original_f',num2str(f),'_mat.csv'))
end


