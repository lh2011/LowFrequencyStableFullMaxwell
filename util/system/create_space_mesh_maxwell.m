function [space, msh, geometry, boundaries] = create_space_mesh_maxwell(geo_name, deg, nsubs)
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
    [knots, zeta] = ...
        kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
    [knots_hcurl, degree_hcurl] = knt_derham (knots, degree, 'Hcurl');

    % Construct msh structure
    rule      = msh_gauss_nodes (nquad);
    [qn, qw]  = msh_set_quad_nodes (zeta, rule);
    msh{iptc} = msh_cartesian (zeta, qn, qw, geometry(iptc));

    % Construct space structure
    scalar_spaces = cell (msh{iptc}.ndim, 1);
    for idim = 1:msh{iptc}.ndim
        scalar_spaces{idim} = sp_bspline (knots_hcurl{idim}, degree_hcurl{idim}, msh{iptc});
    end
    sp{iptc} = sp_vector (scalar_spaces, msh{iptc}, 'curl-preserving');
    clear scalar_spaces
end

%% create mp space and mesh
msh = msh_multipatch (msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
clear sp

end