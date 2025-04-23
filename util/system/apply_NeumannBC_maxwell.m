function rhs = apply_NeumannBC_maxwell(rhs, boundaries, g, space, msh, nmn_sides)

Nbnd = cumsum ([0, boundaries.nsides]);
for iref = nmn_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    gref = @(varargin) g(varargin{:},iref);
    rhs_nmnn = op_f_v_mp (space.boundary, msh.boundary, gref, iref_patch_list);
    rhs(space.boundary.dofs) = rhs(space.boundary.dofs) + rhs_nmnn .* space.boundary.boundary_orientation.';
end
end