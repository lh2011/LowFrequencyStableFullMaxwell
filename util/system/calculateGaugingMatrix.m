function B = calculateGaugingMatrix(space_eqs, space, msh, coeffs)

B = sparse(space_eqs.ndof, space.ndof);
for i=1:numel(coeffs)
    B = B + coeffs(i) * op_v_gradp_mp (space, space_eqs, msh, @(x,y,z) ones(size(x)), i);
end

end