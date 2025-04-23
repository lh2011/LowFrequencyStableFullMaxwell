function [M, r] = create_system_matrices_EQS(space, msh, epsilon_coeffs, sigma_coeffs, f, omega)

% mass mat contribution for every patch
M = sparse(space.ndof, space.ndof);
r = sparse(space.ndof, 1);
for i=1:numel(sigma_coeffs)
    M = M + (sigma_coeffs(i) + 1j * omega * epsilon_coeffs(i)) * op_gradu_gradv_mp (space, space, msh, @(x,y,z) ones(size(x)), i);
    if nargin(f)==3
        f_patch = @(x,y,z) f(x,y,z);
    else
        f_patch = @(x,y,z) f(x,y,z,i);
    end
    r = r + op_f_v_mp (space, msh, f_patch); 
end

   
end
