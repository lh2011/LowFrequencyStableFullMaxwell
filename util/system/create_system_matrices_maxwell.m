function [M, K, r] = create_system_matrices_maxwell(space, msh, epsilon_coeffs, nu_coeffs, sigma_coeffs, f, omega)

% mass mat contribution for every patch
M = sparse(space.ndof, space.ndof);
K = sparse(space.ndof, space.ndof);
r = sparse(space.ndof, 1);
for i=1:numel(sigma_coeffs)
    M = M + (1j * omega * sigma_coeffs(i) - omega^2 * epsilon_coeffs(i)) * op_u_v_mp (space, space, msh, @(x,y,z) ones(size(x)), i);
    K = K + nu_coeffs(i) * op_curlu_curlv_mp (space, space, msh, @(x,y,z) ones(size(x)), i);
    if nargin(f)==3
        f_patch = @(x,y,z) f(x,y,z);
    else
        f_patch = @(x,y,z) f(x,y,z,i);
    end
    r = r + op_f_v_mp (space, msh, f_patch, i);
end

end
