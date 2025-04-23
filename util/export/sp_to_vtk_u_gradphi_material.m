% SP_TO_VTK: Export to VTK format for plotting.

function sp_to_vtk_u_gradphi_material (u, space_u, phi, space_phi, geometry, npts, filename, epsilon_coeff, sigma_coeff, nu_coeff, fieldname, varargin)

  [eu, ~] = sp_eval (u, space_u, geometry, npts, varargin{:});
  [gradphi, F] = sp_eval (phi, space_phi, geometry, npts, 'gradient');

  msh_to_vtk_material (F, eu+gradphi, epsilon_coeff, sigma_coeff, nu_coeff, filename, fieldname);

end