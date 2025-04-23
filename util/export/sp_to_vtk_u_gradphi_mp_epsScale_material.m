% SP_TO_VTK: Export multipatch results to VTK format for plotting.


function sp_to_vtk_u_gradphi_mp_epsScale_material (u, space_u, phi, space_phi, geometry, npts, filename, scaling_coeffs, epsilon_coeffs, sigma_coeffs, nu_coeffs, fieldname, varargin)

  str1 = cat (2,'<?xml version="1.0"?> \n', ...
'<VTKFile type="Collection" version="0.1"> \n', ...
'<Collection> \n');

  str2 = cat (2, '<DataSet part="%d" file="%s.vts"/> \n');

  str3 = cat (2, ...
'</Collection>\n', ...
'</VTKFile> \n');

  if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
    pvd_filename = cat (2, filename, '.pvd');
  else
    pvd_filename = filename;
    filename = filename (1:end-4);
  end

  fid = fopen (pvd_filename, 'w');
  if (fid < 0)
    error ('mp_sp_to_vtk: could not open file %s', pvd_filename);
  end

  fprintf (fid, str1);
  ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
  if (isempty (ind)); ind = 1; end
  for iptc = 1:space_u.npatch
    filename_patch_without_path = cat (2, filename(ind:end), '_', num2str (iptc));
    filename_patch = cat (2, filename, '_', num2str (iptc));
    fprintf (fid, str2, iptc, filename_patch_without_path);
    if (isempty (space_u.dofs_ornt))
      sp_to_vtk_u_gradphi (scaling_coeffs(iptc) * u(space_u.gnum{iptc}), space_u.sp_patch{iptc}, scaling_coeffs(iptc) * phi(space_phi.gnum{iptc}), space_phi.sp_patch{iptc}, geometry(iptc), npts, ...
                           filename_patch, epsilon_coeffs(iptc), sigma_coeffs(iptc), nu_coeffs(iptc), fieldname, varargin{:})
    else
      sp_to_vtk_u_gradphi_material (scaling_coeffs(iptc) * u(space_u.gnum{iptc}) .* space_u.dofs_ornt{iptc}', space_u.sp_patch{iptc}, scaling_coeffs(iptc) * phi(space_phi.gnum{iptc}), space_phi.sp_patch{iptc}, geometry(iptc), npts, ...
                           filename_patch, epsilon_coeffs(iptc), sigma_coeffs(iptc), nu_coeffs(iptc), fieldname, varargin{:})
    end
  end
  fprintf (fid, str3);

  fclose (fid);

end
