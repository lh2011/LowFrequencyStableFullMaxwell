% MSH_TO_VTK: Export to VTK format for plotting.

function msh_to_vtk_material (pts, values, epsilon_coeff, sigma_coeff, nu_coeff, filename, fieldnames)

  if (iscell (values))
    if (numel (values) ~= numel (fieldnames))
      error ('The number of fields and the number of names should be the same')
    end
  else
    values = {values};
    if (~iscell (fieldnames))
      fieldnames = {fieldnames};
    end
  end

  ndim = numel (size (pts)) - 1;
  rdim = size (pts, 1);
  
  str1 = cat (2,'<?xml version="1.0"?> \n', ...
'<VTKFile type="StructuredGrid" version="0.1"> \n', ...
'<StructuredGrid WholeExtent="0 %d 0 %d 0 %d"> \n', ...
'<Piece Extent="0 %d 0 %d 0 %d"> \n', ...
'<PointData>\n');

  str1b = cat (2, ...
'<DataArray type="Float32" Name="%s" format="ascii" NumberOfComponents="%d"> \n');

  str1c = cat (2,'\n</DataArray>');

  str1d = cat(2, '\n<DataArray type="Float32" Name="Eps" format="ascii" NumberOfComponents="3"> \n');

  str1e = cat(2, '\n<DataArray type="Float32" Name="Sigma" format="ascii" NumberOfComponents="3"> \n');

  str1f = cat(2, '\n<DataArray type="Float32" Name="Nu" format="ascii" NumberOfComponents="3"> \n');


  str2 = cat (2,'\n</PointData> \n', ...
'<Points> \n', ...
'<DataArray type="Float32" NumberOfComponents="3"> \n');

  str3 = cat (2, '\n', ...
'</DataArray>\n', ...
'</Points> \n', ...
'</Piece> \n', ...
'</StructuredGrid> \n', ...
'</VTKFile> \n');

% Even for 2D data, everything is saved in 3D 
% ndims (or size) do not work properly for the 1D case. I remove singleton
% dimensions using this trick

  size_pts = size (pts);
  npts = size_pts (2:end);
  
  if (ndim < 3)
    npts (ndim+1:3) = 1;
  end
  if (rdim < 3)
    pts(ndim+1:3,:,:) = 0;
  end

  if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.vts'))
    filename = cat (2, filename, '.vts');
  end

  fid = fopen (filename, 'w');
  if (fid < 0)
    error ('msh_to_vtk: could not open file %s', filename);
  end

  fprintf (fid, str1, ...
           npts(1)-1, npts(2)-1, npts(3)-1, ...
           npts(1)-1, npts(2)-1, npts(3)-1);

  for iopt = 1:numel(values)
    if (sum (size (values{iopt}) > 1) == ndim)
      ncomp = 1;
    elseif (sum (size (values{iopt}) > 1) == ndim + 1)
      ncomp = 3;
    else
      ncomp = 9;
    end
    if (ncomp == 3 && rdim < 3)
      values{iopt}(rdim+1:3,:,:) = 0;
    elseif (ncomp == 9 && rdim < 3)
      values{iopt}(rdim+1:3,rdim+1:3,:,:) = 0;
    end
    fprintf (fid, str1b, fieldnames{iopt}, ncomp);
    fprintf (fid, '%g ', values{iopt}(:));
    fprintf (fid, str1c);
  end

  fprintf (fid, str1d);
  valuesEps = 1/sqrt(3) * epsilon_coeff * ones(size(values{iopt}(:)));
  fprintf (fid, '%g ', valuesEps);
  fprintf (fid, str1c);

  fprintf (fid, str1e);
  valuesSigma = 1/sqrt(3) * sigma_coeff * ones(size(values{iopt}(:)));
  fprintf (fid, '%g ', valuesSigma);
  fprintf (fid, str1c);

  fprintf (fid, str1f);
  valuesNu = 1/sqrt(3) * nu_coeff * ones(size(values{iopt}(:)));
  fprintf (fid, '%g ', valuesNu);
  fprintf (fid, str1c);

  fprintf (fid, str2);

  fprintf (fid, '%g ', pts(:));
  fprintf (fid, str3);

  fclose (fid);

end
