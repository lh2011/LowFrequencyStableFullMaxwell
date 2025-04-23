function [tree_dofs, cotree_dofs] = cotree_multipatch (sp_h1, sp_hcurl, boundaries, drchlt_sides, intrfc_dofs, geometry, isslave)
% Parameter isslave
% 0 dirichlet -> interior                           (no interface)
% 1 interface -> dirichlet -> interior              (stab. in insulator)
% 2 dirichlet -> interface -> interior              (stab. in whole domain)
% 3 dirichlet -> neumann -> interface -> interior
% 4 dirichlet -> neumann, interface -> interior
% 5 dirichlet -> interface -> neumann -> interior   

if (nargin < 7)
    isslave = 0;
end

if (sp_h1.npatch ~= sp_hcurl.npatch)
  error ('Spaces not defined on the same geometry')  
end

nmn_sides = setdiff(1:numel(boundaries), drchlt_sides);

S = zeros (sp_hcurl.ndof, 1);
T = zeros (sp_hcurl.ndof, 1);
for iptc = 1:sp_hcurl.npatch
  sph1_loc = sp_h1.sp_patch{iptc};
  sp_loc = sp_hcurl.sp_patch{iptc};
  
  S_loc = zeros(sp_loc.ndof, 1);
  T_loc = zeros(sp_loc.ndof, 1);
  
  ndof_dir_h1 = sph1_loc.ndof_dir;
  cumsum_ndof = sp_loc.cumsum_ndof;
  for icomp = 1:sp_loc.ncomp_param
    inds = cell (sp_loc.ncomp_param, 1);
    ndof = sp_loc.scalar_spaces{icomp}.ndof;
    ndof_dir = sp_loc.scalar_spaces{icomp}.ndof_dir;
    [inds{:}] = ind2sub (ndof_dir, 1:ndof);
    S_loc(cumsum_ndof(icomp)+(1:ndof)) = sub2ind(ndof_dir_h1, inds{:});
    inds{icomp} = inds{icomp} + 1;
    T_loc(cumsum_ndof(icomp)+(1:ndof)) = sub2ind(ndof_dir_h1, inds{:});
  end
  
  % Without orientation
  S(sp_hcurl.gnum{iptc}) = sp_h1.gnum{iptc}(S_loc);
  T(sp_hcurl.gnum{iptc}) = sp_h1.gnum{iptc}(T_loc);
end


Nbnd = cumsum ([0, boundaries.nsides]);
% Dirichlet boundary dofs
bnd_dofs = [];
for iref = drchlt_sides
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    
  boundary_gnum = sp_hcurl.boundary.gnum;
  bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}]);
end
drchlt_dofs = sp_hcurl.boundary.dofs(bnd_dofs);

% Neumann boundary dofs
bnd_dofs = [];
for iref = nmn_sides
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    
  boundary_gnum = sp_hcurl.boundary.gnum;
  bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}]);
end
nmn_dofs = sp_hcurl.boundary.dofs(bnd_dofs);

nmn_dofs = setdiff(nmn_dofs, drchlt_dofs);

weights = ones (size(S));
weights(drchlt_dofs) = 0;
%weights(nmn_dofs) = 0.5;
%weights(intrfc_dofs) = 2;
if isslave
    weights(intrfc_dofs) = 0;
    weights(drchlt_dofs) = 0.1;
end
if isslave == 2
    weights(intrfc_dofs) = 0.1;
    weights(drchlt_dofs) = 0;
end
if isslave == 3
    weights(intrfc_dofs) = 0.2;
    weights(nmn_dofs) = 0.1;
    weights(drchlt_dofs) = 0;
end
if isslave == 4
    weights(intrfc_dofs) = 0.1;
    weights(nmn_dofs) = 0.1;
    weights(drchlt_dofs) = 0;
end
if isslave == 5
    weights(nmn_dofs) = 0.2;
    weights(intrfc_dofs) = 0.1;
    weights(drchlt_dofs) = 0;
end
if isslave == 6
    weights(nmn_dofs) = 0.2;
    weights(drchlt_dofs) = 0.1;
    weights(intrfc_dofs) = 0;

end

G = graph(S,T,weights);
G2 = minspantree (G, 'Method', 'sparse');
% pp = plot(G,'Layout', 'force3');
% highlight(pp,G2,'EdgeColor','r','LineWidth',1.5)

g2_edges = table2array(G2.Edges);
dofs = G.findedge(g2_edges(:,1),g2_edges(:,2));


% % cut fingers (all except 1) from finger tree
% if isslave == 1
%     idx_inner_nodes = union(S(intrfc_dofs),T(intrfc_dofs));
%     idx_inner_nodes = idx_inner_nodes(2:end); % cut all except 1
%     idx_remaining_edges = setdiff(1:length(g2_edges(:,1)),union(find(ismember(g2_edges(:,1),idx_inner_nodes)),find(ismember(g2_edges(:,2),idx_inner_nodes))));
%     dofs = G.findedge(g2_edges(idx_remaining_edges,1),g2_edges(idx_remaining_edges,2));
% end

% From the numbering in the graph to the numbering in GeoPDEs
space2graph = G.findedge(S,T);
[~,graph2space] = ismember (1:sp_hcurl.ndof, space2graph);

tree_dofs = sort (graph2space(dofs));
cotree_dofs = setdiff(1:sp_hcurl.ndof, tree_dofs);

% if isslave == 1
%    tree_dofs = setdiff(tree_dofs, intrfc_dofs);
%    cotree_dofs = union(cotree_dofs, intrfc_dofs);
% end

%plot_tree (geometry, sp_h1, sp_hcurl, tree_dofs);
end


function plot_tree (geometry, sp_h1, sp_hcurl, tree_dofs)

ndim = numel(sp_h1.sp_patch{1}.degree);
for iptc = [1:sp_h1.npatch]  %1:sp_h1.npatch
  sp_patch = sp_h1.sp_patch{iptc};
  for idim = 1:ndim
    grev_pts{idim} = aveknt(sp_patch.knots{idim}, sp_patch.degree(idim)+1);
  end
  mapped_pts = geometry(iptc).map(grev_pts);
  rdim = size(mapped_pts, 1);
  
  mapped_pts = reshape (mapped_pts, [rdim, sp_patch.ndof_dir]);

  for idim = 1:ndim
  for jdim = idim+1:ndim
    vv = arrayfun(@(x) 1:x(:), sp_patch.ndof_dir, 'UniformOutput', false);
    for idof = 1:sp_patch.ndof_dir(idim)
      vv{idim} = idof;
      for jdof = 1:sp_patch.ndof_dir(jdim)
        vv{jdim} = jdof;
        xpts = mapped_pts(:,vv{:});
        plot3 (xpts(1,:),xpts(2,:),xpts(3,:), 'color', 'k');
        %plot(xpts(1,:),xpts(2,:), 'color', 'k');
        hold on
      end
    end
  end
  end

  [~,~,tree_patch] = intersect (tree_dofs, sp_hcurl.gnum{iptc});
  sp_patch = sp_hcurl.sp_patch{iptc};
  indices = cell(ndim, 1);
  for icomp = 1:sp_patch.ncomp
    dofs = sp_patch.cumsum_ndof(icomp)+1:sp_patch.cumsum_ndof(icomp+1);
    [~,~,tree_comp] = intersect(tree_patch, dofs);
    [indices{:}] = ind2sub(sp_patch.scalar_spaces{icomp}.ndof_dir,tree_comp);
% Plot edges of the tree
    for ii = 1:numel(indices{1})
      vv = cellfun(@(x) x(ii), indices, 'UniformOutput', false);
      vv{icomp} = [indices{icomp}(ii), indices{icomp}(ii)+1];
      xpts = mapped_pts(:,vv{:});
      plot3(xpts(1,:),xpts(2,:),xpts(3,:), 'color', 'r', 'LineWidth', 2);
      %plot(xpts(1,:),xpts(2,:), 'color', 'r', 'LineWidth', 2);
    end
end
  
end

end
