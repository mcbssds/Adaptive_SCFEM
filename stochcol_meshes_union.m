function [xy_union, evt_union] = stochcol_meshes_union(meshes_xy, nocp)
%STOCHCOL_MESHES_UNION generates union mesh for single-level SC

xy_union = meshes_xy{1};
for k = 1:nocp
    [xy_union] = union(xy_union, meshes_xy{k}, 'rows', 'stable');
end

DT = delaunayTriangulation(xy_union);
 %plot_mesh(DT.ConnectivityList, DT.Points, 'union mesh'); % uncomment to plot the union mesh
 %pause
evt_union = DT.ConnectivityList;
return
