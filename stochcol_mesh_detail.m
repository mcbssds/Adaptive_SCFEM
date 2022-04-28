function paras_detail = stochcol_mesh_detail(paras_fem)
%STOCHCOL_MESH_DETAIL detail space generator
%
%  paras_detail = stochcol_mesh_detail(paras_fem)
%
%   TIFISS function: FX 28 November 2019; DJS 24 August 2021
% Copyright (c) 2019 F. Xu

xy  = paras_fem{1};
evt = paras_fem{2};
% Computing edge lenghts/connections
%fprintf('Computing edge lengths/connections...');
%edgeGenTime = tic;
[eex, tve, els] = tedgegen(xy, evt);
%fprintf('done (%.5f sec)\n',toc(edgeGenTime));
% Computing detail space Y
%fprintf('Computing detail space Y info...');
%detailTime = tic;
[evtY, xyY, boundY, Ybasis] = p1grid_detail_space(xy, evt);
%fprintf('done      (%.5f sec)\n',toc(detailTime));
paras_detail = {eex, tve, els, evtY, xyY, boundY, Ybasis};
