function paras_sg = stochcol_sg(X, rule_id)
%STOCHCOL_SG generate attributes of the general sparse grid interpolation
% given an index set and a sequence of 1D collocation nodes
%
% paras_sg = stochcol_sg(X, rule_id)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

% indset is a subset of X, accounts for nonzero c_i
% sum of tensor products
[cset, indset, mindset] = stochcol_gsgterms(X, rule_id); 
gridd = stochcol_getgrid(mindset);
% sum of single-term polynomials
[gridd, clincombiset, indlincombiset, mindlincombiset, multiweights]...
    = stochcol_multilag(gridd, cset, indset, mindset, rule_id); 
% coordinates of underlying grid
coords = stochcol_getcoords(gridd,rule_id);
% outputs
paras_sg = {cset, indset, mindset, gridd, clincombiset, ...
            indlincombiset, mindlincombiset, multiweights, coords};