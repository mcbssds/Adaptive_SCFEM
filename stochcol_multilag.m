function [grid,clincombiset, indlincombiset, mindlincombiset, ...
    multiweights]= stochcol_multilag(grid, cset, indset, mindset, rule_id)
%STOCHCOL_MULTILAG generate sparse grid interpolant
% expressed in terms of the sum of multivariate Lagrange polynomials
% 
% [grid,clincombiset, indlincombiset, mindlincombiset, ...
%    multiweights]= stochcol_multilag(grid, cset, indset, mindset, rule_id)
%
% The interpolant is written as follows
%
%  S_I[u](\bm{y}) = \sum_{\bm{m}} u(\bm{y}_{\bm{m}}) L_{\bm{m}}(\bm{y})
%
% where the multivariate Lagrange polynomials are written as
%
%  L_{\bm{m}}(\bm{y}) = \sum_{\substack{M(\bm{i})\geq \bm{m} \\ \bm{i} 
%                       \in I}}c_{\bm{i}} L_{\bm{m}}^{M(\bm{i})}(\bm{y})
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

[gr, gc] = size(grid);
[ir, ~] = size(indset);
clincombiset = cell(gr,1);
indlincombiset = cell(gr,1);
mindlincombiset = cell(gr,1);
multiweights = zeros(gr,1);
for i=1:gr
    pts = grid(i,:);
    for j = 1:ir
        c  = cset(j);
        index  = indset(j,:);
        mindex = mindset(j,:);
        if min(mindex-pts)>=0
            clincombiset{i} = [clincombiset{i}; c];
            indlincombiset{i} = [indlincombiset{i}; index];
            mindlincombiset{i} = [mindlincombiset{i}; mindex];
            temp = 1;
            for k=1:gc
                if rule_id == 1
                    nodes = stochcol_nodes_leja(index(k));
                elseif rule_id == 2
                    nodes = stochcol_nodes_cc(index(k));
                end
                [~, weight]= stochcol_1Dlagpoly(nodes{1}(pts(k)),nodes{1});
                temp = temp*weight;
            end
            multiweights(i) = multiweights(i) + c*temp;
        end
    end
end
