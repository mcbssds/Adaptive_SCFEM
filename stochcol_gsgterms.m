function [cset, indset, mindset] = stochcol_gsgterms(I, rule_id)
%STOCHCOL_GSGTERMS nonzero terms of general sparse grid interpolant
%
% [cset, indset, mindset] = stochcol_gsgterms(X, rule_id)
%
% The interpolant is written as a sum of full tensor operators
%
% S_I[u](\bm{y}) = \sum_{\bm{j}\in J} c_{\bm{j}}\bigotimes_{n = 1}^{N}
%                  \mathcal{U}^{n, M(j_n)}(\bm{y})
%
% where M is the rule determined by rule_id, and
%  J = \{\bm{j} = \bm{i} -\bm{\alpha}: j_n >0 \text{ for }  
%      n = 1, \dots, N, \forall \bm{i} \in I, 
%      \forall \bm{\alpha} \in \{0, 1\}^N\}, 
%  c_{\bm{j}} = \sum_{\substack{\bm{\alpha} \in \{0, 1\}^N 
%               (\bm{\alpha} + \bm{j}) \in I}} (-1)^{|\bm{\alpha}|}       
%      
% input:
%       X    multi-index set
%      
% output:
%      cset    non-zero c_{\bm{j}}
%    indset    subset of J corresponding cset
%   mindset    M(indset)
%
%   TIFISS function: FX 16 April 2019
% Copyright (c) 2019 F. Xu

N = size(I, 2);  % dimension (i.e., number of parameters)
% generates the set {0,1}^N
alphaset = stochcol_getalphaset(N);
[ar,~] = size(alphaset);
% Generate the set J
J = [];
[ir,~] = size(I);
for i = 1:ir
    for j = 1:ar
        temp = I(i,:) - alphaset(j,:);
        if min(temp)>0
            J = [J; temp];
        end
    end
end
J = unique(J,'rows');
[Jr,~] = size(J);
% Generate the set c
c = zeros(Jr,1);
for j = 1:Jr
    for i = 1:ar
        if ismember(alphaset(i,:) + J(j,:), I, 'rows')
            c(j) = c(j) + (-1)^norm(alphaset(i,:),1);
        end
    end
end
% determine the rule
if rule_id == 1 
    rule = @(i) i; % linear rule
elseif rule_id == 2
    rule = @(i) (i==1).*1 + (i>1).*(2.^(i-1)+1); % doubling rule
end
cset = [];
indset = [];
mindset = [];
for j = 1:Jr
    if c(j) ~= 0
        cset = [cset; c(j)];
        indset = [indset; J(j,:)];
        mindset = [mindset; rule(J(j,:))];
    end
end
end
