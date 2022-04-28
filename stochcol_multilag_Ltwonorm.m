function L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list)
%STOCHCOL_MULTILAG_LTWONORM calculate 2-norm of multivariate Lagrange polynomials
%
%  L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

gridd           = paras_sg{4};
clincombiset    = paras_sg{5};
indlincombiset = paras_sg{6};
G = stochcol_gmatrices(gridd, clincombiset, indlincombiset, list);
L_two_norm = sqrt(diag(G));
