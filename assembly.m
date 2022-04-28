function A = assembly(Ae,mv1,n1,mv2,n2)
%ASSEMBLY assembly global matrix/vector from element contributions
%   A = assembly(Ae,mv1,n1,mv2,n2)
%   input
%       Ae element matrix/vector
%       mv1 mv2 mapping matrix
%       (n1,n2) size of global matrix

switch nargin
    case 5
        A = sparse(n1,n2);
        [~,dim1,dim2] = size(Ae);
        for krow = 1:dim1
            nrow = mv1(:,krow);
            for kcol = 1:dim2
                ncol = mv2(:,kcol);
                A = A + sparse(nrow,ncol,Ae(:,krow,kcol),n1,n2);
            end
        end
    case 3
        A = zeros(n1,1);
        [nel,dim1] = size(Ae);
        for krow = 1:dim1
            nrow = mv1(:,krow);
            for els = 1:nel
                A(nrow(els),1) = A(nrow(els),1) + Ae(els,krow);
            end
        end
end

