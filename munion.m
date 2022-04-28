function C = munion(A)
%MUNION Multiple sets union
% C = mintersect(A)
% input
%     A  cell vector
% output
%     C  vector

C = A{1};
if length(A) > 1
    for k = 2:length(A)
        C = union(C,A{k});
    end
end
