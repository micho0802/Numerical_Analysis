function [a] = hiber(m,n);
for i=m:n
    for j=m:n
    a(i,j) = 1/(i+j-1);
    end
end
%I'm not sure how to compute the inverse of matrix A to become Hilbert
%Matrix. Because the result I computed is already a Hilbert matrix.