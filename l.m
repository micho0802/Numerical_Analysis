%matrix A and vector b
function [x1, x2] = l(A,b)
x1 = A\b;
x2 = inv(A)*b;
end