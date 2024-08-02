function B = avg(A,k) % Average of the adjacent elt in A over k (iters)
if nargin<2, k = 1; end % If #of input arguments < 2 then k = 1 (default)
if size(A,1)==1, A = A'; end % If A is a row vector then A' to make a column vector. Hence average along row
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end % If k > 2, the function call itself recursively with 'k-1' until k<2
if size(A,2)==1, B = B'; end % If A is a column vector then transpose B back to a row vector