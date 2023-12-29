A = [2 -1 1;2 2 2;-1 -1 2];
b = [-1 4 5]';
x0 = [0 0 0]';
%Part a
D = diag(diag(A));
R = A - D;
T = - inv(D) * R
spec_radius = max(abs(eig(T)));
disp('Eigenvalue for matrix A')
spec_radius
%Part b
iter = 25;
N = length(b);
x = zeros(N,1);

for j=1:iter
    for i=1:N
        x(i) = (b(i)/A(i,i)) - (A(i,[1:i-1,i+1:N])*x0([1:i-1,i+1:N]))/A(i,i);
    end
fprintf ('Interation # %d\n', j)
x
x0 = x;
end


