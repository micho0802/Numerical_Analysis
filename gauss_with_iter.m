A = [2 -1 1;2 2 2;-1 -1 2];
b = [-1 4 5]';
x0 = [0 0 0]';
%Part c
D = diag(diag(A));
U = triu(A,1);
L = tril(A,-1);
T = inv(L+D)*U;
spec_radius = max(abs(eig(T)));
disp('Eigenvalue for matrix A')
spec_radius
%Part d
iter = 25;
N = length(b);
x = zeros(N,1);
y = zeros(N,1);
tol=1e-5; %10^-5 

for j=1:iter
    for i=1:N
        x(i) = (b(i)/A(i,i)) - (A(i,[1:i-1,i+1:N])*x0([1:i-1,i+1:N]))/A(i,i);
        x0(i) = x(i);
    end
fprintf ('Interation # %d\n', j)
x
if abs(y-x)< tol 
    break
end    
y=x
end
