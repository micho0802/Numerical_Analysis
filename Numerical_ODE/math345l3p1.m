%Michael Dang - 16257750
%MATH345L
%Assignment 3

%Problem1, complex roots

syms g(t)
clear,clc

%Input

y = input('y='); %input for IC of y
dy = input('dy='); %input for IC of y

a = input('a='); %input for a
b = input('b='); %input for b
c = input('c='); %input for c
tmin = input('tmin='); %input for the min
tmax = input('tmax='); %input for the max

m = [a b c]; 
r = roots(m); %Roots of the char eq

lambda1 = real(r(1));
lambda2 = real(r(2));
mu1 = imag(r(1));
mu2 = imag(r(2));

%print the general sol, with 2 decimal places
fprintf('y(t)= exp^(%.2f*t)*{[C1*cos(%.2f*t)] + [C2*sin(%.2ft)]}\n', lambda1, mu1, -mu2); 

%Find C1 and C2 and substitute in the DE
A = [1 0; lambda1 mu1];
B = [y dy];
x = A/B;
fun = sprintf('y(t) = exp^(%.2d*t)*{[%.2f*cos(%.2f*t)] + [%.2d*sin(%.2ft)]}', lambda1, x(1), mu1, x(2), -mu2);

%Display the function
disp(fun)

%Plot the solution
f = @(t) (x(1) + x(2).*t).*exp(mu2*t) ;
fplot(f,[tmin, tmax]);
grid on;
xlabel('t');
ylabel('y(t)');
title(fun);

