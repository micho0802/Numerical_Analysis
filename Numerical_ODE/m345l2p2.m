%Problem 2
y = input('y='); %input for IC of y
dy = input('dy='); %input for IC of y

a = input('a='); %input for a
b = input('b='); %input for b
c = input('c='); %input for c
tmin = input('tmin='); %input for the min
tmax = input('tmax='); %input for the max

%Part 1
m = [a b c]; 
r = roots(m); %Roots of the char eq

%Part 2
%print the general sol,
fprintf('y(t)= C1*exp^(%.2f*t) + C2*%.2f*exp^(%.2f*t)\n', r(1), r(1), r(2)); 

%Find C1 and C2 and substitute in the DE
A = [1 0; r(1) 1];
B = [y dy];
x = A/B;
fun = sprintf('y(t) = %.2f*exp^(%.2f*t) + %.2f*%.2f*exp^(%.2f*t)\n', x(1), r(1), x(2), r(2), r(2));

%Part4
%Display the function
disp(fun);

%Part5
f = @(t) (x(1) + x(2).*t).*exp(r(2)*t) ;
fplot(f,[tmin, tmax]);
grid on;
xlabel('t');
ylabel('y(t)');
title(fun);