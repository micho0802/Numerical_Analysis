clc
clear
rng('default') % For reproducibility
x = rand(100,1);
y = 20*x + randn(100,1)/5;
lambda = 1e-03;
B = lasso(x,y,'Lambda',lambda)
%B = 1.9825
scatter(x,y)
hold on
xd = 0:0.1:1;
plot(xd,xd*B)
grid on
hold off