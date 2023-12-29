clear all, close all, clc
%load hald;  % Load Portlant Cement dataset
% a
A = rand(8000,10000);
% b
b = rand(1,10000)';
[U,S,V] = svd(A,'econ');
% c
x = V*inv(S)*U'.*b;                    % Solve Ax=b using the SVD

% plot(b,'k','LineWidth',2);  hold on            % Plot data
% plot(A*x,'r-o','LineWidth',1.,'MarkerSize',2); % Plot regression
% l1 = legend('Heat data','Regression')
%x = regress(b,A); %% Alternative 1  (regress)

m = pinv(A); %% Alternative 2  (pinv)
x_appro = m.*b;
% d
% Print first 5 elts of x
first_5_elts_of_x = x(1:5)
% Print entries (1,1) and (1,2) of the pseudoinverse of A
entries_11 = m(1,1)
entries_12 = m(1,2)
% e
% Compute the residual error
residual_error = norm(x - x_appro)
