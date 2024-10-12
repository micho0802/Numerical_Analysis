%Michael Dang - 16257750
%MATH345L
%Final Project
%Problem 2, Euler's method

clear;

%---------------------------------------------------------------
% Define the differential equation y'=f(t,y)
%---------------------------------------------------------------
f=inline('3-2*t-(0.5)*y','t','y');

%---------------------------------------------------------------
% Set initial condition
%---------------------------------------------------------------
t0=0; y0=1;

%---------------------------------------------------------------
% allow user to enter the stepsize h
%---------------------------------------------------------------
h=input('  Enter the stepsize h => ');

%---------------------------------------------------------------
% Here are the 4 t-values to estimate
%---------------------------------------------------------------

T=[.2 .4 .6 .8 1.0];

%---------------------------------------------------------------
% Loop to solve the IVP four times -- once for each entry in T
% Use Euler method.
%---------------------------------------------------------------
for j=1:length(T)
     
    % We need this many iterations to get to T with a stepsize of h.
    steps=round(T(j)/h);         
    
    % Start at initial condition each time through
    y_n=y0;
    t_n=t0;

    % Implement the Euler method
    for i=1:steps
        y_np1=y_n+f(t_n,y_n)*h;
        t_n=t_n+h;
        y_n=y_np1;
    end
    
    y_n_tot{j}=y_n;
    
    disp(sprintf('    For h = %g, y(%g) is approximately %g',h, T(j), y_n));
end

%% Find the exact solution and evaluate function at given T
% This is a First-order Linear ODE. 
syms t
T=[.2 .4 .6 .8 1.0];

t=T;

g_exact = dsolve('Dy = 3-2*t-(0.5)*y','y(0)=1','t'); % using dsovle

y_exact = round(subs(g_exact, t), 4); %plug in value for t in dsovle


%% Find absolute error using |y_exact-y_n_tot|
AbsoluteError = round(abs(y_exact - cell2mat(y_n_tot)),4);


%% Sketch the exact and approximated solutions 
figure(1); 
%plot(t,y_exact, '-ro','Linewidth',2.5); hold on;
plot(T,cell2mat(y_n_tot),'-gs', 'Linewidth',2.5);
legend('Euler','NorthWest');
xlabel('t'); 
ylabel('y'); 
grid on;
p1 = openfig('Rk4_graph.fig');
p2 = openfig('Exact_graph.fig');
%p3 = openfig("Euler_graph.fig");
copyobj(get(gca(p1),'Children'),gca(p2));
