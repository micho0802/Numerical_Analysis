%% MIT18086_NAVIERSTOKES
%    Solves the incompressible Navier-Stokes equations in a
%    rectangular domain with prescribed velocities along the
%    boundary. The solution method is finite differencing on
%    a staggered grid with implicit diffusion and a Chorin
%    projection method for the pressure.
%    Visualization is done by a colormap-isoline plot for
%    pressure and normalized quiver and streamline plot for
%    the velocity field.
%    The standard setup solves a lid driven cavity problem.

% 07/2007 by Benjamin Seibold
%            http://www-math.mit.edu/~seibold/
% Feel free to modify for teaching and learning.
%-----------------------------------------------------------------------
Re = 4000; % Reynolds number
dt = 1e-3; % time step
tf = 4; % final time
lx = 1; % width of box
ly = 1; % height of box
nx = 257; % number of x-gridpoints
ny = 257; % number of y-gridpoints
nsteps = 10; % number of steps with graphic output

%%
nt = ceil(tf/dt); % number of iterations
x = linspace(0, lx, nx+1); % x vector
hx = lx/nx; % spatial step size in x
y = linspace(0, ly, ny+1); % y vector
hy = ly/ny; % spatial step size in y

%%
% initial conditions
U = zeros(nx-1, ny); % initial U matrix
V = zeros(nx, ny-1); % initial V matrix
% boundary conditions
uN = x*0+1; % u - North
vN = avg(x)*0; % v - North
uS = x*0; % u - South
vS = avg(x)*0; % v - South
uW = avg(y)*0; % u - West
vW = y*0; % v - West
uE = avg(y)*0; % u - East
vE = y*0; % v - East

%%
% boundary condition matrices
Ubc = dt/Re * ([2*uS(2:end-1)', zeros(nx-1,ny-2), 2*uN(2:end-1)'] / hx^2 + ...
               [uW; zeros(nx-3,ny); uE] / hy^2);
Vbc = dt/Re * ([vS', zeros(nx,ny-3), vN'] / hx^2 + ...
               [2*vW(2:end-1); zeros(nx-2,ny-1); 2*vE(2:end-1)] / hy^2);

fprintf('initialization')

% Discretize the laplacian
Lu = speye((nx-1)*ny) + dt/Re * (kron(speye(ny), lap1d(nx-1,hx,2)) + ...
                                 kron(lap1d(ny,hy,3), speye(nx-1)));
peru = symamd(Lu);
Ru = chol(Lu(peru,peru));
Rut = Ru';

Lv = speye(nx*(ny-1)) + dt/Re * (kron(speye(ny-1), lap1d(nx,hx,3)) + ...
                                 kron(lap1d(ny-1,hy,2), speye(nx)));
perv = symamd(Lv);
Rv = chol(Lv(perv,perv));
Rvt = Rv';

Lp = kron(speye(ny), lap1d(nx,hx,1)) + kron(lap1d(ny,hy,1), speye(nx));
Lp(1,1) = 3/2 * Lp(1,1);
perp = symamd(Lp);
Rp = chol(Lp(perp,perp));
Rpt = Rp';

Lq = kron(speye(ny-1), lap1d(nx-1,hx,2)) + kron(lap1d(ny-1,hy,2), speye(nx-1));
perq = symamd(Lq);
Rq = chol(Lq(perq,perq));
Rqt = Rq';

fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')


%%
for k = 1:nt
    % treat nonlinear terms
    gamma = min(1.2 * dt * max(max(max(abs(U))) / hx, max(max(abs(V))) / hy), 1);
    Ue = [uW; U; uE];
    Ue = [2*uS'-Ue(:,1), Ue, 2*uN'-Ue(:,end)];
    Ve = [vS', V, vN'];
    Ve = [2*vW-Ve(1,:); Ve; 2*vE-Ve(end,:)];
    Ua = avg(Ue')';
    Ud = diff(Ue')'/2;
    Va = avg(Ve);
    Vd = diff(Ve)/2;
    UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd) / hx;
    UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')' / hy;

    Ua = avg(Ue(:,2:end-1));
    Ud = diff(Ue(:,2:end-1)) / 2;
    Va = avg(Ve(2:end-1,:)')';
    Vd = diff(Ve(2:end-1,:)')' / 2;

    U2x = diff(Ua.^2-gamma*abs(Ua).*Ud) / hx;
    V2y = diff((Va.^2-gamma*abs(Va).*Vd)')' / hy;

    U = U - dt * (UVy(2:end-1,:) + U2x);
    V = V - dt * (UVx(:,2:end-1) + V2y);

    % implicit viscosity
    rhs = reshape(U + Ubc, [], 1);
    u(peru) = Ru \ (Rut \ rhs(peru));
    U = reshape(u, nx-1, ny);

    rhs = reshape(V + Vbc, [], 1);
    v(perv) = Rv \ (Rvt \ rhs(perv));
    V = reshape(v, nx, ny-1);

    % pressure correction
    rhs = reshape(diff([uW;U;uE]) / hx + diff([vS' V vN']')' / hy, [], 1);
    p(perp) = -Rp \ (Rpt \ rhs(perp));
    P = reshape(p, nx, ny);
    U = U - diff(P) / hx;
    V = V - diff(P')' / hy;

    if isnan(U)
        break
    end
    
    U_tensor = zeros(size(U, 1), size(U, 2), nt);
    V_tensor = zeros(size(V, 1), size(V, 2), nt);
    P_tensor = zeros(size(P, 1), size(P, 2), nt);
    % store tensors
    U_tensor(:,:,k) = U;
    V_tensor(:,:,k) = V;
    P_tensor(:,:,k) = P;

    % visualization
    if floor(25*k/nt) > floor(25*(k-1)/nt)
        fprintf("\nTime step: %d\n", k)
    end

    if k == 1 || floor(nsteps * k / nt) > floor(nsteps * (k-1) / nt)
        rhs = reshape(diff(U')' / hy - diff(V) / hx, [], 1);
        q(perq) = Rq \ (Rqt \ rhs(perq));
        Q = zeros(nx+1, ny+1);
        Q(2:end-1,2:end-1) = reshape(q, nx-1, ny-1);

        % visualization
        clf;
        contourf(avg(x), avg(y), P', 20, 'w-');
        hold on;
        contour(x, y, Q', 20, 'k-');

        Ue = [uS' avg([uW;U;uE]')' uN'];
        Ve = [vW; avg([vS' V vN']); vE];
        Len = sqrt(Ue.^2 + Ve.^2 + eps);
        quiver(x, y, (Ue ./ Len)', (Ve ./ Len)', .4, 'k-')
        hold off;
        axis equal, axis([0 lx 0 ly]);

        p = sort(p);
        disp(['min(p): ', num2str(min(p)), ', max(p): ', num2str(max(p))]);
        clim([min(p), max(p)]);
        title(sprintf('Re = %0.1g   t = %0.2g', Re, k * dt))
        drawnow
        %Store the Q tensor
        Q_tensor = zeros(size(Q, 1), size(Q, 2), nt);
        Q_tensor(:,:,k) = Q;
    end
end

% store tensors in h5 file
h5_filename = 'navier_stokes_simulation_Re4000_256256.h5';
h5create(h5_filename, '/U_tensor', [size(U, 1), size(U, 2), nt], 'Datatype', 'double');
h5create(h5_filename, '/V_tensor', [size(V, 1), size(V, 2), nt], 'Datatype', 'double');
h5create(h5_filename, '/P_tensor', [size(P, 1), size(P, 2), nt], 'Datatype', 'double');
h5create(h5_filename, '/Q_tensor', [size(Q, 1), size(Q, 2), nt], 'Datatype', 'double');
h5write(h5_filename, '/U_tensor', U_tensor);
h5write(h5_filename, '/V_tensor', V_tensor);
h5write(h5_filename, '/P_tensor', P_tensor);
h5write(h5_filename, '/Q_tensor', Q_tensor);