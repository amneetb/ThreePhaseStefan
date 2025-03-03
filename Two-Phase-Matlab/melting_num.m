
% =========================================================================
%              Validation of one-dimensional melting Stefan problem
% =========================================================================

% =========================================================================
%   - Using Finite Difference.
%   - Using implicit approach for solving temperature and explicit 
%   discretization for solving convection term.
%   - Using CUI limiter for solving convection term.
%   - Storing Temperature at cell centers and velocity values at cell
%   faces. 
%   - Using the midpoint time stepping.
%   - T and solid_velocity have three columns:
%               column 1 --> t=n-1
%               column 2 --> t=n
%               column 3 --> t=n+1


% =========================================================================

%%

clc;
clear all;
close all;

% Properties of Solid and Liquid phases
params.cp_s = 910;  %(J/kg.K)
params.k_s  = 211;  % (W/m.K)
params.rho_s  = 2698.72; %(kg/m^3)
params.alpha_s = params.k_s/(params.rho_s*params.cp_s);

params.cp_l = 1042.4;  %(J/kg.K)
params.k_l  = 91.0;    % (W/m.K)
params.rho_l  = 2368;  %(kg/m^3)
params.alpha_l = params.k_l/(params.rho_l*params.cp_l);

params.L_m = 383840; %(J/kg)
params.T_melt = 933.6; %melting temperature
params.T_r = 933.6; %reference temperature
params.Rp = params.rho_s/params.rho_l;

% Right and Left boundary temperature
params.T_right = 298;                         
params.T_left = 2200;

%Simulation setting
params.L = 1; % Length of the domain (m)
params.Nx = 1280; % Number of cells
params.dx = params.L/(params.Nx);
params.start_time = 1.0; % Starting time for numerical simulation
params.Tmax = 5.0; 

x_c = linspace(params.dx/2, params.L-params.dx/2, params.Nx);  % cell center position                                
                                
% Compute the initial interface position and temperature 
% distribution in the domain based on the analytical solution.
[T_analytical, dsdt_analytical, s_analytical, ~] = analytical_sol(params.start_time, params);

% Determine a stable dt for explicit convection
params.cfl = 0.005;
params.dt = params.cfl*(params.dx/dsdt_analytical);
time = [params.start_time:params.dt:params.Tmax]';
nt = length(time); % Number of time steps

% Define interface position and initialize it analytically.
interface_pos = zeros(nt, 1);
interface_pos(1) = s_analytical;


% Define interface velocity and initialize it analytically.
interface_vel = zeros(nt, 1);
interface_vel(1) = dsdt_analytical;

% Define velocity field and initialize it analytically. 
u = zeros(params.Nx, 3);
for i = 1:params.Nx
    if (x_c(i) < interface_pos(1))
        u(i,2) = 0.0; 
        u(i,1) = u(i,2);
    elseif (x_c(i) > interface_pos(1))  
        u(i,2) = (1-(1/params.Rp))*interface_vel(1); 
        u(i,1) = u(i,2);
    else 
        u(i,2) = interface_vel(1);
        u(i,1) = u(i,2);
    end
end


% Define temperature matrix and initialize it analytically.
T = zeros(params.Nx, 3);
T(:,1) = T_analytical; % at n-1
T(:,2) = T_analytical; % at n

%Temperature and velocity at n+1/2. This matrix is only used for calculating CUI limiter.
T_star = zeros(params.Nx+2,1);
u_star = zeros(params.Nx+2,1);

% Time integration loop
Tu = [];
Tc = [];
Td = [];
C_L = params.alpha_l*(params.dt / params.dx^2);
C_S = params.alpha_s*(params.dt / params.dx^2);
for t = 2:nt

    A = zeros(params.Nx, params.Nx);
    b = zeros(params.Nx,1);

    % Estimate interface position.
    interface_pos(t) = interface_pos(t-1) + interface_vel(t-1) * params.dt;

    % Forcing point B (in solid phase)
    index_B = find(x_c >= interface_pos(t), 1);

    % Forcing point B (in liquid phase)
    index_A = index_B - 1; 

    % Extrapolate temperature and velocity to midpoint level n+1/2
     T_star(1:params.Nx,1) = 1.5*T(:,2) - 0.5*T(:,1); 
     u_star(1:params.Nx,1) = 1.5*u(:,2) - 0.5*u(:,1); 
     T_star(params.Nx+1,1) = 2*params.T_right - T_star(params.Nx,1);
     T_star(params.Nx+2,1) = 2*params.T_right - T_star(params.Nx-1,1);
     u_star(params.Nx+1,1) = u_star(params.Nx,1); 
     u_star(params.Nx+2,1) = u_star(params.Nx,1); 
     
    % Phase A : Liquid 
    for i = 1:index_A

        if (i > 1 && i<index_A)
    
            A(i, i-1) = -0.5 * C_L;
            A(i, i) = 1 +  C_L;
            A(i, i+1) = -0.5 * C_L;
    
            b(i) = T(i,2) + 0.5 * C_L * (T(i+1,2) -2*T(i,2) + T(i-1,2));
    
        elseif (i == 1)
    
            A(i, i) = 1 + 1.5 * C_L;
            A(i, i+1) = -0.5 * C_L;
    
            b(i) = T(i,2) + 2 * C_L * params.T_left + 0.5 * C_L * (T(i+1,2) -3 * T(i,2));
    
        elseif (i == index_A)
    
            % Calculate temperature at forcing point A. 
            s = interface_pos(t);
            a1 = ((s - x_c(i-1))*(s - x_c(i-2)))/((x_c(i)-x_c(i-1))*(x_c(i)-x_c(i-2)));
            a2 = ((s - x_c(i))*(s - x_c(i-2)))/((x_c(i-1)-x_c(i))*(x_c(i-1)-x_c(i-2)));
            a3 = ((s - x_c(i))*(s - x_c(i-1)))/((x_c(i-2)-x_c(i))*(x_c(i-2)-x_c(i-1)));
            
            %Applying zero Jump condition for temperature in Phase A (liquid)
            %to the nearest interface cell
            A(i, i)   = a1;
            A(i, i-1) = a2;
            A(i, i-2) = a3;
            
            b(i) = params.T_melt;
    
        end
   end

    % Phase B : Solid 
   for i = index_B:params.Nx

        if (i > index_B && i < params.Nx)
           
            ur_face = 0.5*(u_star(i+1,1) + u_star(i,1));
            if (ur_face >= 0.0)
                Tu = T_star(i-1,1);
                Tc = T_star(i,1);
                Td = T_star(i+1,1);
            else
                Tu = T_star(i+2,1);
                Tc = T_star(i+1,1);
                Td = T_star(i,1);
            end    
            Tr = cui_scheme(Tc,Tu,Td);

            ul_face = 0.5*(u_star(i-1,1) + u_star(i,1));
            if (ul_face >= 0.0)
                Tu = T_star(i-2,1);
                Tc = T_star(i-1,1);
                Td = T_star(i,1);
            else
                Tu = T_star(i+1,1);
                Tc = T_star(i,1);
                Td = T_star(i-1,1);
            end    
            Tl = cui_scheme(Tc,Tu,Td);

            A(i, i) = 1 + C_S ;
            A(i, i+1) = -0.5 * C_S ;
            A(i, i-1) = - 0.5 * C_S;
            b(i) = T(i,2) + 0.5 * C_S * (T(i+1,2) -2*T(i,2) + T(i-1,2)) - u_star(i,1)*(params.dt/params.dx)*(Tr - Tl);

      elseif (i == params.Nx)

            ur_face = 0.5*(u_star(i+1,1) + u_star(i,1));
            if (ur_face >= 0.0)
                Tu = T_star(i-1,1);
                Tc = T_star(i,1);
                Td = T_star(i+1,1);
            else
                Tu = T_star(i+2,1);
                Tc = T_star(i+1,1);
                Td = T_star(i,1);
            end    
            Tr = cui_scheme(Tc,Tu,Td);

            ul_face = 0.5*(u_star(i-1,1) + u_star(i,1));
            if (ul_face >= 0.0)
                Tu = T_star(i-2,1);
                Tc = T_star(i-1,1);
                Td = T_star(i,1);
            else
                Tu = T_star(i+1,1);
                Tc = T_star(i,1);
                Td = T_star(i-1,1);
            end    
            Tl = cui_scheme(Tc,Tu,Td);

            A(i,i) = 1 + 1.5*C_S ;
            A(i,i-1) = - 0.5 * C_S;
            b(i) = T(i,2) + 2 * C_S * params.T_right + 0.5 * C_S * (-3*T(i,2) + T(i-1,2)) - u_star(i,1)*(params.dt/params.dx)*(Tr - Tl);

      elseif (i == index_B)
    
           % Calculate temperature at forcing point B. 
            s = interface_pos(t);
            b1 = ((s - x_c(i+1))*(s - x_c(i+2)))/((x_c(i)-x_c(i+1))*(x_c(i)-x_c(i+2)));
            b2 = ((s - x_c(i))*(s - x_c(i+2)))/((x_c(i+1)-x_c(i))*(x_c(i+1)-x_c(i+2)));
            b3 = ((s - x_c(i))*(s - x_c(i+1)))/((x_c(i+2)-x_c(i))*(x_c(i+2)-x_c(i+1)));
            
            % Applying zero jump condition for temperature in Phase B (solid)
            % to the nearest interface cell
            A(i, i)   = b1;
            A(i, i+1) = b2;
            A(i, i+2) = b3;
    
            b(i) = params.T_melt;
        end

   end

   % Solve the linear system. 
   T(:,3) = A \ b;


   % Compute dTdx at the interface in phase A
   s = interface_pos(t);
   x1 = x_c(index_A) - s;
   x2 = x_c(index_A-1) - s;
   x3 = x_c(index_A-2) - s;
   
   c = compute_c_coefs(x1, x2, x3); 
   dTdx_liq = c(1)*T(index_A,3) + c(2)*T(index_A-1,3) + c(3)*T(index_A-2,3);
   
   % Compute dTdx at the interface in phase B
   x1 = x_c(index_B) - s;
   x2 = x_c(index_B+1) - s;
   x3 = x_c(index_B+2) - s;
   
   c = compute_c_coefs(x1, x2, x3);
   dTdx_sol = c(1)*T(index_B,3) + c(2)*T(index_B+1,3) + c(3)*T(index_B+2,3);
%% Solving interface velocity by considering kinetic energy.
%   % calculate interface velocity and correct interface position.
%   f = @(vel) params.rho_l * vel*(params.L_m + 0.5 * (1 - 1 / params.Rp^2) * vel^2) ...
%            - (params.k_s * dTdx_sol - params.k_l * dTdx_liq);
% 
%   % Initial guess for interface_vel(t)
%   initial_guess = interface_vel(t-1); 
% 
%   % options = optimoptions('fsolve', 'Display', 'iter'); % Display iteration details (optional)
%   [dsdt_new, fval] = fsolve(f, initial_guess);
%   interface_vel(t) = dsdt_new;
  %%
  interface_vel(t)=(params.k_s * dTdx_sol - params.k_l * dTdx_liq)/(params.rho_l*((params.cp_l-params.cp_s)*(params.T_melt-params.T_r)+params.L_m));
  interface_pos(t) = interface_pos(t-1) + 0.5*(interface_vel(t)+interface_vel(t-1))*params.dt;

  % Update velocity field
  for i = 1:params.Nx
     if (x_c(i) < interface_pos(t))
        u(i,3) = 0; 
     elseif (x_c(i) > interface_pos(t))
        u(i,3) = (1-(1/params.Rp))*interface_vel(t); 
     else
        u(i,3)=interface_vel(t);
     end
  end

   % Shift the temperature and velocity data for the next time step
   T(:,1) = T(:,2); %time at n-1
   T(:,2) = T(:,3); %time at n

   u(:,1) = u(:,2);
   u(:,2) = u(:,3);

   if (interface_pos(t) > params.L - 3*params.dx)
       break
   end
end

%% Find the analytical solution at the end of the simulation
[T_analytical, dsdt_analytical, s_analytical, lambda] = analytical_sol(time(end), params);
dsdt_analytical_vec = lambda*sqrt(params.alpha_s./time);
s_analytical_vec = 2*lambda*(sqrt(time*params.alpha_s));

%% Validation of temperature profile with analytical solution

figure(1)
plot(x_c, T_analytical, 'LineWidth', 2);
hold on
plot(x_c, T(:,3),'--', 'LineWidth', 2);
xlabel('Domain length (m)');
ylabel('Temperature (K)');
title('Temperature Profile');
legend ('Analytical Solution','Numerical Solution')

%% Validation of interface velocity with analytical solution

figure(2)
plot(time, dsdt_analytical_vec , 'LineWidth', 2);
hold on
plot(time, interface_vel,'--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend ('Analytical Solution','Numerical Solution')


%% Validation of interface position with analytical solution

figure(3)
plot(time, s_analytical_vec , 'LineWidth', 2);
hold on;
plot(time, interface_pos,'--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Interface position (m)');
legend ('Analytical Solution','Numerical Solution')


%% calculate L2 norm of interface position and temperature profile.
 
L2_norm_s =  sqrt(sum((s_analytical_vec - interface_pos).^2*params.dt));
L2_norm_temp = sqrt(sum((T_analytical - T(:,3)).^2*params.dx));

% Display the result
fprintf('The L2 norm error of interface position at time %f is: %f\n', time(end), L2_norm_s);
fprintf('The L2 norm error of temperature profile at time %f is: %f\n', time(end), L2_norm_temp);

%% Storing results
fileName = sprintf('Run%d.mat', params.Nx);
% Save the variables into the file
save(fileName, 'L2_norm_s', 'L2_norm_temp', 'time','interface_pos', 's_analytical_vec','x_c', 'T_analytical','T');


%%
function c = compute_c_coefs(x1,x2,x3)

   coef_mat = [1     1     1;
               x1    x2    x3;
               x1^2  x2^2  x3^2];
   d = [0; 1; 0];
   c = coef_mat \ d;

end