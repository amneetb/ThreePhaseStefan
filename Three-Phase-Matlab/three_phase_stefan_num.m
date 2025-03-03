
% =========================================================================
%        Validation of one-dimensional Three-Phase Stefan problem
% =========================================================================

% =========================================================================
%   - Three phases are as follows:
%               phase A  --> Vapor (abuts the left wall)
%               phase B  --> Liquid 
%               Phase C  --> Solid (abuts the right wall)
%   - Uses second-order finite differences.
%   - Uses implicit discretization of the conduction term and explicit 
%   discretization of the convection term.
%   - Uses CUI limiter for upwinding the convection term.
%   - Stores temperature at cell centers and velocity at cell
%   faces. 
%   - Uses midpoint time stepping scheme.
%   - T and solid_velocity have three columns:
%               column 1 --> t=n-1
%               column 2 --> t=n
%               column 3 --> t=n+1


% =========================================================================

% %%

clc;
clear all;
close all;

% Thermophysical properties of the vapor, liquid and solid phases
params.cp_s = 910;  %(J/kg.K)
params.k_s  = 211;  % (W/m.K)
params.rho_s  = 2698.72; %(kg/m^3)
params.alpha_s = params.k_s/(params.rho_s*params.cp_s);

params.cp_l = 1042.4;  %(J/kg.K)
params.k_l  = 91.0;    % (W/m.K)
params.rho_l  = 2368;  %(kg/m^3)
params.alpha_l = params.k_l/(params.rho_l*params.cp_l);

params.cp_v = 770.69;  %(J/kg.K)
params.k_v  = 115.739;    % (W/m.K)
params.rho_v  = 0.08644;  %(kg/m^3)
params.alpha_v = params.k_v/(params.rho_v*params.cp_v);

params.L_m = 383840; % Latent heat melting (J/kg)
params.L_v = 9462849.518; % Latent heat vaporization (J/kg)
params.T_r = 933.6; % Reference temperature
params.T_m = 933.6; % Melting temperature
params.T_v = 2767; % Boiling temperature

% Density ratios
params.R_rho_LS = params.rho_l/params.rho_s;
params.R_rho_VL = params.rho_v/params.rho_l;
params.R_rho_VS = params.rho_v/params.rho_s;


% Right and Left boundary temperature
params.T_0  = 298;  %(K)                       
params.T_inf  = 5000; %(K)

%Simulation setting
params.L = 1;  % Length of the domain (m)
params.Nx = 1280;  % Number of cells
params.dx = params.L/(params.Nx);
params.start_time = 1;  % Start time
params.Tmax = 5.0;  % End time

x_c = linspace(params.dx/2, params.L-params.dx/2, params.Nx);  % cell center position                                
                                
% Initialize temperature, interface position and velocity using analytical solution
% at start_time
[T_analytical, dsdt_analytical_LV,dsdt_analytical_LS, s_analytical_LV,s_analytical_LS, ~] = analytical_sol(params.start_time, params);


% dt for explicit convection
params.cfl = 1e-4;
params.dt = params.cfl*(params.dx/max(dsdt_analytical_LS,dsdt_analytical_LV));
time = [params.start_time:params.dt:params.Tmax]';
nt = length(time);                           

% Define interface position and initialize it analytically.
interface1_pos = zeros(nt, 1);
interface1_pos(1) = s_analytical_LS;

interface2_pos = zeros(nt, 1);
interface2_pos(1) = s_analytical_LV;


% Define interface velocity and initialize it analytically.
interface1_vel = zeros(nt, 1);
interface1_vel(1) = dsdt_analytical_LS;

interface2_vel = zeros(nt, 1);
interface2_vel(1) = dsdt_analytical_LV;


% Define velocity field and initialize it analytically. 
u = zeros(params.Nx, 3);
for i = 1:params.Nx
    if (x_c(i) < interface2_pos(1))
        u(i,2) = 0.0; 
        u(i,1) = u(i,2);
    elseif (x_c(i) > interface2_pos(1) && x_c(i) < interface1_pos(1))  
        u(i,2) = (1-params.R_rho_VL)*interface2_vel(1); 
        u(i,1) = u(i,2);
    elseif (x_c(i) > interface1_pos(1))
        u(i,2) = (params.R_rho_LS-params.R_rho_VS)*interface2_vel(1) + (1-params.R_rho_LS)*interface1_vel(1);
        u(i,1) = u(i,2);
    elseif x_c(i) == interface1_pos(1)
        u(i,2)=interface1_vel(1);
        u(i,1) = u(i,2);
    elseif x_c(i) == interface2_pos(1)
        u(i,2)=interface2_vel(1);
        u(i,1) = u(i,2);    
    end
end


% Define temperature matrix and initialize it analytically.
T = zeros(params.Nx, 3);
T(:,1) = T_analytical; % at n-1
T(:,2) = T_analytical; % at n

%Temperature and velocity at n+1/2. This matrix is only used for evaluating the CUI limiter.
T_star = zeros(params.Nx+2,1);
u_star = zeros(params.Nx+2,1);

% Time integration loop
Tu = [];
Tc = [];
Td = [];
C_L = params.alpha_l*(params.dt / params.dx^2);
C_S = params.alpha_s*(params.dt / params.dx^2);
C_V = params.alpha_v*(params.dt / params.dx^2);

for t = 2:nt

    A = zeros(params.Nx, params.Nx);
    b = zeros(params.Nx,1);

    % Estimate the two interface positions.
    interface2_pos(t) = interface2_pos(t-1) + interface2_vel(t-1) * params.dt/2;
    interface1_pos(t) = interface1_pos(t-1) + interface1_vel(t-1) * params.dt/2;


    % Forcing point B in the liquid phase across the LV interface
    index_B = find(x_c >= interface2_pos(t), 1);
    % Forcing point A in the vapor phase before the LV interface
    index_A = index_B - 1; 

    % Forcing point D in the solid phase across the LS interface
    index_D = find(x_c >= interface1_pos(t), 1);
    % Forcing point C in the liquid phase before the LS interface
    index_C = index_D - 1; 


    % Extrapolate temperature and velocity to midpoint level n+1/2
     T_star(1:params.Nx,1) = 1.5*T(:,2) - 0.5*T(:,1); 
     u_star(1:params.Nx,1) = 1.5*u(:,2) - 0.5*u(:,1); 
     T_star(params.Nx+1,1) = 2*params.T_0 - T_star(params.Nx,1);
     T_star(params.Nx+2,1) = 2*params.T_0 - T_star(params.Nx-1,1);
     u_star(params.Nx+1,1) = u_star(params.Nx,1); 
     u_star(params.Nx+2,1) = u_star(params.Nx,1); 
     
    % Phase A : Vapor 
    for i = 1:index_A

        if (i > 1 && i<index_A)
    
            A(i, i-1) = -0.5 * C_V;
            A(i, i) = 1 +  C_V;
            A(i, i+1) = -0.5 * C_V;
    
            b(i) = T(i,2) + 0.5 * C_V * (T(i+1,2) -2*T(i,2) + T(i-1,2));
    
        elseif (i == 1)
    
            A(i, i) = 1 + 1.5 * C_V;
            A(i, i+1) = -0.5 * C_V;
    
            b(i) = T(i,2) + 2 * C_V * params.T_inf + 0.5 * C_V * (T(i+1,2) -3 * T(i,2));
    
        elseif (i == index_A)
    
            % Calculate temperature at forcing point A. 
            s = interface2_pos(t);
            a1 = ((s - x_c(i-1))*(s - x_c(i-2)))/((x_c(i)-x_c(i-1))*(x_c(i)-x_c(i-2)));
            a2 = ((s - x_c(i))*(s - x_c(i-2)))/((x_c(i-1)-x_c(i))*(x_c(i-1)-x_c(i-2)));
            a3 = ((s - x_c(i))*(s - x_c(i-1)))/((x_c(i-2)-x_c(i))*(x_c(i-2)-x_c(i-1)));
             
            % Apply zero jump condition for temperature in Phase A (vapor)
            % to the nearest interface cell
            A(i, i)   = a1;
            A(i, i-1) = a2;
            A(i, i-2) = a3;
            
            b(i) = params.T_v;
        end
    end

    % Phase B : Liquid 
    for i = index_B:index_C

        if (i > index_B && i < index_C)
           
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

            A(i, i) = 1 + C_L ;
            A(i, i+1) = -0.5 * C_L ;
            A(i, i-1) = - 0.5 * C_L;
            b(i) = T(i,2) + 0.5 * C_L * (T(i+1,2) -2*T(i,2) + T(i-1,2)) - u_star(i,1)*(params.dt/params.dx)*(Tr - Tl);

      elseif (i == index_B)
    
           % Calculate temperature at forcing point B. 
            s = interface2_pos(t);
            b1 = ((s - x_c(i+1))*(s - x_c(i+2)))/((x_c(i)-x_c(i+1))*(x_c(i)-x_c(i+2)));
            b2 = ((s - x_c(i))*(s - x_c(i+2)))/((x_c(i+1)-x_c(i))*(x_c(i+1)-x_c(i+2)));
            b3 = ((s - x_c(i))*(s - x_c(i+1)))/((x_c(i+2)-x_c(i))*(x_c(i+2)-x_c(i+1)));
            
            % Applying zero jump condition for temperature in Phase B (liquid)
            % to the nearest interface cell
            A(i, i)   = b1;
            A(i, i+1) = b2;
            A(i, i+2) = b3;
    
            b(i) = params.T_v;

       elseif (i == index_C)
    
           % Calculate temperature at forcing point C. 
            s = interface1_pos(t);
            a1 = ((s - x_c(i-1))*(s - x_c(i-2)))/((x_c(i)-x_c(i-1))*(x_c(i)-x_c(i-2)));
            a2 = ((s - x_c(i))*(s - x_c(i-2)))/((x_c(i-1)-x_c(i))*(x_c(i-1)-x_c(i-2)));
            a3 = ((s - x_c(i))*(s - x_c(i-1)))/((x_c(i-2)-x_c(i))*(x_c(i-2)-x_c(i-1)));
            
            %Applying zero jump condition for temperature in Phase B (liquid)
            %to the nearest interface cell
            A(i, i)   = a1;
            A(i, i-1) = a2;
            A(i, i-2) = a3;
    
            b(i) = params.T_m;
        end

    end

    for i = index_D:params.Nx

        if (i > index_D && i < params.Nx)
           
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
            b(i) = T(i,2) + 2 * C_S * params.T_0 + 0.5 * C_S * (-3*T(i,2) + T(i-1,2)) - u_star(i,1)*(params.dt/params.dx)*(Tr - Tl);

      elseif (i == index_D)
    
           % Calculate temperature at forcing point D. 
            s = interface1_pos(t);
            b1 = ((s - x_c(i+1))*(s - x_c(i+2)))/((x_c(i)-x_c(i+1))*(x_c(i)-x_c(i+2)));
            b2 = ((s - x_c(i))*(s - x_c(i+2)))/((x_c(i+1)-x_c(i))*(x_c(i+1)-x_c(i+2)));
            b3 = ((s - x_c(i))*(s - x_c(i+1)))/((x_c(i+2)-x_c(i))*(x_c(i+2)-x_c(i+1)));
            
            % Applying zero jump condition for temperature in Phase C (solid)
            % to the nearest interface cell
            A(i, i)   = b1;
            A(i, i+1) = b2;
            A(i, i+2) = b3;
    
            b(i) = params.T_m;
        end

   end

   % Solve the linear system. 
   T(:,3) = A \ b;


   % Compute dTdx at the interface 2 in phase A
   s = interface2_pos(t);
   x1 = x_c(index_A) - s;
   x2 = x_c(index_A-1) - s;
   x3 = x_c(index_A-2) - s;
   
   c = compute_c_coefs(x1, x2, x3); 
   dTdx_vap = c(1)*T(index_A,3) + c(2)*T(index_A-1,3) + c(3)*T(index_A-2,3);

   % Compute dTdx at the interface 2 in phase B
   s = interface2_pos(t);
   x1 = x_c(index_B) - s;
   x2 = x_c(index_B+1) - s;
   x3 = x_c(index_B+2) - s;
   
   c = compute_c_coefs(x1, x2, x3);
   dTdx_liq_interface2 = c(1)*T(index_B,3) + c(2)*T(index_B+1,3) + c(3)*T(index_B+2,3);

   % Compute dTdx at the interface 1 in phase C
   s = interface1_pos(t);
   x1 = x_c(index_C) - s;
   x2 = x_c(index_C-1) - s;
   x3 = x_c(index_C-2) - s;
   
   c = compute_c_coefs(x1, x2, x3); 
   dTdx_liq_interface1 = c(1)*T(index_C,3) + c(2)*T(index_C-1,3) + c(3)*T(index_C-2,3);
   
   % Compute dTdx at the interface 1 in phase D
   s = interface1_pos(t);
   x1 = x_c(index_D) - s;
   x2 = x_c(index_D+1) - s;
   x3 = x_c(index_D+2) - s;
   
   c = compute_c_coefs(x1, x2, x3);
   dTdx_sol = c(1)*T(index_D,3) + c(2)*T(index_D+1,3) + c(3)*T(index_D+2,3);


   %% Calculation of liquid-vapor (interface 2) and liquid-solid (interface 1) without adding kinetic energy.
   interface2_vel(t) =(params.k_l*dTdx_liq_interface2-params.k_v*dTdx_vap)/(params.rho_v*((params.cp_v-params.cp_l)*(params.T_v-params.T_r)+params.L_v));
   fac = (params.cp_l-params.cp_s)*(params.T_m-params.T_r)+params.L_m;
   interface1_vel(t) =(1/params.rho_l)*((params.k_s*dTdx_sol-params.k_l*dTdx_liq_interface1)/fac+(params.rho_l-params.rho_v)*interface2_vel(t));

  interface1_pos(t) = interface1_pos(t-1) + interface1_vel(t)*params.dt;
  interface2_pos(t) = interface2_pos(t-1) + interface2_vel(t)*params.dt;

  for i = 1:params.Nx
    if (x_c(i) < interface2_pos(t))
        u(i,3) = 0.0; 
    elseif (x_c(i) > interface2_pos(t) && x_c(i) < interface1_pos(t))  
        u(i,3) = (1-params.R_rho_VL)*interface2_vel(t); 
    elseif (x_c(i) > interface1_pos(t))
        u(i,3) = (params.R_rho_LS-params.R_rho_VS)*interface2_vel(t) + (1-params.R_rho_LS)*interface1_vel(t);
    elseif x_c(i) == interface1_pos(t)
        u(i,3)=interface1_vel(t);
    elseif x_c(i) == interface2_pos(t)
        u(i,3)=interface2_vel(t);
    end
 end


   % Shift the temperature and velocity data for the next time step
   T(:,1) = T(:,2); %time at n-1
   T(:,2) = T(:,3); %time at n

   u(:,1) = u(:,2);
   u(:,2) = u(:,3);


end

%% Find the analytical solution at the end of the simulation
[T_analytical, dsdt_analytical_LV,dsdt_analytical_LS, s_analytical_LV,s_analytical_LS, lambda, beta] = analytical_sol(time(end), params);

dsdt_analytical_vec_LS = beta*sqrt(params.alpha_s./time);
s_analytical_vec_LS = 2*beta*(sqrt(time*params.alpha_s));

dsdt_analytical_vec_LV = lambda*sqrt(params.alpha_l./time);
s_analytical_vec_LV = 2*lambda*(sqrt(time*params.alpha_l));

%% Validation of temperature profile with analytical solution

figure(1)
plot(x_c, T_analytical, 'LineWidth', 2);
hold on
plot(x_c, T(:,3),'--', 'LineWidth', 2);
xlabel('Domain length (m)');
ylabel('Temperature (K)');
title('Temperature Profile');
legend ('Analytical Solution','Numerical Solution')
print('Temperature','-dpng','-r1000');  


%% Validation of interface velocity with analytical solution

figure(2)
hold on

plot(time, dsdt_analytical_vec_LS , 'LineWidth', 2, 'DisplayName', 'U^*_1 (Analytical solution)');
plot(time, dsdt_analytical_vec_LV , 'LineWidth', 2, 'DisplayName', 'U^*_2 (Analytical solution)');

plot(time, interface1_vel,'--', 'LineWidth', 2, 'DisplayName', 'U^*_1 (Numerical solution)');
plot(time, interface2_vel,'--', 'LineWidth', 2, 'DisplayName', 'U^*_2 (Numerical solution)');

hold off
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Location', 'northwest');
box on;

print('Interface_velocity','-dpng','-r1000');  


%% Validation of interface position with analytical solution

figure(3)
hold on;

plot(time, s_analytical_vec_LS , 'LineWidth', 2, 'DisplayName', 'X^*_1 (Analytical solution)');
plot(time, s_analytical_vec_LV , 'LineWidth', 2, 'DisplayName', 'X^*_2 (Analytical solution)');

plot(time, interface1_pos,'--', 'LineWidth', 2, 'DisplayName', 'X^*_1 (Numerical solution)');
plot(time, interface2_pos,'--', 'LineWidth', 2, 'DisplayName', 'X^*_2 (Numerical solution)');

hold off
xlabel('Time (s)');
ylabel('Interface position (m)');
legend('Location', 'northwest');
box on;

print('Interface_position','-dpng','-r1000');  


%% Calculate L2 norm of interface position 2 & 1 and temperature profile.
 
L2_norm_s_LV =  sqrt(sum((s_analytical_vec_LV - interface2_pos).^2*params.dt));
L2_norm_s_LS =  sqrt(sum((s_analytical_vec_LS - interface1_pos).^2*params.dt));
L2_norm_temp = sqrt(sum((T_analytical - T(:,3)).^2*params.dx));

% Display the result
fprintf('The L2 norm error of interface position 2 at time %f is: %f\n', time(end), L2_norm_s_LV);
fprintf('The L2 norm error of interface position 1 at time %f is: %f\n', time(end), L2_norm_s_LS);
fprintf('The L2 norm error of temperature profile at time %f is: %f\n', time(end), L2_norm_temp);

%% Storing results
fileName = sprintf('Run%d.mat', params.Nx);
% Save the variables into the file
save(fileName, 'L2_norm_s_LV','L2_norm_s_LS', 'L2_norm_temp', 'time','interface2_pos','interface1_pos', 's_analytical_vec_LV','s_analytical_vec_LS','x_c', 'T_analytical','T');


%%
function c = compute_c_coefs(x1,x2,x3)

   coef_mat = [1     1     1;
               x1    x2    x3;
               x1^2  x2^2  x3^2];
   d = [0; 1; 0];
   c = coef_mat \ d;

end

