function [T, dsdt_analytical_LV, dsdt_analytical_LS, s_analytical_LV, s_analytical_LS, lambda_sol, beta_sol] = analytical_sol(time, params)

T = zeros(params.Nx,1);

diffusivity_ratio_LS = params.alpha_l/params.alpha_s;
diffusivity_ratio_LV = params.alpha_l/params.alpha_v;


dT_vinf = params.T_v - params.T_inf;
dT_m0 = params.T_m - params.T_0;
dT_vm = params.T_v - params.T_m;

% Define the system of equations as a function
    func = @(x) [
        % Equation 1: Liquid-vapor interface
        x(1) * sqrt(params.alpha_l) * params.rho_v * ((params.cp_v - params.cp_l) * (params.T_v - params.T_r) + params.L_v) + ...
        params.k_v * (dT_vinf * exp(-(x(1)^2 * diffusivity_ratio_LV)) / (sqrt(pi * params.alpha_v) * erf(x(1) * sqrt(params.alpha_l / params.alpha_v)))) - ...
        params.k_l * (dT_vm * exp(-(x(1) * params.R_rho_VL)^2) / (sqrt(pi * params.alpha_l) * (erf(x(1) * params.R_rho_VL) - erf(x(2) * sqrt(params.alpha_s / params.alpha_l) - x(1) * (1 - params.R_rho_VL)))));

        % Equation 2: Liquid-solid interface
        params.k_l * (dT_vm * exp(-(x(2) * sqrt(params.alpha_s / params.alpha_l) - x(1) * (1 - params.R_rho_VL))^2) / ...
        (sqrt(pi * params.alpha_l) * (erf(x(1) * params.R_rho_VL) - erf(x(2) * sqrt(params.alpha_s / params.alpha_l) - x(1) * (1 - params.R_rho_VL))))) + ...
        params.k_s * (dT_m0 * exp(-(x(2) * params.R_rho_LS - x(1) * sqrt(diffusivity_ratio_LS) * (params.R_rho_LS - params.R_rho_VS))^2) / ...
        (sqrt(pi * params.alpha_s) * (erfc(x(2) * params.R_rho_LS - x(1) * sqrt(diffusivity_ratio_LS) * (params.R_rho_LS - params.R_rho_VS))))) + ...
        (params.rho_l * x(2) * sqrt(params.alpha_s) - (params.rho_l - params.rho_v) * x(1) * sqrt(params.alpha_l)) * ...
        ((params.cp_l - params.cp_s) * (params.T_m - params.T_r) + params.L_m);
    ];

    % Initial guesses for lambda and beta
    initialGuess = [0.0012, 0.4304];

    % Solve the system of equations using fsolve
    options = optimset('Display', 'off');  % Optional: to suppress fsolve output
    sol = fsolve(func, initialGuess, options);

    % Extract the solutions for lambda and beta
    lambda_sol = sol(1);
    beta_sol = sol(2);

    %% Calculation of LV_interface and LS_interface velocity with the respected position

    dsdt_analytical_LS = beta_sol*sqrt(params.alpha_s/time);
    s_analytical_LS = 2*beta_sol*(sqrt(time*params.alpha_s));

    dsdt_analytical_LV = lambda_sol*sqrt(params.alpha_l/time);
    s_analytical_LV = 2*lambda_sol*(sqrt(time*params.alpha_l));

    %%

% Temperature profile in all 3 phases
for i = 1:params.Nx
    x = params.dx/2 + (i-1)*params.dx;    
    if (x < s_analytical_LV) 
         T(i,1) = params.T_inf + (params.T_v - params.T_inf) * (erf(x / (2 * sqrt(params.alpha_v * time))) / erf(lambda_sol * sqrt(diffusivity_ratio_LV)));

    elseif (x >= s_analytical_LV && x < s_analytical_LS)
         A_2 = (params.T_m-params.T_v)/(erf(beta_sol*sqrt(params.alpha_s/params.alpha_l)-lambda_sol*(1-params.R_rho_VL))-erf(lambda_sol*params.R_rho_VL));
         B_2 = params.T_v-A_2*erf(lambda_sol*params.R_rho_VL);
         T(i,1) = B_2 + A_2 * erf(x / (2 * sqrt(params.alpha_l * time))-lambda_sol*(1-params.R_rho_VL)); 

    elseif (x > s_analytical_LS)
         A_1 = (params.T_m-params.T_0)/(erfc(beta_sol*params.R_rho_LS - lambda_sol*sqrt(diffusivity_ratio_LS)*(params.R_rho_LS-params.R_rho_VS)));
         T(i,1) = params.T_0 +  A_1 * erfc(x / (2 * sqrt(params.alpha_s * time))-lambda_sol*sqrt(diffusivity_ratio_LS)*(params.R_rho_LS-params.R_rho_VS)-beta_sol*(1-params.R_rho_LS));

    elseif (x == s_analytical_LS)
         T(i,1) = params.T_m;
    else 
         T(i,1) = params.T_v;
    end
end