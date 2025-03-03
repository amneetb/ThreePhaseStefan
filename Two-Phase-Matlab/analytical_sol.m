function [T, dsdt_analytical, s_analytical, lambda_sol] = analytical_sol(time,params)


T = zeros(params.Nx,1);

dT_right = params.T_melt - params.T_right;
dT_left = params.T_melt - params.T_left;
diffusivity_ratio = params.alpha_s/params.alpha_l;

    
% Left-hand side calculation
func = @(lambda) lambda * sqrt(params.alpha_s)*params.rho_l * ((params.cp_l-params.cp_s)*(params.T_melt-params.T_r)+params.L_m) +...
params.k_s*dT_right *exp(-(lambda / params.Rp)^2) / (erfc(lambda / params.Rp)*sqrt(pi * params.alpha_s)) +...
params.k_l*dT_left*exp(-lambda^2 * diffusivity_ratio)/ (erf(lambda * sqrt(diffusivity_ratio))*sqrt(pi * params.alpha_l));
    
initialGuess = 2.0;
lambda_sol = fzero(func,initialGuess);
   
dsdt_analytical = lambda_sol*sqrt(params.alpha_s/time);
s_analytical = 2*lambda_sol*(sqrt(time*params.alpha_s));

% Temperature profile in solid phase
for i = 1:params.Nx
    x = params.dx/2 + (i-1)*params.dx;    
    if (x < s_analytical) 
         T(i,1) = params.T_left + (params.T_melt - params.T_left)* ...
             (erf(x/(2 * sqrt(params.alpha_l * time))) / erf(lambda_sol * sqrt(diffusivity_ratio)));
    elseif (x > s_analytical)
         T(i,1) = params.T_right + (params.T_melt - params.T_right)* ...
            (erfc(x/(2*sqrt(params.alpha_s*time)) - ...
            lambda_sol*(1 - 1 / params.Rp))/erfc(lambda_sol/params.Rp));     
    else
        T(i,1) = params.T_melt;
    end
end
