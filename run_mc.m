%% ============================================================
%  HELPER: Single Monte Carlo run for given epsilon, Pmax, gamma_val
%% ============================================================
function [mae_val, eta_val, rmse_attack] = run_mc(epsilon_vals_in, Delta_f_in, N_SU, N_PU, P_max_in, area_size, num_runs, gamma_val)
    n_eps = length(epsilon_vals_in);
    mae_val      = zeros(1, n_eps);
    eta_val      = zeros(1, n_eps);
    rmse_attack  = zeros(1, n_eps);

    for ei = 1:n_eps
        eps = epsilon_vals_in(ei);
        b   = Delta_f_in / eps;

        mae_runs   = zeros(num_runs, 1);
        eta_runs   = zeros(num_runs, 1);
        rmse_runs  = zeros(num_runs, 1);

        for r = 1:num_runs
            su_locs = rand(N_SU, 2) * area_size;
            pu_locs = rand(N_PU, 2) * area_size;

            MTP_true = zeros(N_SU, 1);
            for i = 1:N_SU
                dists = sqrt(sum((su_locs(i,:) - pu_locs).^2, 2));
                d = min(dists);
                MTP_true(i) = P_max_in * exp(-d / gamma_val);
            end

            % Laplace noise generation (safeguarded against b=0 / eps=inf)
            if b == 0
                noise = zeros(N_SU, 1);
            else
                U = rand(N_SU, 1) - 0.5;
                % Prevent exact -0.5 or 0.5 to avoid log(0) = -inf
                U(U == 0.5) = 0.499999; 
                U(U == -0.5) = -0.499999;
                noise = -b * sign(U) .* log(1 - 2*abs(U));
            end
            
            MTP_DP = max(0, min(P_max_in, MTP_true + noise));

            mae_runs(r)  = mean(abs(MTP_DP - MTP_true));
            eta_runs(r)  = sum(MTP_DP > 0) / N_SU;

            % Adversary RMSE: adversary tries to infer MTP_true from MTP_DP
            rmse_runs(r) = sqrt(mean((MTP_DP - MTP_true).^2));
        end

        mae_val(ei)     = mean(mae_runs);
        eta_val(ei)     = mean(eta_runs);
        rmse_attack(ei) = mean(rmse_runs);
    end
end