%% ============================================================
%  HELPER: Single Monte Carlo run for given epsilon, Pmax, gamma_val
%% ============================================================
function [mae_val, eta_val, rmse_attack_m] = run_mc(epsilon_vals_in, Delta_f_in, N_SU, N_PU, P_max_in, area_size, num_runs, gamma_val)
    n_eps = length(epsilon_vals_in);
    mae_val        = zeros(1, n_eps);
    eta_val        = zeros(1, n_eps);
    rmse_attack_m  = zeros(1, n_eps);
    
    P_th = P_max_in * exp(-area_size / gamma_val); % آستانه حساسیت گیرنده
    eps_power = 1e-15;

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
            dists_true_min = zeros(N_SU, 1);
            for i = 1:N_SU
                dists = sqrt(sum((su_locs(i,:) - pu_locs).^2, 2));
                dists_true_min(i) = min(dists);
                MTP_true(i) = P_max_in * exp(-dists_true_min(i) / gamma_val);
            end

            if b == 0
                noise = zeros(N_SU, 1);
            else
                U = rand(N_SU, 1) - 0.5;
                U(U == 0.5) = 0.499999; 
                U(U == -0.5) = -0.499999;
                noise = -b * sign(U) .* log(1 - 2*abs(U));
            end
            
            % افزودن نویز قبل از محدودسازی (برای محاسبه صحیح MAE)
            MTP_noisy_unclipped = MTP_true + noise;
            mae_runs(r) = mean(abs(MTP_noisy_unclipped - MTP_true));
            
            % اعمال محدودیت توان 
            MTP_DP = max(0, min(P_max_in, MTP_noisy_unclipped));

            % محاسبه راندمان (eta) بر اساس آستانه منطقی
            connected_true = MTP_true >= P_th;
            connected_DP   = MTP_DP >= P_th;
            if sum(connected_true) > 0
                eta_runs(r) = sum(connected_DP & connected_true) / sum(connected_true);            else
                eta_runs(r) = 1;
            end

            % محاسبه خطای مکانی مهاجم (RMSE مکانی به متر)
            d_est = -gamma_val * log(max(MTP_DP, eps_power) / P_max_in);
            rmse_runs(r) = sqrt(mean((d_est - dists_true_min).^2));
        end

        mae_val(ei)       = mean(mae_runs);
        eta_val(ei)       = min(1, mean(eta_runs));
        rmse_attack_m(ei) = mean(rmse_runs);
    end
end