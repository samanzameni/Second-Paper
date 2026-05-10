%% ============================================================
%  HELPER: Baseline mechanism simulations
%% ============================================================
function [rmse_lo_m, rmse_kanon_m, rmse_cr_m, rmse_base_m, mae_lo, mae_kanon, mae_cr, eta_lo, eta_kanon, eta_cr] = ...
        run_baselines(N_SU, N_PU, P_max, area_size, num_runs, gamma_val, Delta_f)

    rmse_lo_r    = zeros(num_runs,1);
    rmse_kanon_r = zeros(num_runs,1);
    rmse_cr_r    = zeros(num_runs,1);
    rmse_base_r  = zeros(num_runs,1);
    mae_lo_r     = zeros(num_runs,1);
    mae_kanon_r  = zeros(num_runs,1);
    mae_cr_r     = zeros(num_runs,1);
    margin_dB      = 0.5;
    eta_lo_r     = zeros(num_runs,1);
    eta_kanon_r  = zeros(num_runs,1);
    eta_cr_r     = zeros(num_runs,1);

    sigma_lo  = 20;   
    k_anon    = 3;    
    cr_size   = 100;  
    
    P_th = P_max * exp(-area_size / gamma_val); % آستانه حساسیت توان
    eps_power = 1e-15; % برای جلوگیری از لگاریتم صفر

    for r = 1:num_runs
        su_locs = rand(N_SU, 2) * area_size;
        pu_locs = rand(N_PU, 2) * area_size;

        MTP_true = zeros(N_SU, 1);
        dists_true_min = zeros(N_SU, 1);
        for i = 1:N_SU
            dists = sqrt(sum((su_locs(i,:) - pu_locs).^2, 2));
            dists_true_min(i) = min(dists);
            MTP_true(i) = P_max * exp(-dists_true_min(i) / gamma_val);
        end
        
        connected_true = MTP_true >= P_th;
        num_connected_true = max(1, sum(connected_true));

        % 1. No protection (Base)
        gps_error_m = 3;
        su_locs_base = su_locs + (gps_error_m / sqrt(2)) * randn(N_SU, 2);
        MTP_base = zeros(N_SU, 1);
        for i = 1:N_SU
            dists_base = sqrt(sum((su_locs_base(i,:) - pu_locs).^2, 2));
            MTP_base(i) = P_max * exp(-min(dists_base) / gamma_val);
        end
        MTP_base_clipped = max(0, min(P_max, MTP_base));
        d_est_base = -gamma_val * log(max(MTP_base_clipped, eps_power) / P_max);
        rmse_base_r(r) = sqrt(mean((d_est_base - dists_true_min).^2));

        % 2. LO
        su_locs_lo = su_locs + sigma_lo * randn(N_SU, 2);
        MTP_lo = zeros(N_SU, 1);
        for i = 1:N_SU
            dists = sqrt(sum((su_locs_lo(i,:) - pu_locs).^2, 2));
            MTP_lo(i) = P_max * exp(-min(dists) / gamma_val);
        end
        MTP_lo_clipped = max(0, min(P_max, MTP_lo));
        d_est_lo = -gamma_val * log(max(MTP_lo_clipped, eps_power) / P_max);
        rmse_lo_r(r) = sqrt(mean((d_est_lo - dists_true_min).^2));
        mae_lo_r(r)  = mean(abs(MTP_lo_clipped - MTP_true));
        eta_lo_r(r) = sum(MTP_lo_clipped >= (MTP_true - margin_dB)) / length(MTP_true);
        % 3. K-Anonymity
        dists_su = pdist2(su_locs, su_locs);
        dists_su(logical(eye(N_SU))) = inf;
        MTP_kanon = zeros(N_SU, 1);
        for i = 1:N_SU
            [~, idx_k] = sort(dists_su(i,:));
            neighbors = idx_k(1:k_anon-1);
            MTP_kanon(i) = mean(MTP_true([i; neighbors']));
        end
        MTP_kanon_clipped = max(0, min(P_max, MTP_kanon));
        d_est_kanon = -gamma_val * log(max(MTP_kanon_clipped, eps_power) / P_max);
        rmse_kanon_r(r) = sqrt(mean((d_est_kanon - dists_true_min).^2));
        mae_kanon_r(r)  = mean(abs(MTP_kanon_clipped - MTP_true));
        eta_kanon_r(r) = sum(MTP_kanon_clipped >= (MTP_true - margin_dB)) / length(MTP_true);
        % 4. Cloaking Regions
        cell_x = floor(su_locs(:,1) / cr_size);
        cell_y = floor(su_locs(:,2) / cr_size);
        MTP_cr = zeros(N_SU, 1);
        for i = 1:N_SU
            same_cell = (cell_x == cell_x(i)) & (cell_y == cell_y(i));
            MTP_cr(i) = mean(MTP_true(same_cell));
        end
        MTP_cr_clipped = max(0, min(P_max, MTP_cr));
        d_est_cr = -gamma_val * log(max(MTP_cr_clipped, eps_power) / P_max);
        rmse_cr_r(r) = sqrt(mean((d_est_cr - dists_true_min).^2));
        mae_cr_r(r)  = mean(abs(MTP_cr_clipped - MTP_true));
        eta_cr_r(r) = sum(MTP_cr_clipped >= (MTP_true - margin_dB)) / length(MTP_true);
    end

    rmse_lo_m    = mean(rmse_lo_r);
    rmse_kanon_m = mean(rmse_kanon_r);
    rmse_cr_m    = mean(rmse_cr_r);
    rmse_base_m  = mean(rmse_base_r);
    mae_lo       = mean(mae_lo_r);
    mae_kanon    = mean(mae_kanon_r);
    mae_cr       = mean(mae_cr_r);
    eta_lo       = min(1, mean(eta_lo_r));
    eta_kanon    = min(1, mean(eta_kanon_r));
    eta_cr       = min(1, mean(eta_cr_r));
end