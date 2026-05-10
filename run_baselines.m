%% ============================================================
%  HELPER: Baseline mechanism simulations
%% ============================================================
function [rmse_lo, rmse_kanon, rmse_cr, rmse_base, mae_lo, mae_kanon, mae_cr, eta_lo, eta_kanon, eta_cr] = ...
        run_baselines(N_SU, N_PU, P_max, area_size, num_runs, gamma_val, Delta_f)

    rmse_lo_r    = zeros(num_runs,1);
    rmse_kanon_r = zeros(num_runs,1);
    rmse_cr_r    = zeros(num_runs,1);
    rmse_base_r  = zeros(num_runs,1);
    mae_lo_r     = zeros(num_runs,1);
    mae_kanon_r  = zeros(num_runs,1);
    mae_cr_r     = zeros(num_runs,1);

    eta_lo_r     = zeros(num_runs,1);
    eta_kanon_r  = zeros(num_runs,1);
    eta_cr_r     = zeros(num_runs,1);

    sigma_lo  = 20;   % LO sigma=20m
    k_anon    = 3;    % K-Anonymity k=3
    cr_size   = 100;  % Cloaking region 100m

    for r = 1:num_runs
        su_locs = rand(N_SU, 2) * area_size;
        pu_locs = rand(N_PU, 2) * area_size;

        MTP_true = zeros(N_SU, 1);
        for i = 1:N_SU
            dists = sqrt(sum((su_locs(i,:) - pu_locs).^2, 2));
            MTP_true(i) = P_max * exp(-min(dists) / gamma_val);
        end

        % 1. No protection (Realistic Baseline with ~3m GPS error)
                gps_error_m = 3;
                su_locs_base = su_locs + (gps_error_m / sqrt(2)) * randn(N_SU, 2);
                MTP_base = zeros(N_SU, 1);
                for i = 1:N_SU
                    dists_base = sqrt(sum((su_locs_base(i,:) - pu_locs).^2, 2));
                    MTP_base(i) = P_max * exp(-min(dists_base) / gamma_val);
                end
                MTP_base_clipped = max(0, min(P_max, MTP_base));
                rmse_base_r(r) = sqrt(mean((MTP_base_clipped - MTP_true).^2));

        % LO
        su_locs_lo = su_locs + sigma_lo * randn(N_SU, 2);
        MTP_lo = zeros(N_SU, 1);
        for i = 1:N_SU
            dists = sqrt(sum((su_locs_lo(i,:) - pu_locs).^2, 2));
            MTP_lo(i) = P_max * exp(-min(dists) / gamma_val);
        end
        % شبیه‌سازی شرایط محدودیت توان (کلیپ کردن مقادیر غیرمجاز مشابه DP)
        MTP_lo_clipped = max(0, min(P_max, MTP_lo));
        rmse_lo_r(r) = sqrt(mean((MTP_lo_clipped - MTP_true).^2));
        mae_lo_r(r)  = mean(abs(MTP_lo_clipped - MTP_true));
        eta_lo_r(r)  = sum(MTP_lo_clipped > 0) / N_SU;

        % K-Anonymity
        dists_su = pdist2(su_locs, su_locs);
        dists_su(logical(eye(N_SU))) = inf;
        MTP_kanon = zeros(N_SU, 1);
        for i = 1:N_SU
            [~, idx_k] = sort(dists_su(i,:));
            neighbors = idx_k(1:k_anon-1);
            MTP_kanon(i) = mean(MTP_true([i; neighbors']));
        end
        MTP_kanon_clipped = max(0, min(P_max, MTP_kanon));
        rmse_kanon_r(r) = sqrt(mean((MTP_kanon_clipped - MTP_true).^2));
        mae_kanon_r(r)  = mean(abs(MTP_kanon_clipped - MTP_true));
        eta_kanon_r(r)  = sum(MTP_kanon_clipped > 0) / N_SU;

        % Cloaking Regions
        cell_x = floor(su_locs(:,1) / cr_size);
        cell_y = floor(su_locs(:,2) / cr_size);
        MTP_cr = zeros(N_SU, 1);
        for i = 1:N_SU
            same_cell = (cell_x == cell_x(i)) & (cell_y == cell_y(i));
            MTP_cr(i) = mean(MTP_true(same_cell));
        end
        MTP_cr_clipped = max(0, min(P_max, MTP_cr));
        rmse_cr_r(r) = sqrt(mean((MTP_cr_clipped - MTP_true).^2));
        mae_cr_r(r)  = mean(abs(MTP_cr_clipped - MTP_true));
        eta_cr_r(r)  = sum(MTP_cr_clipped > 0) / N_SU;
    end

    rmse_lo    = mean(rmse_lo_r);
    rmse_kanon = mean(rmse_kanon_r);
    rmse_cr    = mean(rmse_cr_r);
    rmse_base  = mean(rmse_base_r);
    mae_lo     = mean(mae_lo_r);
    mae_kanon  = mean(mae_kanon_r);
    mae_cr     = mean(mae_cr_r);
    eta_lo     = mean(eta_lo_r);
    eta_kanon  = mean(eta_kanon_r);
    eta_cr     = mean(eta_cr_r);
end