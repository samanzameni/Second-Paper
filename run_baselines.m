%% ============================================================
% =============================================================
%
% FUNCTION 2:
% BASELINE MECHANISM SIMULATIONS
%
% UPDATED WITH MEASUREMENT NOISE
%
% =============================================================
% ============================================================


function [ ...
rmse_lo_m, ...
rmse_kanon_m, ...
rmse_cr_m, ...
rmse_base_m, ...
mae_lo, ...
mae_kanon, ...
mae_cr, ...
eta_lo, ...
eta_kanon, ...
eta_cr] = ...
run_baselines( ...
N_SU, ...
N_PU, ...
P_max, ...
area_size, ...
num_runs, ...
gamma_val, ...
Delta_f, ...
sigma_m)



    %% ========================================================
    % INITIALIZATION
    %% ========================================================


    rmse_lo_r = zeros(num_runs,1);

    rmse_kanon_r = zeros(num_runs,1);

    rmse_cr_r = zeros(num_runs,1);

    rmse_base_r = zeros(num_runs,1);


    mae_lo_r = zeros(num_runs,1);

    mae_kanon_r = zeros(num_runs,1);

    mae_cr_r = zeros(num_runs,1);


    eta_lo_r = zeros(num_runs,1);

    eta_kanon_r = zeros(num_runs,1);

    eta_cr_r = zeros(num_runs,1);



    margin_dB = 1;
    sigma_lo = 20;
    k_anon = 3;
    cr_size = 100;
    N_particles = 1000;                 
    d_max = area_size * sqrt(2);

    %% ========================================================
    % MONTE CARLO LOOP
    %% ========================================================


    for r = 1:num_runs



        %% ====================================================
        % USER LOCATIONS
        %% ====================================================


        su_locs = ...
        rand(N_SU,2) * area_size;


        pu_locs = ...
        rand(N_PU,2) * area_size;



        MTP_true = zeros(N_SU,1);

        dists_true_min = zeros(N_SU,1);



        %% ====================================================
        % TRUE POWER
        %
        % P_true = P_max exp(-d/gamma)
        %% ====================================================


        for i = 1:N_SU


            dists = sqrt( ...
            sum( ...
            (su_locs(i,:) - pu_locs).^2, ...
            2));


            dists_true_min(i) = ...
            min(dists);


            MTP_true(i) = ...
            P_max * ...
            exp( ...
            -dists_true_min(i) / gamma_val);


        end



        %% ====================================================
        % OBSERVED POWER
        %
        % P_obs = clip(P_true + Gaussian noise)
        %% ====================================================


        measurement_noise = ...
        sigma_m * randn(N_SU,1);


        MTP_obs = ...
        MTP_true + measurement_noise;


        MTP_obs = max( ...
        0, ...
        min( ...
        P_max, ...
        MTP_obs));



        %% ====================================================
        % 1. NO PROTECTION BASELINE
        %
        % Adversary receives observed power
        %% ====================================================


       d_est_base = particle_attack(MTP_obs, P_max, gamma_val, sigma_m, 0, N_particles, d_max);



        rmse_base_r(r) = ...
        sqrt( ...
        mean( ...
        (d_est_base - dists_true_min).^2));



        %% ====================================================
        % 2. LOCATION OBFUSCATION
        %% ====================================================


        su_locs_lo = ...
        su_locs + ...
        sigma_lo * randn(N_SU,2);


        MTP_lo = zeros(N_SU,1);



        for i = 1:N_SU


            dists = sqrt( ...
            sum( ...
            (su_locs_lo(i,:) - pu_locs).^2, ...
            2));


            MTP_lo(i) = ...
            P_max * ...
            exp( ...
            -min(dists) / gamma_val);


        end



        % Add measurement noise


        MTP_lo = ...
        MTP_lo + ...
        sigma_m * randn(N_SU,1);


        MTP_lo_clipped = ...
        max( ...
        0, ...
        min( ...
        P_max, ...
        MTP_lo));



        d_est_lo = particle_attack(MTP_lo_clipped, P_max, gamma_val, sigma_m, 0, N_particles, d_max);



        rmse_lo_r(r) = ...
        sqrt( ...
        mean( ...
        (d_est_lo - dists_true_min).^2));


        mae_lo_r(r) = ...
        mean( ...
        abs( ...
        MTP_lo_clipped - MTP_obs));


        eta_lo_r(r) = ...
        sum( ...
        MTP_lo_clipped >= ...
        (MTP_true - margin_dB)) / N_SU;



      %% ====================================================
% 3. K-ANONYMITY 
%% ====================================================

dists_su = pdist2(su_locs, su_locs);
dists_su(logical(eye(N_SU))) = inf;

MTP_kanon = zeros(N_SU,1);

for i = 1:N_SU

    [~,idx_k] = sort(dists_su(i,:));
    neighbors = idx_k(1:k_anon-1);
    group_idx = [i, neighbors];

    group_centroid = mean(su_locs(group_idx,:), 1);

    dists_centroid = sqrt(sum((group_centroid - pu_locs).^2, 2));
    MTP_true_centroid = P_max * exp(-min(dists_centroid) / gamma_val);

    MTP_kanon(i) = MTP_true_centroid + sigma_m * randn();

end

MTP_kanon_clipped = max(0, min(P_max, MTP_kanon));

d_est_kanon = particle_attack(MTP_kanon_clipped, P_max, gamma_val, sigma_m, 0, N_particles, d_max);

rmse_kanon_r(r) = sqrt(mean((d_est_kanon - dists_true_min).^2));

mae_kanon_r(r) = mean(abs(MTP_kanon_clipped - MTP_obs));

eta_kanon_r(r) = sum(MTP_kanon_clipped >= (MTP_true - margin_dB)) / N_SU;



       %% ====================================================
% 4. CLOAKING REGIONS (Location Generalization — Fixed)
%% ====================================================

cell_x = floor(su_locs(:,1) / cr_size);
cell_y = floor(su_locs(:,2) / cr_size);

MTP_cr = zeros(N_SU,1);

for i = 1:N_SU

    same_cell = (cell_x == cell_x(i)) & (cell_y == cell_y(i));

    cell_centroid = mean(su_locs(same_cell,:), 1);

    dists_centroid = sqrt(sum((cell_centroid - pu_locs).^2, 2));
    MTP_true_centroid = P_max * exp(-min(dists_centroid) / gamma_val);

    MTP_cr(i) = MTP_true_centroid + sigma_m * randn();

end

MTP_cr_clipped = max(0, min(P_max, MTP_cr));

d_est_cr = particle_attack(MTP_cr_clipped, P_max, gamma_val, sigma_m, 0, N_particles, d_max);

rmse_cr_r(r) = sqrt(mean((d_est_cr - dists_true_min).^2));

mae_cr_r(r) = mean(abs(MTP_cr_clipped - MTP_obs));

eta_cr_r(r) = sum(MTP_cr_clipped >= (MTP_true - margin_dB)) / N_SU;

    end


    %% ========================================================
    % FINAL AVERAGES
    %% ========================================================


    rmse_lo_m = ...
    mean(rmse_lo_r);


    rmse_kanon_m = ...
    mean(rmse_kanon_r);


    rmse_cr_m = ...
    mean(rmse_cr_r);


    rmse_base_m = ...
    mean(rmse_base_r);



    mae_lo = ...
    mean(mae_lo_r);


    mae_kanon = ...
    mean(mae_kanon_r);


    mae_cr = ...
    mean(mae_cr_r);



    eta_lo = ...
    min(1,mean(eta_lo_r));


    eta_kanon = ...
    min(1,mean(eta_kanon_r));


    eta_cr = ...
    min(1,mean(eta_cr_r));


end

%% ============================================================
% =============================================================
%
%  FUNCTION 3:
%  PARTICLE FILTER ATTACK (used for baselines )
%
% =============================================================
% ============================================================
function d_est = particle_attack( ...
    P_dp_vec, ...
    P_max_in, ...
    gamma_val, ...
    sigma_m, ...
    b, ...
    N_particles, ...
    d_max)

% Improved particle / importance-sampling attack
% - Analytic proposal from the inverse path-loss model
% - Mixture of analytic + uniform proposals
% - Residual resampling
% - Bias correction for pure Gaussian noise

N = length(P_dp_vec);
d_est = zeros(N,1);

% Safety & defaults
sigma_eff = max(sigma_m, 1e-6);
N_part    = max(N_particles, 1000);   % حداقل ۱۲۰۰ ذره
alpha     = 0.85;                     % وزن پیشنهاد تحلیلی

for i = 1:N
    P_obs = P_dp_vec(i);
    
    % -------------------------------------------------
    % 1. پیشنهاد تحلیلی
    % -------------------------------------------------
    P_clip = max(min(P_obs, P_max_in*0.999), P_max_in*1e-6);
    d_anal = -gamma_val * log(P_clip / P_max_in);
    d_anal = max(0, min(d_max, d_anal));
    
    sigma_d = sigma_eff * gamma_val / max(P_clip, 1e-3);
    sigma_d = max(sigma_d, 2);
    
    % -------------------------------------------------
    % 2. تولید ذرات (مخلوط)
    % -------------------------------------------------
    n_anal = round(alpha * N_part);
    n_unif = N_part - n_anal;
    
    d_particles = [ ...
        d_anal + sigma_d * randn(n_anal,1);
        rand(n_unif,1) * d_max ];
    
    d_particles = max(0, min(d_max, d_particles));
    
    % -------------------------------------------------
    % 3. قدرت پیش‌بینی‌شده و باقی‌مانده
    % -------------------------------------------------
    p_pred = P_max_in * exp(-d_particles / gamma_val);
    z      = P_obs - p_pred;
    
    % -------------------------------------------------
    % 4. وزن likelihood
    % -------------------------------------------------
    if b == 0
        weights = exp(-0.5*(z./sigma_eff).^2);
    else
        term1 = exp(-z./b);
        term2 = exp( z./b);
        Phi1  = 0.5*(1 + erf( ( z./sigma_eff - sigma_eff/b)/sqrt(2) ));
        Phi2  = 0.5*(1 + erf( (-z./sigma_eff - sigma_eff/b)/sqrt(2) ));
        weights = (1/(2*b)) * exp(sigma_eff^2/(2*b^2)) .* ...
                  (term1.*Phi1 + term2.*Phi2);
        weights(~isfinite(weights) | weights < 0) = 0;
    end
    
    % -------------------------------------------------
    % 5. نرمال‌سازی + Residual Resampling
    % -------------------------------------------------
    sw = sum(weights);
    if sw < 1e-300 || any(~isfinite(weights))
        weights = ones(N_part,1)/N_part;
    else
        weights = weights / sw;
    end
    
    ess = 1 / sum(weights.^2);
    if ess < 0.4 * N_part
        N_res = floor(N_part * weights);
        residual = weights - N_res/N_part;
        residual = residual / sum(residual + eps);
        
        idx = [];
        for k = 1:N_part
            if N_res(k) > 0
                idx = [idx; repmat(k, N_res(k), 1)];
            end
        end
        n_remain = N_part - length(idx);
        if n_remain > 0
            cdf = cumsum(residual);
            u = rand(n_remain,1);
            extra = zeros(n_remain,1);
            for uu = 1:n_remain
                extra(uu) = find(cdf >= u(uu), 1, 'first');
            end
            idx = [idx; extra];
        end
        d_particles = d_particles(idx);
        weights = ones(N_part,1)/N_part;
    end
    
    % -------------------------------------------------
    % 6. تخمین نهایی + تصحیح بایاس
    % -------------------------------------------------
    d_mean = sum(weights .* d_particles);
    
    if b == 0
        bias = (sigma_eff^2 * gamma_val) / (2 * max(P_clip,1e-3)^2);
        d_mean = d_mean - bias;
    end
    
    d_est(i) = max(0, min(d_max, d_mean));
end
end