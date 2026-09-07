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
