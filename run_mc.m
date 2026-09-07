%% ============================================================
% =============================================================
%
%  FUNCTION 1:
%  MAIN MONTE CARLO SIMULATION
%
% =============================================================
% ============================================================


function [mae_val,eta_val,rmse_attack_m] = ...
run_mc( ...
epsilon_vals_in, ...
Delta_f_in, ...
N_SU, ...
N_PU, ...
P_max_in, ...
area_size, ...
num_runs, ...
gamma_val, ...
sigma_m)



    n_eps = length(epsilon_vals_in);


    mae_val = zeros(1,n_eps);

    eta_val = zeros(1,n_eps);

    rmse_attack_m = zeros(1,n_eps);



    %% ========================================================
    % PARTICLE FILTER PARAMETERS
    %% ========================================================

    N_particles = 1000;

    d_max = area_size * sqrt(2);

    margin_dB = 1;


    %% ========================================================
    % LOOP OVER EPSILON
    %% ========================================================


    for ei = 1:n_eps


        eps_val = epsilon_vals_in(ei);


        if isinf(eps_val)

            b = 0;

        else

            b = Delta_f_in / eps_val;

        end



        mae_runs = zeros(num_runs,1);

        eta_runs = zeros(num_runs,1);

        rmse_runs = zeros(num_runs,1);



        %% ====================================================
        % MONTE CARLO LOOP
        %% ====================================================


        for r = 1:num_runs



            %% =================================================
            % GENERATE SU AND PU LOCATIONS
            %% =================================================


            su_locs = ...
            rand(N_SU,2) * area_size;


            pu_locs = ...
            rand(N_PU,2) * area_size;



            MTP_true = zeros(N_SU,1);

            MTP_obs = zeros(N_SU,1);

            MTP_DP = zeros(N_SU,1);

            dists_true_min = zeros(N_SU,1);



            %% ===============================================
            % CALCULATE TRUE DISTANCE AND TRUE POWER
            %% ===============================================


            for i = 1:N_SU


                dists = sqrt( ...
                sum( ...
                (su_locs(i,:) - pu_locs).^2, ...
                2));


                dists_true_min(i) = min(dists);



                % =============================================
                % PATH-LOSS / PROPAGATION MODEL
                %
                % P_true = P_max * exp(-d/gamma)
                % =============================================


                MTP_true(i) = ...
                P_max_in * ...
                exp( ...
                -dists_true_min(i) / gamma_val);


            end



            %% ===============================================
            % MEASUREMENT MODEL
            %
            % P_obs = clip(P_true + N(0,sigma_m^2))
            %% ===============================================


            measurement_noise = ...
            sigma_m * randn(N_SU,1);


            MTP_obs_unclipped = ...
            MTP_true + measurement_noise;


            MTP_obs = max( ...
            0, ...
            min( ...
            P_max_in, ...
            MTP_obs_unclipped));



            %% ===============================================
            % DIFFERENTIAL PRIVACY NOISE
            %
            % P_DP = clip(P_obs + Lap(0,b))
            %
            % b = Delta_f / epsilon
            %% ===============================================


            if b == 0


                dp_noise = zeros(N_SU,1);


            else


                U = rand(N_SU,1) - 0.5;


                U(U >= 0.5) = 0.499999;

                U(U <= -0.5) = -0.499999;


                dp_noise = ...
                -b .* ...
                sign(U) .* ...
                log( ...
                1 - 2*abs(U));


            end



            MTP_DP_unclipped = ...
            MTP_obs + dp_noise;



            %% ===============================================
            % MAE OF DP MECHANISM
            %
            % IMPORTANT:
            %
            % MAE = |P_DP - P_obs|
            %
            % Therefore:
            %
            % E[|Lap(0,b)|] = b = Delta_f/epsilon
            %% ===============================================


            mae_runs(r) = ...
            mean( ...
            abs( ...
            MTP_DP_unclipped - MTP_obs));


            %% ===============================================
            % CLIP FINAL DP POWER
            %% ===============================================


            MTP_DP = max( ...
            0, ...
            min( ...
            P_max_in, ...
            MTP_DP_unclipped));



            %% ===============================================
            % SPECTRUM UTILIZATION EFFICIENCY
            %
            % Compare protected power with true power
            %% ===============================================


            connected_DP = ...
            MTP_DP >= ...
            (MTP_true - margin_dB);


            eta_runs(r) = ...
            sum(connected_DP) / N_SU;



     %% ===============================================
% PARTICLE FILTER ATTACK (Improved)
%
% - Analytic proposal from inverse path-loss
% - Mixture proposal (analytic + uniform)
% - Residual resampling
% - Bias correction for pure Gaussian case
%% ===============================================
d_est = particle_attack(MTP_DP, P_max_in, gamma_val, sigma_m, b, N_particles, d_max);


            %% ===============================================
            % ATTACK RMSE
            %% ===============================================


            rmse_runs(r) = ...
            sqrt( ...
            mean( ...
            (d_est - dists_true_min).^2));





        %% ====================================================
        % AVERAGE OVER MONTE CARLO RUNS
        %% ====================================================


        mae_val(ei) = ...
        mean(mae_runs);


        eta_val(ei) = ...
        min(1,mean(eta_runs));


        rmse_attack_m(ei) = ...
        mean(rmse_runs);


    end


end
