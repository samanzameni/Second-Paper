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

    N_particles = 400;

    d_max = area_size * sqrt(2);

    margin_dB = 1;


    %% ========================================================
    % LOOP OVER EPSILON
    %% ========================================================


    for ei = 1:n_eps


        eps = epsilon_vals_in(ei);


        if isinf(eps)

            b = 0;

        else

            b = Delta_f_in / eps;

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
            % PARTICLE FILTER ATTACK
            %
            % Observation:
            %
            % P_DP = P_true + Gaussian + Laplace
            %
            % The likelihood is modeled using the
            % Normal-Laplace distribution.
            %% ===============================================


            d_est = zeros(N_SU,1);



            for i = 1:N_SU



                %% ============================================
                % 1. PARTICLE INITIALIZATION
                %
                % Uniform prior over possible distance
                %% ============================================


                d_particles = ...
                rand(N_particles,1) * d_max;



                %% ============================================
                % 2. PREDICTED TRUE POWER
                %
                % P_pred = P_max exp(-d/gamma)
                %% ============================================


                p_pred = ...
                P_max_in * ...
                exp( ...
                -d_particles / gamma_val);



                %% ============================================
                % OBSERVATION ERROR
                %% ============================================


                z = ...
                MTP_DP(i) - p_pred;



                %% ============================================
                % 3. LIKELIHOOD
                %% ============================================


                if b == 0


                    % ========================================
                    % ONLY GAUSSIAN MEASUREMENT NOISE
                    %
                    % P_DP = P_obs
                    % ========================================


                    sigma_eff = max(sigma_m,1e-6);


                    weights = ...
                    exp( ...
                    -0.5 * ...
                    (z ./ sigma_eff).^2);



                else


                    % ========================================
                    % NORMAL-LAPLACE LIKELIHOOD
                    %
                    % Z = Gaussian(0,sigma_m^2)
                    %     + Laplace(0,b)
                    %
                    % Closed-form PDF
                    % ========================================


                    sigma_eff = ...
                    max(sigma_m,1e-6);



                    term1 = ...
                    exp( ...
                    -z ./ b);


                    term2 = ...
                    exp( ...
                    z ./ b);



                    Phi1 = ...
                    0.5 * ...
                    ( ...
                    1 + erf( ...
                    (z./sigma_eff - ...
                    sigma_eff/b) / sqrt(2)));


                    Phi2 = ...
                    0.5 * ...
                    ( ...
                    1 + erf( ...
                    (-z./sigma_eff - ...
                    sigma_eff/b) / sqrt(2)));



                    weights = ...
                    (1/(2*b)) * ...
                    exp( ...
                    sigma_eff^2 / ...
                    (2*b^2)) .* ...
                    ( ...
                    term1 .* Phi1 + ...
                    term2 .* Phi2);



                    % ========================================
                    % Numerical stability
                    % ========================================


                    weights(~isfinite(weights)) = 0;

                    weights(weights < 0) = 0;


                end



                %% ============================================
                % 4. NORMALIZE PARTICLE WEIGHTS
                %% ============================================


                sum_w = sum(weights);


                if sum_w < 1e-300 || ...
                   any(isnan(weights)) || ...
                   any(isinf(weights))


                    weights = ...
                    ones(N_particles,1) / N_particles;


                else


                    weights = weights / sum_w;


                end



                %% ============================================
                % 5. ESTIMATED DISTANCE
                %
                % Weighted posterior mean
                %% ============================================


                d_est(i) = ...
                sum( ...
                weights .* d_particles);


            end



            %% ===============================================
            % ATTACK RMSE
            %% ===============================================


            rmse_runs(r) = ...
            sqrt( ...
            mean( ...
            (d_est - dists_true_min).^2));


        end



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
