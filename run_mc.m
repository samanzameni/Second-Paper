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
d_est = zeros(N_SU,1);

alpha_mix   = 0.85;          % وزن پیشنهاد تحلیلی

for i = 1:N_SU
    P_obs = MTP_DP(i);
    
    % -------------------------------------------------
    % 1. پیشنهاد تحلیلی (معکوس مدل مسیر)
    % -------------------------------------------------
    P_clip = max(min(P_obs, P_max_in*0.999), P_max_in*1e-6);
    d_anal = -gamma_val * log(P_clip / P_max_in);
    d_anal = max(0, min(d_max, d_anal));
    
    % انحراف معیار تقریبی فاصله ناشی از نویز اندازه‌گیری
    sigma_d = max(sigma_m,1e-6) * gamma_val / max(P_clip, 1e-3);
    sigma_d = max(sigma_d, 2);   % کف برای جلوگیری از فروپاشی
    
    % -------------------------------------------------
    % 2. تولید ذرات با ترکیب پیشنهاد تحلیلی + یکنواخت
    % -------------------------------------------------
    n_anal = round(alpha_mix * N_particles);
    n_unif = N_particles - n_anal;
    
    d_particles = [ ...
        d_anal + sigma_d * randn(n_anal,1);   % پیشنهاد تحلیلی
        rand(n_unif,1) * d_max ];             % پیشنهاد یکنواخت (ایمنی)
    
    d_particles = max(0, min(d_max, d_particles));
    
    % -------------------------------------------------
    % 3. قدرت پیش‌بینی‌شده و باقی‌مانده
    % -------------------------------------------------
    p_pred = P_max_in * exp(-d_particles / gamma_val);
    z      = P_obs - p_pred;
    
    % -------------------------------------------------
    % 4. وزن‌های likelihood
    % -------------------------------------------------
    if b == 0
        % فقط نویز گاوسی
        sigma_eff = max(sigma_m, 1e-6);
        weights = exp(-0.5 * (z ./ sigma_eff).^2);
    else
        % چگالی Normal-Laplace بسته
        sigma_eff = max(sigma_m, 1e-6);
        term1 = exp(-z ./ b);
        term2 = exp( z ./ b);
        Phi1  = 0.5 * (1 + erf( ( z./sigma_eff - sigma_eff/b ) / sqrt(2) ));
        Phi2  = 0.5 * (1 + erf( (-z./sigma_eff - sigma_eff/b ) / sqrt(2) ));
        weights = (1/(2*b)) * exp(sigma_eff^2 / (2*b^2)) .* ...
                  (term1 .* Phi1 + term2 .* Phi2);
        weights(~isfinite(weights) | weights < 0) = 0;
    end
    
    % -------------------------------------------------
    % 5. نرمال‌سازی + Residual Resampling
    % -------------------------------------------------
    sum_w = sum(weights);
    if sum_w < 1e-300 || any(~isfinite(weights))
        weights = ones(N_particles,1) / N_particles;
    else
        weights = weights / sum_w;
    end
    
    % Effective Sample Size
    ess = 1 / sum(weights.^2);
    if ess < 0.4 * N_particles
        % Residual resampling
        N_res = floor(N_particles * weights);
        residual = weights - N_res / N_particles;
        residual = residual / sum(residual + eps);
        
        idx = [];
        for k = 1:N_particles
            if N_res(k) > 0
                idx = [idx; repmat(k, N_res(k), 1)];
            end
        end
        n_remain = N_particles - length(idx);
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
        weights = ones(N_particles,1) / N_particles;
    end
    
    % -------------------------------------------------
    % 6. تخمین نهایی (میانگین وزن‌دار) + تصحیح بایاس
    % -------------------------------------------------
    d_mean = sum(weights .* d_particles);
    
    if b == 0
        % تصحیح بایاس مرتبه اول برای تبدیل لگاریتمی
        bias = (sigma_m^2 * gamma_val) / (2 * max(P_clip,1e-3)^2);
        d_mean = d_mean - bias;
    end
    
    d_est(i) = max(0, min(d_max, d_mean));
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
