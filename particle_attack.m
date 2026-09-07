function d_est = particle_attack(P_dp_vec, P_max_in, gamma_val, sigma_m, b, N_particles, d_max)
% Improved particle / importance-sampling attack
% - Analytic proposal from the inverse path-loss model
% - Mixture of analytic + uniform proposals
% - Residual resampling
% - Bias correction for pure Gaussian noise

N = length(P_dp_vec);
d_est = zeros(N,1);

% Safety & defaults
sigma_eff = max(sigma_m, 1e-6);
N_part    = max(N_particles, 1000);   % حداقل ۱۰۰۰ ذره
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
        log_w = -0.5*(z./sigma_eff).^2;
    else
        log_Phi1 = log_normcdf( ( z./sigma_eff - sigma_eff/b) );
        log_Phi2 = log_normcdf( (-z./sigma_eff - sigma_eff/b) );
        log_A = -z./b + sigma_eff^2/(2*b^2) + log_Phi1;
        log_B =  z./b + sigma_eff^2/(2*b^2) + log_Phi2;
        m = max(log_A, log_B);
        m(~isfinite(m)) = -Inf;
        log_w = m + log(exp(log_A - m) + exp(log_B - m)) - log(2*b);
        log_w(~isfinite(log_w)) = -Inf;
    end

    log_w = log_w - max(log_w(isfinite(log_w)));
    weights = exp(log_w);
    weights(~isfinite(weights) | weights < 0) = 0;

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

% -------------------------------------------------
% تابع کمکی (داخل همان فایل)
% -------------------------------------------------
function y = log_normcdf(x)
% log(Phi(x)) با پایداری عددی برای x خیلی منفی
    y = zeros(size(x));
    idx_pos = x > -20;
    idx_neg = ~idx_pos;

    y(idx_pos) = log(0.5*(1+erf(x(idx_pos)/sqrt(2))));

    % برای x خیلی منفی، از تقریب مجانبی استفاده کن تا log(0) نشه
    xn = x(idx_neg);
    y(idx_neg) = -0.5*xn.^2 - log(-xn*sqrt(2*pi)) - log(1 - 1./xn.^2);
end