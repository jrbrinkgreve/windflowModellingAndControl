%implementation of Mohammadi et al.:
%https://www.mdpi.com/1996-1073/15/23/9135
%GIF Animation sweeping beta from 0 to 30 degrees
close all
clear all

function params = load_params(beta_deg)
    params.numxpoints = 100;
    params.numypoints = 100;
    params.numzpoints = 1;
    params.R = 63; %rotor radius
    params.D = params.R * 2;
    params.z_hub = 90; %hub height
    params.xlim_min = 1.*params.D;
    params.xlim_max = 10.*params.D;
    params.ylim_min = -1.5.*params.D;
    params.ylim_max = 1.5.*params.D;
    params.zlim_min = params.z_hub;
    params.zlim_max = params.z_hub;
    params.xrange = linspace(params.xlim_min, params.xlim_max, params.numxpoints);
    params.yrange = linspace(params.ylim_min, params.ylim_max, params.numypoints);
    params.zrange = linspace(params.zlim_min, params.zlim_max, params.numzpoints);
    params.CT = 0.66; %torque coeff.
    params.beta = deg2rad(beta_deg); %yaw angle (now parametric)
    params.u_hub = 8; %m./s, hub inflow vel.
    params.u_in = linspace(8,8, params.numzpoints);
    params.u_star = 0.4; %friction velocity
    params.lambda = 7.5; %tip-speed ratio
    params.k = 0.6 .* params.u_star ./ params.u_hub;
    params.alpha_gradient = 0;
    params.z_relative = params.zrange - params.z_hub;
    params.alpha = deg2rad(params.alpha_gradient * params.z_relative);
    params.angle_tolerance = 1e-10;
end

% Animation parameters
beta_values = linspace(0, 30, 31); % 0 to 30 degrees, 31 frames
gif_filename = 'wake_yaw_sweep.gif';
frame_delay = 0.1; % seconds between frames

% Create figure
fig = figure('Color','w', 'Position', [100 100 800 600]);

for frame_idx = 1:length(beta_values)
    beta_current = beta_values(frame_idx);

    % Load parameters with current beta
    params = load_params(beta_current);

    % Generate grid
    [X, Y] = ndgrid(params.xrange, params.yrange);
    X_vec = X(:); 
    Y_vec = Y(:);

    % Preallocate flowfield
    flowfield = zeros(params.numxpoints, params.numypoints, params.numzpoints);

    % Compute flowfield (z loop)
    for z_idx = 1:params.numzpoints
        z = params.zrange(z_idx);
        alpha = params.alpha(z_idx);
        gamma = params.beta + alpha;

        % eq 3
        A_star = (1 + sqrt(1- params.CT .* cos(gamma).^2) ) ./ (2.*sqrt(1- params.CT .* cos(gamma).^2));

        % eq 4
        xi_0_hat = params.R .* sqrt(A_star);

        % Build rotation matrix
        rotmtx = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];

        % Veered directions
        XYv_vec = [X_vec Y_vec] * rotmtx;
        Xv_vec = XYv_vec(:,1);
        Yv_vec = XYv_vec(:,2);

        % eq 5: determine dimless time 
        t_hat = -1.44 .* (params.u_in(z_idx) ./ params.u_star) .* (params.R ./ xi_0_hat) .* ...
                 (params.CT .* cos(gamma)^2 .* sin(gamma)) .* ...
                 (1 - exp(-0.35 .* (params.u_star ./ params.u_in(z_idx)) .* (Xv_vec ./ params.R)) );

        % eq 6: normalized wake center
        sgn_t_hat = sign(t_hat);
        abs_t_hat = abs(t_hat);
        num = (pi - 1) .* abs_t_hat.^3 + 2 .* sqrt(3).*pi^2 .* t_hat.^2 + 48.*(pi - 1)^2 .* abs_t_hat;
        den = 2.*pi.*(pi - 1).*t_hat.^2 + 4 .* sqrt(3) .*pi^2 .* abs_t_hat + 96.*(pi - 1)^2;
        y_hat_c = (num ./ den) .* sgn_t_hat - (2./pi) .* t_hat ./ (((z + params.z_hub)./xi_0_hat).^2 - 1);

        % eq 7: lateral wake center
        y_c = y_hat_c .* xi_0_hat;

        % eq 8: polar angle theta calc
        theta = atan((z - params.z_hub) ./ (Yv_vec - y_c));

        % eq 9: initial wake shape
        xi_0 = xi_0_hat .* abs(cos(gamma)) ./ sqrt(1 - sin(gamma).^2 .* sin(theta).^2 );

        if gamma < params.angle_tolerance
            chi = 0;
            xi_hat = ones(size(theta));
        else
            % eq 12: chi rotation rate calc
            chi = 1 ./ (params.lambda * sin(gamma));

            % eq 11: analytical params
            a = 1.263.*cos(0.33.*chi);
            c1 = 0.5 .* tanh(t_hat.^2 ./ (4.*a) );
            c2 = (-1/3) .* tanh( t_hat.^3 ./ (8.*a) );
            c3 = (-1/4) .* tanh( t_hat.^3 ./ (8.*a) );
            c4 = (-1/6) .* tanh( t_hat.^4 ./ (16.*a) );
            c5 = (5/16) .* tanh( t_hat.^4 ./ (16.*a) );
            c6 = (-5/48) .* tanh( t_hat.^4 ./ (16.*a) );
            c7 = (7/48) .* tanh( t_hat.^4 ./ (16.*a) );
            xi_hat = 1 - a*( c1.*cos(2.*theta)... 
                            + c2.*chi.*sin(2.*theta)...
                            + c3.*cos(3.*theta)...
                            + c4.*chi.^2.*cos(2.*theta)...
                            + c5.*chi.*sin(3.*theta)...
                            + c6.*cos(2.*theta)...
                            + c7.*cos(4.*theta) );
        end

        % eq 10: wake shape function
        xi = xi_0 .* xi_hat;

        % eq 13: wake width calculation
        sigma = params.k .* Xv_vec + 0.4 .* xi;
        sigma_hat_squared = (params.k .* Xv_vec + 0.4 .* xi_0_hat ) .* (params.k .* Xv_vec + 0.4 .* xi_0_hat .* cos(gamma));
        ratio = (params.R.^2 .* params.CT .* cos(gamma).^3) ./ (2 .* sigma_hat_squared);
        ratio = min(ratio, 0.9999);       
        C = 1 - sqrt(1 - ratio);
        dU_vec = C .* params.u_in(z_idx) .* exp(-(  (Yv_vec - y_c).^2   + (z-params.z_hub).^2  ) ./ (2 .* sigma.^2   ));
        dU = reshape(dU_vec, size(X));
        U = params.u_in(z_idx) - dU;
        flowfield(:,:,z_idx) = U;
    end

    % Plot current frame
    clf(fig);
    z_idx = 1;
    surf(X./params.D, Y./params.D, flowfield(:,:,z_idx) ./ params.u_hub);
    shading interp; view(2);
    xlabel('x/D'); ylabel('y/D');
    title(sprintf('Normalized velocity u/u_h | Yaw angle β = %.1f°', beta_current));
    colorbar; axis equal tight;
    caxis([0.4 1.0]); % Fixed colorbar range for consistency
    grid on;

    % Capture frame
    drawnow;
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % Write to GIF
    if frame_idx == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end

    fprintf('Frame %d/%d complete (β = %.1f°)\n', frame_idx, length(beta_values), beta_current);
end

fprintf('\nGIF saved as: %s\n', gif_filename);
close(fig);