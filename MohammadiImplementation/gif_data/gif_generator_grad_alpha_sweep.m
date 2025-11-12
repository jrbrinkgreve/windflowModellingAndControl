% This script sweeps the alpha_gradient parameter from -0.2 to 0.2 deg/m
% and generates an animated GIF showing how the wake evolves

close all; clear all; clc;

function params = load_params(gradient_alpha_val)
    params.numxpoints = 10;
    params.numypoints = 100;
    params.numzpoints = 100;
    params.R = 63; %rotor radius
    params.D = params.R * 2;
    params.z_hub = 90; %hub height
    params.xlim_min = 1.*params.D;
    params.xlim_max = 10.*params.D;
    params.ylim_min = -1.5.*params.D;
    params.ylim_max = 1.5.*params.D;
    params.zlim_min = 0;
    params.zlim_max = params.D*2;
    params.xrange = linspace(params.xlim_min, params.xlim_max, params.numxpoints);
    params.yrange = linspace(params.ylim_min, params.ylim_max, params.numypoints);
    params.zrange = linspace(params.zlim_min, params.zlim_max, params.numzpoints);
    params.CT = 0.66; %torque coeff.
    params.beta = deg2rad(-25); %yaw angle
    params.u_hub = 8; %m./s, hub inflow vel.
    params.u_in = linspace(8,8, params.numzpoints); %distribution of inflow vel over z
    params.u_star = 0.4; %friction velocity
    params.lambda = 7.5; %tip-speed ratio
    params.k = 0.6 .* params.u_star ./ params.u_hub;
    params.alpha_gradient = gradient_alpha_val; % PARAMETER BEING SWEPT
    params.z_relative = params.zrange - params.z_hub;
    params.alpha = deg2rad(params.alpha_gradient * params.z_relative);
    params.angle_tolerance = 1e-10;
end

%% Main computation and GIF generation
tic;

% Define gradient_alpha sweep range
gradient_alpha_values = linspace(-0.2, 0.2, 50); % 20 frames
num_frames = length(gradient_alpha_values);

% Create figure for capturing frames
fig = figure('Color','w', 'Position', [100, 100, 1000, 700]);
gif_filename = 'gradient_alpha_sweep.gif';

% Preallocate for GIF writing
first_frame = true;

for frame_idx = 1:num_frames
    fprintf('Processing frame %d of %d (gradient_alpha = %.3f deg/m)\n', ...
            frame_idx, num_frames, gradient_alpha_values(frame_idx));

    % Load parameters with current gradient_alpha value
    params = load_params(gradient_alpha_values(frame_idx));

    % Create grid
    [X, Y] = ndgrid(params.xrange, params.yrange);

    % Preallocate flowfield
    flowfield = zeros(params.numxpoints, params.numypoints, params.numzpoints);
    delta_flowfield = zeros(params.numxpoints, params.numypoints, params.numzpoints);

    % Vectorize for faster operations
    X_vec = X(:);
    Y_vec = Y(:);

    % Compute flowfield for all z heights
    for z_idx = 1:params.numzpoints
        z = params.zrange(z_idx);
        alpha = params.alpha(z_idx);
        gamma = params.beta + alpha;

        % eq 3
        A_star = (1 + sqrt(1- params.CT .* cos(gamma).^2) ) ./  (2.*sqrt(1- params.CT .* cos(gamma).^2));

        % eq 4
        xi_0_hat = params.R .* sqrt(A_star);

        % build rotation matrix for x,y at each z
        rotmtx = [cos(alpha), sin(alpha); ...
                 -sin(alpha), cos(alpha)];

        % veered directions
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
        y_hat_c = (num ./ den) .* sgn_t_hat ...
                  - (2./pi) .* t_hat ./ (((z + params.z_hub)./xi_0_hat).^2 - 1);

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

            % eq 11: analytical params ~ fourier expansion
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
        delta_flowfield(:,:,z_idx) = dU;
    end

    % Extract YZ plane slice at x_index = 6
    [Ygrid, Zgrid] = ndgrid(params.yrange, params.zrange);
    x_index = 6;
    dUz = squeeze(delta_flowfield(x_index, :, :)) ./ params.u_hub;

    % Clear previous plot
    clf;

    % Create contour plot
    nLevels = 14;
    levels = linspace(min(dUz,[],'all'), max(dUz,[],'all'), nLevels);

    hold on;
    contourf(Ygrid./params.D, Zgrid./params.D, dUz, levels, 'LineStyle','none');
    colormap(turbo);
    clim([levels(1) levels(end)]);
    cbar = colorbar;
    cbar.Label.String = 'ΔU/u_h';

    contour(Ygrid./params.D, Zgrid./params.D, dUz, levels, 'LineColor','k', 'LineWidth',0.8);

    axis equal tight; grid on;
    xlabel('y/D', 'FontSize', 12);
    ylabel('z/D', 'FontSize', 12);

    % Title with current gradient_alpha value
    t = sprintf('Normalized velocity deficit ΔU/u_h at x/D = %.1f', ...
                params.xrange(x_index)/params.D);
    s = sprintf('Veer gradient: %.3f °/m | Yaw: %.1f°', ...
                gradient_alpha_values(frame_idx), rad2deg(params.beta));
    title({t, s}, 'FontSize', 14);

    % Capture frame for GIF
    drawnow;
    frame = getframe(fig);
    img = frame2im(frame);
    [img_indexed, img_map] = rgb2ind(img, 256);

    % Write to GIF
    if first_frame
        imwrite(img_indexed, img_map, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
        first_frame = false;
    else
        imwrite(img_indexed, img_map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end

toc;
fprintf('\nGIF saved as: %s\n', gif_filename);
fprintf('Total frames: %d\n', num_frames);
fprintf('Gradient alpha range: [%.2f, %.2f] deg/m\n', min(gradient_alpha_values), max(gradient_alpha_values));
