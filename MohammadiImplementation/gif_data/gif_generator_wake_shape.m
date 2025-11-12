 %Script to generate animated GIF of normalized velocity deficit contours
% Creates a GIF showing how wake deficit evolves at different downstream positions

% Setup for GIF creation
gif_filename = 'wake_velocity_deficit.gif';
delay_time = 0.5;  % seconds between frames

% Loop through x-slices and create frames
for idx = 1:length(x_slices)
    xi = x_slices(idx);
    
    % Normalized delta U slice at fixed x
    dUz = squeeze(delta_flowfield(xi, :, :)) ./ params.u_hub;  % Ny x Nz
    
    % Choose smooth levels (adjust nLevels to taste)
    nLevels = 14;
    levels = linspace(min(dUz, [], 'all'), max(dUz, [], 'all'), nLevels);
    
    % Create figure
    figure('Color', 'w', 'Position', [100, 100, 800, 600]); 
    hold on;
    
    % Filled contours (no edges) for smooth color bands
    contourf(Ygrid./params.D, Zgrid./params.D, dUz, levels, 'LineStyle', 'none');
    colormap(turbo); 
    clim([levels(1) levels(end)]); 
    colorbar;
    
    % Crisp black isolines on top
    contour(Ygrid./params.D, Zgrid./params.D, dUz, levels, 'LineColor', 'k', 'LineWidth', 0.8);
    
    axis equal tight; 
    grid on;
    xlabel('y/D', 'FontSize', 12); 
    ylabel('z/D', 'FontSize', 12);
    
    % Two-line title (main + subtitle)
    t = sprintf('Normalized velocity deficit ΔU/u_h at x/D = %.1f', ...
                params.xrange(xi)/params.D);
    s = sprintf('Veer: %.2f °/m | Yaw: %.1f°', ...
                params.alpha_gradient, rad2deg(params.beta));
    title({t, s}, 'FontSize', 12);
    
    % Capture frame for GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to GIF file
    if idx == 1
        % Create new GIF file with first frame
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        % Append subsequent frames
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
    
    % Close figure to save memory
    close(gcf);
    
    % Display progress
    fprintf('Processed frame %d of %d (x/D = %.1f)\n', idx, length(x_slices), params.xrange(xi)/params.D);
end

fprintf('\nGIF creation complete! Saved as: %s\n', gif_filename);

%% Optional: Display the GIF in MATLAB (requires Image Processing Toolbox)
% figure('Name', 'Velocity Deficit Animation');
% imshow(gif_filename);