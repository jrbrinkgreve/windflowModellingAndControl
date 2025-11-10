%implementation of Mohammadi et al.:
%https:././www.mdpi.com./1996-1073./15./23./9135

%parameter definition:


function params = load_params()
    params.numypoints = 100;
    params.numxpoints = 100;
    params.CT = 0.66; %torque coeff.
    params.beta = deg2rad(20); %yaw angle
    %yaw angle is assumed to be positive if the rotor is misaligned in the anticlockwise direction, seen from above.
    params.R = 63; %rotor radius
    params.u_hub = 8; %m./s, hub inflow vel.
    params.u_in = params.u_hub; %linspace(6,10, params.numypoints); %distribution of inflow vel over z
    params.u_star = 0.4; %friction velocity
    %(x,y,z) locations from turbine and from ground
    params.z_hub = 90; %hub height
    params.lambda = 7.5; %tip-speed ratio
    params.k = 0.6 .* params.u_star ./ params.u_hub;
end



params = load_params();

xrange = linspace(2*params.R, 10.*params.R, params.numxpoints);
yrange = linspace(-3.*params.R, 3.*params.R, params.numypoints);

[X, Y] = meshgrid(xrange, yrange);
z = params.z_hub;



%eq 3
A_star = (1 + sqrt(1- params.CT .* cos(params.beta).^2)) ./  (2.*sqrt(1- params.CT .* cos(params.beta).^2));

%eq 4
xi_0_hat = params.R .* sqrt(A_star);

%eq 5
t_hat = -1.44 .* (params.u_hub ./ params.u_star) .* (params.R ./ xi_0_hat) .* ...
         (params.CT .* cos(params.beta)^2 .* sin(params.beta)) .* ...
         (1 - exp(-0.35 .* (params.u_star ./ params.u_in) .* (X ./ params.R)));


%eq 6
sgn_t_hat = sign(t_hat);
abs_t_hat = abs(t_hat);

num = (pi - 1) .* abs_t_hat.^3 + 2 .* sqrt(3).*pi^2 .* t_hat.^2 + 48.*(pi - 1)^2 .* abs_t_hat;
den = 2.*pi.*(pi - 1).*t_hat.^2 + 4 .* sqrt(3) .*pi^2 .* abs_t_hat + 96.*(pi - 1)^2;


y_hat_c = (num ./ den) .* sgn_t_hat ...
          - (2./pi) .* t_hat ./ (((z + params.z_hub)./xi_0_hat).^2 - 1);


%eq 7
y_c = y_hat_c .* xi_0_hat;




%skip eq 8-12: we use the wake shape approximation for now

sigma_hat_squared = (params.k .* X + 0.4 .* xi_0_hat ) .* (params.k .* X + 0.4 .* xi_0_hat .* cos(params.beta));
ratio = (params.R.^2 .* params.CT .* cos(params.beta).^3) ./ (2 .* sigma_hat_squared);
ratio = min(ratio, 0.9999);     
C = 1 - sqrt(1 - ratio);

dU = C .* params.u_hub .* exp(-(  (Y - y_c).^2   + (z-params.z_hub).^2  ) ./ (2 .* sigma_hat_squared   ));

U = params.u_hub - dU;


figure('Color','w');
surf(X./params.R, Y./params.R, dU ./ params.u_hub);
shading interp; view(2);
xlabel('x/R'); ylabel('y/R');
title(['Normalized velocity deficit \Deltau./u_h (yaw = ', num2str(rad2deg(params.beta)), '°)']);
colorbar; axis equal tight;
grid on;

