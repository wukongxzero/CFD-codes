% Parameters
nx = 100;           % Number of grid points in x-direction
ny = 100;           % Number of grid points in y-direction
Lx = 1;             % Length of the domain in x-direction
Ly = 1;             % Length of the domain in y-direction
D = 0.1;            % Diffusion coefficient
u = 0.1;            % Velocity in x-direction
v = 0.2;            % Velocity in y-direction
bc = 0;             % Boundary condition value

% Grid
dx = Lx / (nx - 1);     % Grid spacing in x-direction
dy = Ly / (ny - 1);     % Grid spacing in y-direction
x = linspace(0, Lx, nx);% x-coordinate array
y = linspace(0, Ly, ny);% y-coordinate array

% Initialize solution matrix
T_cd = zeros(ny, nx);   % Solution using central difference scheme
T_upwind = zeros(ny, nx); % Solution using upwind scheme
%BOUNDARY CONDITIONS
T_cd(:, 1) = 100;        % Left boundary
T_cd(:, nx) = 200;       % Right boundary
T_cd(1, :) = 50;         % Bottom boundary
T_cd(ny, :) = 150;       % Top boundary

T_upwind(:, 1) = 100;    % Left boundary
T_upwind(:, nx) = 200;   % Right boundary
T_upwind(1, :) = 50;     % Bottom boundary
T_upwind(ny, :) = 150;   % Top boundary

% Central difference scheme
for j = 2:(ny-1)
    for i = 2:(nx-1)
        % Calculate the convective terms
        conv_x = u * (T_cd(j, i+1) - T_cd(j, i-1)) / (2 * dx);
        conv_y = v * (T_cd(j+1, i) - T_cd(j-1, i)) / (2 * dy);
        
        % Calculate the diffusion term
        diff = (T_cd(j, i+1) + T_cd(j, i-1) + T_cd(j+1, i) + T_cd(j-1, i) - 4 * T_cd(j, i)) / (dx^2 + dy^2);
        
        % Update the temperature at (i, j)
        T_cd(j, i) = T_cd(j, i) + D * (diff - conv_x - conv_y);
    end
end

% Upwind scheme
for j = 2:(ny-1)
    for i = 2:(nx-1)
        % Calculate the convective terms using upwind scheme
        if u >= 0
            conv_x = u * (T_upwind(j, i) - T_upwind(j, i-1)) / dx;
        else
            conv_x = u * (T_upwind(j, i+1) - T_upwind(j, i)) / dx;
        end
        
        if v >= 0
            conv_y = v * (T_upwind(j, i) - T_upwind(j-1, i)) / dy;
        else
            conv_y = v * (T_upwind(j+1, i) - T_upwind(j, i)) / dy;
        end
        
        % Calculate the diffusion term
        diff = (T_upwind(j, i+1) + T_upwind(j, i-1) + T_upwind(j+1, i) + T_upwind(j-1, i) - 4 * T_upwind(j, i)) / (dx^2 + dy^2);
        
        % Update the temperature at (i, j)
        T_upwind(j, i) = T_upwind(j, i) + D * (diff - conv_x - conv_y);
    end
end



% Plot the solutions
[X, Y] = meshgrid(x, y);

subplot(1, 2, 1);
surf(X, Y, T_cd);
xlabel('x');
ylabel('y');
zlabel('Temperature');
title('Central Difference Scheme');

subplot(1, 2, 2);
surf(X, Y, T_upwind);
xlabel('x');
ylabel('y');
zlabel('Temperature');
title('Upwind Scheme');
