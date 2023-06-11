% Parameters
nx = 100;        % Number of grid points in x-direction
ny = 100;        % Number of grid points in y-direction
Lx = 1;          % Length of the domain in x-direction
Ly = 1;          % Length of the domain in y-direction
D = 0.1;         % Diffusion coefficient
bc = 0;          % Boundary condition value

% Grid
dx = Lx / (nx - 1);     % Grid spacing in x-direction
dy = Ly / (ny - 1);     % Grid spacing in y-direction
x = linspace(0, Lx, nx);% x-coordinate array
y = linspace(0, Ly, ny);% y-coordinate array

% Initialize solution matrix
T = zeros(ny, nx);

bc_left = -5;     % Gradient on the left boundary
bc_right = 10;    % Gradient on the right boundary
bc_bottom = 100;    % Gradient on the bottom boundary
bc_top = 0;       % Gradient on the top boundary

T(:, 1) = T(:, 2) - dx * bc_left;              % Left boundary
T(:, nx) = T(:, nx-1) + dx * bc_right;         % Right boundary
T(1, :) = T(2, :) - dy * bc_bottom;             % Bottom boundary
T(ny, :) = T(ny-1, :) + dy * bc_top;            % Top boundary
% Finite difference method
for j = 2:(ny-1)
    for i = 2:(nx-1)
        % Calculate the Laplacian
        laplacian = (T(j, i+1) + T(j, i-1) + T(j+1, i) + T(j-1, i) - 4 * T(j, i)) / (dx^2 + dy^2);
        
        % Update the temperature at (i, j)
        T(j, i) = T(j, i) + D * laplacian;
    end
end

% Plot the solution
[X, Y] = meshgrid(x, y);
surf(X, Y, T);
xlabel('x');
ylabel('y');
zlabel('Temperature');
title('Steady Diffusion in 2D with Neumann Boundary Conditions');
