% Parameters
nx = 100;           % Number of grid points in x-direction
L = 1;              % Length of the domain
u = 0.5;            % Velocity
bc_left = 2;        % Boundary condition value at the left boundary
bc_right = 0;       % Boundary condition value at the right boundary

% Grid
dx = L / (nx - 1);   % Grid spacing in x-direction
x = linspace(0, L, nx); % x-coordinate array

% Time parameters
dt = 0.001;         % Time step size
tEnd = 1;           % End time
t = 0;              % Initial time

% Initialize solution vector
T = zeros(1, nx);

% Boundary conditions
T(1) = bc_left;     % Left boundary
T(nx) = bc_right;   % Right boundary

% Time stepping loop
while t < tEnd
    % Create a copy of the solution vector for the next time step
    T_new = T;
    
    % Loop over interior points
    for i = 2:(nx-1)
        % Calculate the convection term
        convection = u * (T(i) - T(i-1)) / dx;
        
        % Update the temperature at (i)
        T_new(i) = T(i) - dt * convection;
    end
    
    % Update the solution vector for the next time step
    T = T_new;
    
    t = t + dt;     % Update time
end

% Plot the final solution
plot(x, T, 'b', 'LineWidth', 2);
xlabel('x');
ylabel('Temperature');
title('Unsteady 1D Convection');
