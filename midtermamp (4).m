% Create the thermal model
thermalmodel = createpde("thermal", "transient");

% Define the powder bed geometry with 11 layers
L = 0.004; % Length in meters (4 mm)
W = 0.0005; % Width in meters (0.5 mm)
H_total = 0.00022; % Total height in meters (0.22 mm, 11 layers of 20 microns each)

% Create the geometry for the subset powder bed
geometry = multicuboid(L, W, H_total);

% Assign geometry to the model
thermalmodel.Geometry = geometry;

% Generate a mesh with a smaller mesh size for increased accuracy
generateMesh(thermalmodel, 'Hmax', 0.0001); % Refine mesh size for better resolution in depth direction % Generate mesh with 20,000+ isoparametric linear hexahedron elements

% Define thermal properties for Nickel Alloy 625
k_solid = @(region, state) 5.331 + 0.015 * state.u; % Thermal conductivity for T ≤ TS
Cp_solid = @(region, state) 338.98 + 0.2437 * state.u; % Specific heat for T ≤ TS
k_liquid = 30.05; % Thermal conductivity for T ≥ TL
Cp_liquid = 735; % Specific heat for T ≥ TL
density = 8440 * 0.55; % Mass density in kg/m^3 adjusted for 55% powder packing density
TL = 1623; % Liquidus Temperature
TS = 1563; % Solidus Temperature
Lf = 227e3; % Latent heat of fusion [J/kg]

% Apply temperature-dependent thermal properties
thermalProperties(thermalmodel, ...
    'ThermalConductivity', @(region, state) (state.u <= TS) .* k_solid(region, state) + ...
                                          (state.u > TS & state.u < TL) .* ((k_liquid - k_solid(region, state)) .* ((state.u - TS) / (TL - TS)) + k_solid(region, state)) + ...
                                          (state.u >= TL) .* k_liquid, ...
    'MassDensity', density, ...
    'SpecificHeat', @(region, state) (state.u <= TS) .* Cp_solid(region, state) + ...
                                     (state.u > TS & state.u < TL) .* ((Cp_liquid - Cp_solid(region, state)) .* ((state.u - TS) / (TL - TS)) + Cp_solid(region, state) + Lf / (TL - TS)) + ...
                                     (state.u >= TL) .* Cp_liquid);

% Set initial temperature condition (ambient temperature)
initialTemperature = 293; % Initial temperature in K (20°C)
thermalIC(thermalmodel, initialTemperature); % Apply initial condition

% Assign boundary conditions
% Fixed temperature at the bottom face (Dirichlet boundary condition)
thermalBC(thermalmodel, 'Face', 1, 'Temperature', 353); % 80°C

% Convection heat loss at other surfaces (Neumann boundary condition)
convection_coefficient = 15; % Updated convection coefficient for better accuracy
thermalBC(thermalmodel, 'Face', 2:thermalmodel.Geometry.NumFaces, ...
    'ConvectionCoefficient', convection_coefficient, ...
    'AmbientTemperature', 293); % 20°C

% Define the Gaussian beam heat source function based on given formula
function Qflux = movingHeatSource(region, state, P, L, vs)
    % Inputs
    R = 0.7; % Reflectivity
    wo = 0.1e-3; % Waist size of the laser beam [m] (0.1 mm radius)
    
    % Define the movement of the laser beam along x-axis for a single pass
    % Updated to consider the time step size and scanning velocity
    x_center = vs * state.time; % Ensure laser moves continuously based on velocity
    x_center = min(x_center, L); % Ensure the beam does not exceed the powder bed length
    
    % Calculate distance from beam center
    r = sqrt((region.x - x_center).^2 + region.y.^2);
    
    % Gaussian beam heat flux based on given formula
    Qflux = (1 - R) * (2 * P / (pi * wo^2)) * exp(-2 * (r.^2 / wo^2));
end

% Solver settings for better convergence
thermalmodel.SolverOptions.RelativeTolerance = 1e-2; % Relax relative tolerance
thermalmodel.SolverOptions.AbsoluteTolerance = 1e-3; % Increase absolute tolerance

% Define process parameters for low, medium, and high density
energy_levels = {'Low', 'Medium', 'High'};
power_levels = [169, 195, 195]; % Laser powers in [W]
scan_velocities = [975e-3, 800e-3, 725e-3]; % Scanning velocities in [m/s]
melt_pool_dimensions = [];

% Calculate appropriate time step size to avoid skipping elements
dx = 0.03e-3; % Smallest element dimension in the scanning axis [m]

for i = 1:length(power_levels)
    % Set parameters for current energy level
    P = power_levels(i);
    vs = scan_velocities(i);
    dt = dx / vs; % Calculate time step size [s]
    
    % Define the time vector for simulation
    simulation_time = 0.2; % Increase simulation time for better heat penetration % Total simulation time in seconds
    tlist = linspace(0, simulation_time, round(simulation_time / dt)); % Create a row vector of linearly spaced time steps
    
    % Update the heat source function for the given power level
    thermalBC(thermalmodel, 'Face', 2, 'HeatFlux', @(region, state) movingHeatSource(region, state, P, L, vs), 'Vectorized', 'on');
    
    % Solve the thermal model for the given power level
    result = solve(thermalmodel, tlist); % Use refined time steps to avoid skipping elements
    
    % Extract temperature data at the final time step
    temp = result.Temperature(:, end);
    nodes = thermalmodel.Mesh.Nodes;
    x = nodes(1, :); % x-coordinates
    y = nodes(2, :); % y-coordinates
    z = nodes(3, :); % z-coordinates
    
    % Find nodes where temperature exceeds solidus temperature (TS)
    melt_pool_nodes = temp > TS;
    x_melt = x(melt_pool_nodes);
    y_melt = y(melt_pool_nodes);
    z_melt = z(melt_pool_nodes);
    
    % Calculate melt pool length (MPL), width (MPW), and depth (MPD)
    if isempty(x_melt)
        MPL = 0;
        MPW = 0;
        MPD = 0;
    else
        MPL = max(x_melt) - min(x_melt); % Length along x-direction
        MPW = max(y_melt) - min(y_melt); % Width along y-direction
        MPD = max(z_melt) - min(z_melt); % Depth along z-direction
    end
    
    melt_pool_dimensions = [melt_pool_dimensions; {energy_levels{i}}, MPL, MPW, MPD];
    
    % Plot temperature distribution at final time step for visualization
    figure;
    pdeplot3D(thermalmodel, 'Mesh', 'on', 'ColorMap', temp);
    caxis([500 TL]); % Set temperature scale to limit to T ≤ TL for better visualization
    colorbar;
    title(['Temperature Distribution - ', energy_levels{i}, ' Energy Density']);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
end

% Display melt pool dimensions
disp('Energy Level | Melt Pool Length (mm) | Melt Pool Width (mm) | Melt Pool Depth (mm)');
disp(melt_pool_dimensions);

% Display temperature along scanning direction for validation purposes
figure;
scatter3(x, y, z, 10, temp, 'filled');
colorbar;
caxis([500, TL]); % Limit the temperature scale
title('Temperature Distribution Along Scanning Direction');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

% Validate melt pool length, width, depth
melt_pool_nodes = temp > TS;
x_melt = x(melt_pool_nodes);


if isempty(x_melt)
    MPL = 0; MPW = 0; MPD = 0;
else
    MPL = max(x_melt) - min(x_melt);
    MPW = max(y_melt) - min(y_melt);
    MPD = max(z_melt) - min(z_melt);
end

fprintf('Melt Pool Dimensions: \n');
fprintf('Length (mm): %f \n', MPL * 1e3);
fprintf('Width (mm): %f \n', MPW * 1e3);
fprintf('Depth (mm): %f \n', MPD * 1e3);

% Validate results with experimental values and compare the thermal fields
% Table summarizing the melt pool results

disp('Energy Density [J/mm^3] | Melt Pool Length (mm) | Melt Pool Width (mm) | Melt Pool Depth (mm)');
disp(melt_pool_dimensions);
