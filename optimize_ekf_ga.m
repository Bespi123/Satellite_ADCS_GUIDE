function optimal_params = optimize_ekf_ga(~, disturbances, simParameters, time)
    % optimize_ekf_ga - Orchestrates the optimization of EKF parameters
    %                     using a single-objective Genetic Algorithm (GA).
    %
    % This script configures and runs the MATLAB 'ga' optimizer. Its goal
    % is to find the set of EKF parameters that minimizes the scalar cost
    % returned by the 'ekfObjectiveFunction'.

    % --- 1. Genetic Algorithm Configuration ---
    
    % Define the fitness function. The 'ga' optimizer will call this
    % function for each evaluation, passing it a new vector of parameters.
    fitnessFcn = @(params) ekfObjectiveFunction(params, disturbances, simParameters, time);
    
    % Number of variables to optimize (the EKF parameters).
    if simParameters.sensors.star.enable == 1
        nvars = 12; % 3 for gyro, 3 for acc, 3 for mag, 3 for star sensor
    else
        nvars = 9;  % 3 for gyro, 3 for acc, 3 for mag
    end
    
    % Lower Bounds and Upper Bounds.
    % These define the search space for the parameters.
    lb = 1e-9 * ones(1, nvars);  % Lower bound (very small values, but not zero)
    ub = 1e-1 * ones(1, nvars);  % Upper bound. If the optimal results are "stuck" at this
                                 % limit, consider increasing it.
    
    % Initial Conditions (Initial Population).
    % Providing a good starting point can accelerate convergence.
    gyro = simParameters.ekf.gyro.std;
    acc  = simParameters.ekf.acc.std;
    mag  = simParameters.ekf.mag.std;
    if simParameters.sensors.star.enable == 1
        star = simParameters.ekf.star.std;
        initial_population = [gyro; acc; mag; star]';
    else
        initial_population = [gyro; acc; mag]';
    end
    
    % Options for the 'ga' optimizer.
    % THIS IS THE MOST IMPORTANT PART FOR GETTING GOOD RESULTS!
    options = optimoptions('ga', ...
        'PopulationSize', 5, ... % RECOMMENDATION: Increase. A value like 5 is too low to explore. Try 50 or 100.
        'MaxGenerations', 10, ... % RECOMMENDATION: Increase. 10 generations is very few. Try 50 or 100.
        'Display', 'iter', ...    % Displays progress in the console.
        'PlotFcn', {@gaplotbestf, @gaplotstopping}, ... % Displays real-time progress plots.
        'UseParallel', true, ...  % Speeds up optimization if you have the Parallel Computing Toolbox.
        'InitialPopulationMatrix', initial_population); % Uses your current parameters as a starting point.
        
    % --- 2. Run the Genetic Algorithm ---
    disp('Starting EKF parameter optimization...');
    
    % Main call to the optimization library.
    [optimal_params, min_cost] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], options);
    
    % --- 3. Display and Apply Results ---
    disp('Optimization finished.');
    fprintf('Minimum cost found: %.4f\n', min_cost);
    disp('Optimal parameters found:');
    disp(optimal_params);
    
end