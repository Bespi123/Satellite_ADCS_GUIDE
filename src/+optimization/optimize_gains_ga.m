function optimal_gains = optimize_gains_ga(~, disturbances, simParameters, time)
    import optimization.*

    % optimize_gains_ga - Orchestrates the optimization of controller gains
    %                     using a Genetic Algorithm (GA).
    %
    % This function configures and runs the MATLAB 'ga' optimizer to find
    % the set of gains (P and K) that minimizes a defined cost function.
    %
    % Inputs:
    %   app           - Handle to the main application object (GUI).
    %   disturbances  - Structure with the disturbances for the simulation.
    %   simParameters - Structure with the base simulation parameters.
    %   time          - Structure with the simulation time parameters.

    %Disable AHRS 
    simParameters.ahrs.enable = 0;

    % --- 1. Genetic Algorithm Configuration ---

    % Define the objective function (or fitness function).
    % A function handle (@) is used to pass the 'objectiveFunction'
    % to the optimizer. The optimizer will call this function repeatedly.
    fitnessFcn = @(gains) objectiveFunction(gains, disturbances, simParameters, time);

    % Number of variables to optimize.
    if simParameters.controller.selector == 1
        % In this case, 6 gains: 3 for the
        % P matrix (proportional) and 3 for the K matrix (derivative).
        nvars = 6;
        % Lower Bounds and Upper Bounds for the gains.
        % These limits define the "search space". It is crucial to define them
        % wisely to guide the algorithm.
        lb = [0, 0, 0, 0, 0, 0]; % Minimum allowed values for each gain.
        ub = [10, 10, 10, 10, 10, 10]; % Maximum allowed values for each gain.
        % --- Define Initial Population ---
        % A starting point can be provided to the GA. This is very useful if you
        % already have an idea of what good gains are.
        % Each row of 'initial_population' represents an initial "individual".
        P = simParameters.controller.feedback.Peye;
        K = simParameters.controller.feedback.Keye;
        initial_population = [P(1,1), P(2,2), P(3,3), K(1,1), K(2,2), K(3,3)]; % Use the current gains from the app as a starting point.
    else
        % In this case, delta, gamma, k0
        nvars = 3; 
        % Lower Bounds and Upper Bounds for the gains.
        % These limits define the "search space". It is crucial to define them
        % wisely to guide the algorithm.
        lb = [0, 0, 0]; % Minimum allowed values for each gain.
        ub = [10, 1, 10]; % Maximum allowed values for each gain.
        % --- Define Initial Population ---
        % A starting point can be provided to the GA. This is very useful if you
        % already have an idea of what good gains are.
        % Each row of 'initial_population' represents an initial "individual".
        delta = simParameters.controller.boskController.delta;
        gamma = simParameters.controller.boskController.gamma;
        k0    = simParameters.controller.boskController.k0;
        initial_population = [delta, gamma, k0]; % Use the current gains from the app as a starting point.
    end


    % GA Options: Detailed configuration of the algorithm's behavior.
    options = optimoptions('ga', ...
        'PopulationSize', 10, ...       % Number of individuals (solutions) in each generation.
        'MaxGenerations', 20, ...      % Stopping criterion: maximum number of generations.
        'Display', 'iter', ...          % Displays progress in the command window at each iteration.
        'PlotFcn', {@gaplotbestf, @gaplotstopping}, ... % Generates real-time plots:                                               
        'UseParallel', false, ...        
        'InitialPopulationMatrix', initial_population); % Provides the previously defined initial population.

    % --- 2. Run the Genetic Algorithm ---
    msgbox('Starting optimization with Genetic Algorithm...','Run optimization');
    
    % Temporarily turn off the specific warning
    warning('off', 'MATLAB:nearlySingularMatrix');

    % Call to the 'ga' function. This is the main execution of the optimizer.
    % It returns the best gains found and the associated cost.
    [optimal_gains, min_cost] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], options);

    % Turn the warning back on
    warning('on', 'MATLAB:nearlySingularMatrix');

    % --- 3. Display Results ---
    disp('Optimization finished.');
    fprintf('Minimum cost found: %.4f\n', min_cost);
    disp('Optimal gains found:');
    disp(optimal_gains);
    
    % Notifies the user that the process has finished.
    msgbox('Optimal gains found and loaded into the app!', 'Optimization Complete');
end