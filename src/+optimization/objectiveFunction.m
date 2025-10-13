function cost = objectiveFunction(gains, disturbances, simParameters, time)
    import core.*
    % objectiveFunction - Evaluates the quality (cost) of a set of controller gains.
    %
    % This function is designed to be used by an optimizer (such as a
    % genetic algorithm) that searches for the best controller gains.
    % A lower cost means better controller performance.
    %
    % Inputs:
    %   gains         - Row vector (1x6) with the gains to be optimized.
    %                 Ex: [Kp1, Kp2, Kp3, Kd1, Kd2, Kd3]
    %   app           - Handle to the main application object (GUI), necessary
    %                 to pass to the simulator.
    %   disturbances  - Structure or matrix with external disturbances for the simulation.
    %   simParameters - Structure with all the base simulation parameters.
    %                 This structure will be temporarily modified with the new 'gains'.
    %   time          - Structure with the simulation time parameters (dt, t_final).
    %
    % Output:
    %   cost          - A single scalar value representing the total "cost"
    %                 of the controller's performance. A lower value is better.

    %Disable AHRS 
    simParameters.ahrs.enable = 0;

    if simParameters.controller.selector == 1
        % --- 1. Update Simulation Parameters with New Gains ---
        % The optimizer passes a 'gains' vector. Here, we extract those values and
        % place them in the controller gain matrices within the
        % simulation parameters structure.
        % It is assumed that the gains are for the diagonals of the P and K matrices.
        simParameters.feedback.Peye = diag(gains(1:3)); % Proportional gain matrix (P)
        simParameters.feedback.Keye = diag(gains(4:6)); % Derivative gain matrix (K)
    else
        % --- 1. Update Simulation Parameters with New Gains ---
        % The optimizer passes a 'gains' vector. Here, we extract those values and
        % place them in the controller gain matrices within the
        % simulation parameters structure.
        % It is assumed that the gains are for the diagonals of the P and K matrices.
        simParameters.boskController.delta = gains(1); % delta
        simParameters.boskController.gamma = gains(2); % gamma
        simParameters.boskController.k0    = gains(3); % k0
    end

    % --- 2. Run the Simulation ---
    % A try-catch block is used to handle unexpected errors. If a
    % combination of gains causes the simulation to fail (e.g., due to
    % numerical instability), the optimizer will not stop.
    try
        %Dump variable
        App = NaN;
        % Call to the main simulation function.
        [~, ~, ~, ~, indicators, ~, error_flag] = simulation_rk4(App, disturbances, simParameters, time);
        
        % Check if the simulation itself reported an error (e.g., divergence).
        if error_flag == 1
            cost = 1e10; % Assign a very high penalty to discard this solution.
            return;      % Ends the execution of this function.
        end
    catch
        % If any other error occurs during the simulation execution.
        cost = 1e10; % Assign the same high penalty.
        return;
    end

    % --- 3. Calculate the Cost from Performance Indicators ---
    % This is the most critical part. Several performance metrics are combined
    % into a single "cost" number. The weighting of each metric determines
    % which aspect of performance is more important.
    
    % Weights: You must adjust these values according to your priorities!
    % The sum of the weights does not have to be 1, but it helps with interpretation.
    w1 = 0.5;  % Weight for the settling time (ts). Prioritizes speed.
    w2 = 0.3;  % Weight for the control effort (ASCCT). Prioritizes energy efficiency.
    w3 = 0.2;  % Weight for the integrated error (EULERINT). Prioritizes accuracy over time.
    
    % Extract the final values of the indicators from the simulation.
    cost_ts = indicators.ts;               % Time it takes for the system to stabilize.
    cost_ascct = indicators.ascct(end);    % Accumulated control effort at the end.
    cost_euler = indicators.eulerInt(end); % Integrated Euler error at the end.
    
    % Cost Formula: Weighted combination of the metrics.
    cost = w1 * cost_ts + w2 * cost_ascct + w3 * cost_euler;
    
    % --- 4. Sanity Check ---
    % Ensure the calculated cost is not an invalid value like NaN (Not-a-Number).
    if isnan(cost)
        cost = 1e10; % Penalty if the cost is invalid
    end
    
    % --- 5. Visualization (Optional) ---
    % Print the result of this evaluation to the console. It is useful for
    % monitoring the optimizer's progress in real-time.
    %fprintf('Gains: [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f] -> Cost: %.4f\n', gains, cost);
end