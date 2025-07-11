function costs_vector = ekfObjectiveFunction(params, app, disturbances, simParameters, time)
    % ekfObjectiveFunction_per_state - Evalúa el costo de CADA estado por separado.
    %
    % Output:
    %   costs_vector - Un vector fila con el RMSE promedio para cada estado.

    num_simulations = 5;
    
    % Inicializar una matriz para guardar los costos de cada estado en cada simulación
    % El número de columnas debe coincidir con el número de estados en indicators.RSME
    % Asumimos un número de estados (ej. 4 para quaternions), ajústalo si es necesario.
    num_states = 1; 
    total_costs_matrix = zeros(num_simulations, num_states);

    % --- Extraer parámetros ---
    simParameters.ekf.gyro.std = params(1:3)';
    simParameters.ekf.acc.std  = params(4:6)';
    simParameters.ekf.mag.std  = params(7:9)';
    if simParameters.sensors.star.enable == 1
        simParameters.ekf.star.std  = params(10:12)';
    end
    simParameters.ekf.enable = 0;
    simParameters.ekf.equalModel = 0;

    % --- Bucle de Montecarlo ---
    for i = 1:num_simulations
        try
            [~, ~, ~, ~, indicators, ~, error_flag] = simulation_rk4(app, disturbances, simParameters, time);
            
            if error_flag
                % Si hay un error, penaliza todos los estados por igual
                total_costs_matrix(i, :) = 1e10;
            else
                % Guarda el vector de RMSE del final de la simulación
                current_run_costs = indicators.RMSE;
                current_run_costs(isnan(current_run_costs)) = 1e10; % Penaliza si hay NaN
                total_costs_matrix(i, :) = current_run_costs;
            end
        catch
            % Penalización si la simulación falla por completo
            total_costs_matrix(i, :) = 1e10;
        end
    end
    
    % --- Calcular el costo promedio para cada estado ---
    % 'mean' sobre la primera dimensión promedia las 5 simulaciones
    avg_costs_per_state = mean(total_costs_matrix, 1);
    
    % El costo final es un vector fila
    costs_vector = avg_costs_per_state * 100;
    
    fprintf('Costos Promedio por Estado (RMSE): %s\n', mat2str(costs_vector, 4));
end