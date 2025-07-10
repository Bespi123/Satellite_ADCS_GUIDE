function optimize_ekf_ga(app, disturbances, simParameters, time)
    % --- 2. Configuración del Algoritmo Genético Multi-Objetivo ---
    
    % El fitnessFcn ahora apunta a la función que devuelve múltiples objetivos
    fitnessFcn = @(params) ekfObjectiveFunction(params, app, disturbances, simParameters, time);
    
    % Número de variables: Diagonales de Q y R
    if simParameters.sensors.star.enable == 1
        nvars = 12; % 3 para Q y 9 para R
    else
        nvars = 9;  % 3 para Q y 6 para R
    end
    
    % Límites inferiores (Lower Bounds) y superiores (Upper Bounds)
    lb = 1e-9 * ones(1, nvars);  % Límite inferior
    ub = 1e-1 * ones(1, nvars);  % Límite superior (AUMENTADO para dar más espacio de búsqueda)
    
    % Condiciones Iniciales
    gyro = simParameters.ekf.gyro.std;
    acc  = simParameters.ekf.acc.std;
    mag  = simParameters.ekf.mag.std;
    if simParameters.sensors.star.enable == 1
        star = simParameters.ekf.star.std * 100;
        initial_population = [gyro; acc; mag; star]';
    else
        initial_population = [gyro; acc; mag]';
    end
    
    % Opciones para 'gamultiobj'
    % Se ajustan los parámetros para proteger y refinar las buenas soluciones
    options = optimoptions('gamultiobj', ...
        'PopulationSize', 5, ... % AUMENTADO: Más individuos para una mejor exploración
        'MaxGenerations', 10, ...  % AUMENTADO: Más tiempo para converger
        'ParetoFraction', 0.7, ... % AÑADIDO: Conserva una gran parte (70%) de las mejores soluciones (élite) en cada generación
        'CrossoverFraction', 0.7, ... % AÑADIDO: Reduce la agresividad del cruce para no destruir buenas soluciones
        'Display', 'iter', ...
        'PlotFcn', {@gaplotpareto}, ... % Gráfica específica para frentes de Pareto
        'UseParallel', true, ...
        'InitialPopulationMatrix', initial_population);
        
    % --- 3. Ejecutar el Algoritmo Genético Multi-Objetivo ---
    disp('Iniciando optimización multi-objetivo de parámetros del EKF...');
    
    % La llamada ahora es a gamultiobj
    % x: Matriz de soluciones en el frente de Pareto (cada fila es una solución)
    % fval: Matriz con los valores de los objetivos para cada solución en x
    [x, fval] = gamultiobj(fitnessFcn, nvars, [], [], [], [], lb, ub, [], options);
    
    % --- 4. Procesar y Mostrar Resultados ---
    disp('Optimización del EKF finalizada.');
    
    % ESTRATEGIA DE SELECCIÓN:
    % De todas las soluciones en el frente de Pareto, elegimos la que tiene
    % la menor SUMA de los RMSEs de todos los estados.
    
    % fval tiene N+1 columnas: [rmse1, rmse2, ..., rmseN, complejidad]
    % Sumamos solo las columnas de RMSE (todas menos la última)
    sum_of_rmses = sum(fval(:, 1:end-1), 2);
    
    [~, idx] = min(sum_of_rmses); % Encontramos el índice de la mejor solución
    
    optimal_params = x(idx, :);       % Parámetros de esa solución
    optimal_costs = fval(idx, :);     % Costos de esa solución
    
    fprintf('Se encontraron %d soluciones en el frente de Pareto.\n', size(x, 1));
    fprintf('Seleccionando la solución con la mejor suma de RMSEs.\n');
    disp('Parámetros óptimos seleccionados:');
    disp(optimal_params);
    fprintf('Costos de la solución seleccionada (RMSEs y Complejidad): %s\n', mat2str(optimal_costs, 4));

    % Actualiza los parámetros en la aplicación con la solución encontrada
    app.simulationParameters.ekf.equalModel = 0;
    app.simulationParameters.ekf.gyro.std = optimal_params(1:3);
    app.simulationParameters.ekf.acc.std = optimal_params(4:6);
    app.simulationParameters.ekf.mag.std = optimal_params(7:9);
    if simParameters.sensors.star.enable == 1
        app.simulationParameters.ekf.star.std = optimal_params(10:12);
    end
end
