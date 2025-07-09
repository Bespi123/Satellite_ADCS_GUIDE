function cost = objectiveFunction(gains, app, disturbances, simParameters, time)
    % objectiveFunction - Evalúa el costo de un conjunto de ganancias del controlador.
    %
    % Inputs:
    %   gains         - Vector con las ganancias a optimizar.
    %   app           - Objeto de la app (para pasar al simulador).
    %   disturbances  - Perturbaciones.
    %   simParameters - Parámetros base de la simulación.
    %   time          - Parámetros de tiempo.
    %
    % Output:
    %   cost          - Un único valor escalar que representa el "costo". Menor es mejor.

    % --- 1. Actualizar los parámetros de simulación con las nuevas ganancias ---
    % Asumimos que optimizaremos el controlador Feedback.
    % 'gains' es un vector: [P11, P22, P33, K11, K22, K33]
    simParameters.feedback.Peye = diag(gains(1:3));
    simParameters.feedback.Keye = diag(gains(4:6));

    % --- 2. Ejecutar la simulación ---
    try
        [~, ~, ~, ~, indicators, ~, error_flag] = simulation_rk4(app, disturbances, simParameters, time);

        % Si la simulación falla (se vuelve inestable), asigna un costo muy alto.
        if error_flag == 1
            cost = 1e10; % Penalización alta
            return;
        end
    catch
        % Si ocurre cualquier otro error
        cost = 1e10;
        return;
    end

    % --- 3. Calcular el costo a partir de los indicadores ---
    % Ponderamos cada indicador para crear un costo total.
    % ¡Estos pesos (w1, w2, w3) son muy importantes y debes ajustarlos!
    w1 = 0.5;  % Peso para el tiempo de estabilización (ts)
    w2 = 0.3;  % Peso para el gasto de control (ASCCT)
    w3 = 0.2;  % Peso para el error integrado (EULERINT)

    cost_ts = indicators.ts;
    cost_ascct = indicators.ascct(end); % Valor final del gasto de control
    cost_euler = indicators.eulerInt(end); % Valor final del error integrado

    % Combinar en un único costo
    cost = w1 * cost_ts + w2 * cost_ascct + w3 * cost_euler;
    
    % Asegurarse de que el costo no sea NaN
    if isnan(cost)
        cost = 1e10; % Penalización si el costo es inválido
    end
    
    % Opcional: Mostrar el progreso en la consola
    fprintf('Gains: [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f] -> Cost: %.4f\n', gains, cost);
end