function optimize_gains_ga(app,disturbances,simParameters,time)

% --- 2. Configuración del Algoritmo Genético ---
% Define la función objetivo
fitnessFcn = @(gains) objectiveFunction(gains, app, disturbances, simParameters, time);

% Número de variables a optimizar (6 para el controlador Feedback: 3 para P, 3 para K)
nvars = 6; 

% Límites inferiores (Lower Bounds) y superiores (Upper Bounds) para las ganancias
% ¡Estos límites son cruciales! Debes definirlos con sensatez.
lb = [0, 0, 0, 0, 0, 0]; % Valores mínimos para cada ganancia
ub = [10,   10,   10,   10,   10,   10];   % Valores máximos para cada ganancia

% %%% --- DEFINIR LA POBLACIÓN INICIAL --- %%%
% Define aquí uno o más conjuntos de ganancias que creas que son buenos.
% Cada fila es un individuo.
P = simParameters.feedback.Peye;
K = simParameters.feedback.Keye;

initial_population = [P(1,1), P(2,2), P(3,3), K(1,1), K(2,2), K(3,3)]; % Un buen punto de partida conocido

% Opciones del GA
options = optimoptions('ga', ...
    'PopulationSize', 10, ...       % Individuos por generación
    'MaxGenerations', 10, ...      % Número de generaciones
    'Display', 'iter', ...          % Muestra el progreso en cada iteración
    'PlotFcn', {@gaplotbestf, @gaplotstopping}, ... % Gráfica el mejor fitness y el criterio de parada
    'UseParallel', true, ...
    'InitialPopulationMatrix', initial_population);           % ¡MUY RECOMENDADO para acelerar el proceso!

% --- 3. Ejecutar el Algoritmo Genético ---
disp('Iniciando optimización con Algoritmo Genético...');
[optimal_gains, min_cost] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], options);

% --- 4. Mostrar Resultados ---
disp('Optimización finalizada.');
fprintf('Costo mínimo encontrado: %.4f\n', min_cost);
disp('Ganancias óptimas encontradas (P y K diagonales):');
disp(optimal_gains);

% Actualiza los parámetros en tu app con los resultados óptimos
app.simulationParameters.feedback.Peye = diag(optimal_gains(1:3));
app.simulationParameters.feedback.Keye = diag(optimal_gains(4:6));
msgbox('¡Ganancias óptimas encontradas y cargadas en la app!', 'Optimización Completa');