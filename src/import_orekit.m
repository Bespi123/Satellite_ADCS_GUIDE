%% 1. Configuración del Entorno
disp('Configurando el entorno de Orekit...');

% Definir las rutas a las carpetas (ajústalas si es necesario)
lib_folder = 'orekit';
data_folder = 'orekit/orekit-data';

% --- Añadir todas las librerías .jar al path de Java ---
% Librería principal de Orekit
javaaddpath(fullfile(lib_folder, 'orekit-13.0.jar'));
% Módulos de la librería Hipparchus (versión 3.0)
javaaddpath(fullfile(lib_folder, 'hipparchus-clustering-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-core-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-fft-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-filtering-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-fitting-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-geometry-4.0.2.jar'));
%javaaddpath(fullfile(lib_folder, 'hipparchus-migration-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-ode-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-optim-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-samples-4.0.2.jar'));
javaaddpath(fullfile(lib_folder, 'hipparchus-stat-4.0.2.jar'));

% --- Configurar el proveedor de datos de Orekit ---
import org.orekit.data.*;
import java.io.File;
manager = DataContext.getDefault().getDataProvidersManager();
manager.clearProviders(); % Buena práctica para evitar rutas duplicadas
manager.addProvider(DirectoryCrawler(File(data_folder)));

disp('Entorno listo.');
fprintf('\n');