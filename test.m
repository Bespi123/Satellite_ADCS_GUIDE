clear 
clc
syms gx gy gz mx my mz q0 q1 q2 q3
% Restringir las variables simbólicas a ser flotantes
assume([gx, gy, gz, mx, my, mz, q0, q1, q2, q3], 'real');

% Definir los campos gravitacional y magnético en el sistema inercial
g_I = [gx; gy; gz];
m_I = [mx; my; mz];
%g_I = [0; 0; -1];
%m_I = [0; 1; 0];

% Definir el cuaternión de actitud
q = [q0; q1; q2; q3];
q_conj = [q0; -q1; -q2; -q3];

% Función para convertir el cuaternión en una matriz de rotación
R = my_quat2rot(q_conj);

% Medición del campo gravitacional en el sistema del satélite
g_B = R * g_I;

% Medición del campo magnético en el sistema del satélite (debe rotarse con R)
m_B = R * m_I;

% Vector de mediciones (gravedad y magnético)
y = [g_B; m_B];

% Calcular el Jacobiano de las mediciones con respecto a los cuaterniones
J = jacobian(y, q);  % Ahora pasamos 'q' directamente como un vector simbólico
disp('Jacobiano:');
disp(J);

%subs_J = subs(J, {gx, gy, gz, mx, my, mz}, {0, 0, -1, 0, 1, 0});
subs_J = subs(J, {q0, q1, q2,q3, gx, gy, gz, mx, my, mz}, {1, 0, 0, 0, 0, 0, -1, 0, 1, 0});

% Mostrar el Jacobiano después de la sustitución
disp('Jacobiano después de la sustitución:');
disp(subs_J);

J1=measurement_jacobian([1;0;0;0], [0; 0; -1], [0; 1; 0]);

disp(J1);



function R = my_quat2rot(q)
    % Normalizar el cuaternión primero
    q = q / norm(q);

    % Extraer componentes
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    % Construir matriz de rotación
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end


function C = measurement_jacobian(q, g_I, m_I)
% Measurement Jacobian matrix for EKF
% Inputs:
%   q - quaternion [q0, q1, q2, q3] (scalar first)
%   g_I - gravity vector in inertial frame [gx, gy, gz]
%   m_I - magnetic field vector in inertial frame [mx, my, mz]
%
% Output:
%   C - 6x7 measurement Jacobian matrix

% Extract quaternion components
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

% Extract gravity components
gx = g_I(1); gy = g_I(2); gz = g_I(3);

% Extract magnetic field components
rx = m_I(1); ry = m_I(2); rz = m_I(3);

% Create the Jacobian matrix
C = 2*[gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, -gx*q2 + gy*q1 - gz*q0, -gx*q3 + gy*q0 + gz*q1, 0, 0, 0;
    -gx*q3 + gy*q0 + gz*q1,  gx*q2 - gy*q1 + gz*q0,  gx*q1 + gy*q2 + gz*q3, -gx*q0 - gy*q3 + gz*q2, 0, 0, 0;
     gx*q2 - gy*q1 + gz*q0,  gx*q3 - gy*q0 - gz*q1,  gx*q0 + gy*q3 - gz*q2,  gx*q1 + gy*q2 + gz*q3, 0, 0, 0;
     rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, -rx*q2 + ry*q1 - rz*q0, -rx*q3 + ry*q0 + rz*q1, 0, 0, 0;
    -rx*q3 + ry*q0 + rz*q1,  rx*q2 - ry*q1 + rz*q0,  rx*q1 + ry*q2 + rz*q3, -rx*q0 - ry*q3 + rz*q2, 0, 0, 0;
     rx*q2 - ry*q1 + rz*q0,  rx*q3 - ry*q0 - rz*q1,  rx*q0 + ry*q3 - rz*q2,  rx*q1 + ry*q2 + rz*q3, 0, 0, 0];
end