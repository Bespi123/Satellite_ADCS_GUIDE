function dx = BrushelessModel_simu(x, u, param)
%BrushelessModel_simu Models the dynamics of a brushless DC motor for simulation.
%
%   This function computes the state derivatives of a brushless DC motor,
%   making it suitable for use with an ordinary differential equation (ODE)
%   solver like ode45.
%
%   Syntax:
%       dx = BrushelessModel_simu(x, u, param)
%
%   Inputs:
%       x     - State vector (column) [2x1]:
%               x(1) = w: Motor angular velocity [rad/s]
%               x(2) = i: Motor current [A]
%
%       u     - Input voltage applied to the motor (scalar) [V].
%
%       param - Structure containing the physical parameters of the motor:
%               .kt = Torque constant [N*m/A]
%               .J  = Rotor inertia [kg*m^2]
%               .B  = Viscous friction coefficient [N*m*s]
%               .kc = Coulomb friction coefficient [N*m]
%               .L  = Inductance [H]
%               .R  = Resistance [Ω]
%               .ke = Back-EMF constant [V*s/rad]
%
%   Outputs:
%       dx    - Column vector [2x1] containing the state derivatives:
%               dx(1) = w_dot: Angular acceleration [rad/s^2]
%               dx(2) = i_dot: Rate of change of current [A/s]

    %% Extract states, parameters, and input
    % System states
    w = x(1); % [rad/s] Angular velocity
    i = x(2); % [A] Current

    % Motor parameters
    kt = param.kt; % [N*m/A] Torque constant
    J  = param.Jrw;  % [kg*m^2] Rotor inertia
    B  = param.b;  % [N*m*s] Viscous friction
    Kc = param.c; % [N*m] Coulomb friction
    L  = param.L;  % [H] Inductance
    R  = param.R;  % [Ω] Resistance
    Ke = param.ke; % [V*s/rad] Back-EMF constant

    %% Compute the state derivatives
    
    % --- Mechanical Equation ---
    % Computes the angular acceleration (w_dot) based on the sum of torques:
    % Motor torque (kt*i) - Viscous friction torque (B*w) - Coulomb friction torque (Kc*sign(w))
    w_dot = (1/J) * (kt * i - B * w - Kc * sign(w));

    % --- Electrical Equation ---
    % Computes the rate of change of current (i_dot) using Kirchhoff's Voltage Law:
    % Input voltage (u) - Resistive drop (R*i) - Back-EMF (Ke*w)
    i_dot = (1/L) * (u - R * i - Ke * w);

    %% Assemble the state derivative vector
    dx = [w_dot; i_dot];
end