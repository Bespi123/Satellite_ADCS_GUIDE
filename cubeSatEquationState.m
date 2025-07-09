function [x_dot] = cubeSatEquationState(Td, I, U, x)
%CUBESATEQUATIONSTATE Computes the state derivative for CubeSat attitude dynamics.
%
%   This function models the rotational dynamics and kinematics of a rigid
%   body, specifically a CubeSat. It uses a state-space representation
%   where the state vector 'x' consists of the attitude quaternion and the
%   angular velocity vector. The output is the time derivative of this
%   state vector, which can be used in an ODE solver like ode45.
%
%   The model is based on the following equations:
%   1.  **Kinematics**: Describes the rate of change of the attitude
%       quaternion as a function of the angular velocity.
%   2.  **Dynamics**: Euler's moment equation, which describes the rate of
%       change of the angular velocity as a function of inertia, applied
%       torques, and gyroscopic effects.
%
%   Syntax:
%       x_dot = cubeSatEquationState(Td, I, U, x)
%
%   Inputs:
%       Td - 3x1 vector of external disturbance torques [Nm].
%       I  - 3x3 inertia matrix of the rigid body [kg*m^2].
%       U  - 3x1 vector of control input torques [Nm].
%       x  - 7x1 state vector, where x = [q; w]:
%            q - 4x1 attitude quaternion [q1, q2, q3, q0].
%            w - 3x1 angular velocity vector [wx, wy, wz] in rad/s.
%
%   Outputs:
%       x_dot - 7x1 state derivative vector, where x_dot = [q_dot; w_dot].
%
%   References:
%       - Markley, F. L., & Crassidis, J. L. (2014).
%         *Fundamentals of Spacecraft Attitude Determination and Control*.
%         Springer.
%
%   Author: Bespi123

% --- 1. Unpack State Vector ---
% Extract the quaternion and angular velocity from the input state vector.
q = x(1:4);
w = x(5:7);

% --- 2. Kinematic and Dynamic Equations ---
% Based on "Fundamentals of Spacecraft Attitude Determination and Control"
% by F. Landis Markley & John L. Crassidis.

% Quaternion Kinematics - See Equation (3.21)
% First, construct the Xi matrix from the quaternion components.
Xi = [-q(2), -q(3), -q(4);
       q(1), -q(4),  q(3);
       q(4),  q(1), -q(2);
      -q(3),  q(2),  q(1)];

% Calculate the time derivative of the quaternion.
q_dot = 0.5 * Xi * w;

% Rotational Dynamics (Euler's Moment Equation) - See Equation (3.147)
% This equation calculates the angular acceleration.
w_dot = I \ (Td + U - cross(w, I * w));

% --- 3. Assemble State Derivative Vector ---
% Combine the derivatives of the quaternion and angular velocity into a
% single output vector for use with an ODE solver.
x_dot = [q_dot; w_dot];

end

