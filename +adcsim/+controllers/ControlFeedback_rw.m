function [U] = ControlFeedback_rw(I, x, dq, Wr, Wr_dot, P, K)
    import adcsim.utils.skew
    % Author: bespi123
    % Control law carried out by Lyapunov with feedback control law
    % This function implements the control law for a rigid body using feedback
    % based on the Lyapunov stability criterion to stabilize the angular
    % velocity and quaternion attitude.
    %
    % The function calculates the control torque (U) required to bring the system 
    % to a desired state based on the current state and control gains.
    %
    % Inputs:
    %   I         : Inertia tensor of the rigid body (3x3 matrix).
    %   x         : State vector containing the current state of the system.
    %             x(1:3) = position, x(4:6) = orientation (quaternion),
    %             x(5:7) = angular velocity (angular rates in body frame).
    %   dq        : Quaternion error (4x1 vector), representing the attitude error.
    %   Wr        : Desired angular velocity (3x1 vector) in body frame.
    %   Wr_dot    : Time derivative of the desired angular velocity (3x1 vector).
    %   P         : Proportional gain matrix (3x3 matrix).
    %   K         : Derivative gain matrix (3x3 matrix).
    %
    % Output:
    %   U         : Control input (torque vector, 3x1).
    
    % Calculate the angular velocity vector from the state vector 'x'
    W = x(5:7);  % Angular rate (angular velocity) in body frame from state vector
    
    % Extract quaternion error (dq)
    dq13 = [dq(2), dq(3), dq(4)]';  % The error in the quaternion components (3x1 vector)
    
    % Control law:
    % The goal is to minimize the angular velocity error (W - Wr) and quaternion error (dq13)
    % to stabilize the attitude and angular velocity of the system.
    
    dW = W - Wr;  % Angular rate error (difference between current and desired angular velocity)
    
    % Compute the control torque input (U)
    U = -P * dW - K * dq13 + skew(W) * I * W + I * (Wr_dot - skew(W) * Wr);
end