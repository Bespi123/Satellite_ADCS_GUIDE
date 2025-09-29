classdef Satellite
    % SATELLITE Models the state, physical properties, and dynamics of a satellite.
    %
    %   This class encapsulates the attitude state (quaternion and angular velocity)
    %   and the rotational dynamics of a satellite. Its primary responsibility is
    %   to propagate the state over time based on the external torques applied.
    %
    %   Properties:
    %       State   - The satellite's current state vector [q (4x1); w (3x1)].
    %       Inertia - The satellite's 3x3 inertia matrix [kg*m^2].
    %
    %   Methods:
    %       Satellite   - The constructor to create an instance of the class.
    %       updateState - Updates the satellite's state for a given time step 'dt'.
    %
    %   Usage Example:
    %       initial.q0 = [1; 0; 0; 0];
    %       initial.Wo = [0.1; -0.1; 0.05];
    %       I = diag([0.1, 0.2, 0.15]);
    %       mySatellite = Satellite(initial, I);
    %       mySatellite = mySatellite.updateState(T_ctrl, T_dist, 0.1);

    properties (SetAccess = public)
        State       % State vector [q (4x1); w (3x1)]
        Inertia     % 3x3 Inertia matrix [kg*m^2]
    end

    methods
        % --- CONSTRUCTOR ---
        function obj = Satellite(initial_conditions, inertia_matrix)
            % SATELLITE Class constructor.
            %
            %   Syntax:
            %       obj = Satellite(initial_conditions, inertia_matrix)
            %
            %   Inputs:
            %       initial_conditions - A struct with the initial conditions:
            %                            .q0 -> Initial attitude quaternion (4x1).
            %                            .Wo -> Initial angular velocity (3x1) [rad/s].
            %       inertia_matrix     - 3x3 Inertia Matrix [kg*m^2].
            %
            %   Outputs:
            %       obj - A new instance of the Satellite object.
            
            obj.State = [initial_conditions.q0; initial_conditions.Wo];
            obj.Inertia = inertia_matrix;
        end
        
        % --- MAIN PUBLIC METHOD ---
        function obj = updateState(obj, T_control, T_disturbance, dt)
            % updateState Integrates the satellite's dynamics over a single time step.
            %
            %   This method uses a 4th-order Runge-Kutta (RK4) integrator to
            %   propagate the satellite's state. This is the primary public method
            %   that should be called from the main simulation loop.
            %
            %   Syntax:
            %       obj = obj.updateState(T_control, T_disturbance, dt)
            %
            %   Inputs:
            %       obj           - The current instance of the Satellite object.
            %       T_control     - The applied control torque (3x1) [Nm].
            %       T_disturbance - The external disturbance torque (3x1) [Nm].
            %       dt            - The integration time step [s].
            %
            %   Outputs:
            %       obj - The Satellite object with its .State property updated.

            % Use the object's current state for the first RK4 step.
            g1 = dt * obj.stateDerivative(T_control, T_disturbance, obj.State);
            
            % Use temporary states for the intermediate RK4 steps.
            g2 = dt * obj.stateDerivative(T_control, T_disturbance, obj.State + 0.5 * g1);
            g3 = dt * obj.stateDerivative(T_control, T_disturbance, obj.State + 0.5 * g2);
            g4 = dt * obj.stateDerivative(T_control, T_disturbance, obj.State + g3);
            
            % Update the object's actual state using the RK4 formula.
            obj.State = obj.State + (1/6) * (g1 + 2*g2 + 2*g3 + g4);
            
            % Normalize the quaternion to prevent numerical drift and ensure
            % it remains a valid unit rotation.
            obj.State(1:4) = obj.State(1:4) / norm(obj.State(1:4));
        end
    end

    methods (Access = private)
        % --- PRIVATE METHOD ---
        function x_dot = stateDerivative(obj, U, Td, current_state)
            % stateDerivative Calculates the state derivative [q_dot; w_dot].
            %
            %   This method implements the rotational kinematics and dynamics
            %   equations for a rigid body. It is a private method because it
            %   represents the internal logic of the satellite's model and does not
            %   need to be accessed from outside the class.
            %
            %   Inputs:
            %       obj           - The current instance of the Satellite object.
            %       U             - Control torque [Nm].
            %       Td            - Disturbance torque [Nm].
            %       current_state - The state vector [q; w] to be evaluated.
            %
            %   Outputs:
            %       x_dot - The derivative of the state vector [q_dot; w_dot].
            
            % 1. Unpack the current state vector
            q = current_state(1:4);
            w = current_state(5:7);
            
            % 2. Quaternion Kinematics
            % Calculates the rate of change of the quaternion.
            Xi = [-q(2), -q(3), -q(4);
                   q(1), -q(4),  q(3);
                   q(4),  q(1), -q(2);
                  -q(3),  q(2),  q(1)];
            q_dot = 0.5 * Xi * w;
            
            % 3. Rotational Dynamics (Euler's Equation for a Rigid Body)
            % Calculates the angular acceleration.
            w_dot = obj.Inertia \ (Td + U - cross(w, obj.Inertia * w));
            
            % 4. Assemble the state derivative vector
            x_dot = [q_dot; w_dot];
        end
    end
end