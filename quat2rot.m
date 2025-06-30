function R = quat2rot(q)
    % quat2rot converts a quaternion to a corresponding rotation matrix.
    % 
    % Input:
    %   q - A quaternion represented as a 4x1 vector [q0, q1, q2, q3], where:
    %       q0 is the scalar part, and
    %       q1, q2, q3 are the vector components.
    %
    % Output:
    %   R - A 3x3 rotation matrix corresponding to the input quaternion.
    %
    % Steps:
    %   1. The quaternion is normalized to ensure it is a unit quaternion.
    %   2. Extracts the scalar and vector components of the quaternion.
    %   3. Constructs the 3x3 rotation matrix using the formula for quaternion-to-rotation conversion.
    
    % Normalize the quaternion first
    q = q / norm(q);
    
    % Extract components
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    % Construct the rotation matrix
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end