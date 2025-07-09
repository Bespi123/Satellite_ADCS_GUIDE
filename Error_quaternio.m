function [q_err] = Error_quaternio(qd, q)
%ERROR_QUATERNIO Calculates the error quaternion between two quaternions.
%
%   This function computes the error quaternion (dq) that represents the
%   rotation needed to go from the current attitude 'q' to the desired
%   attitude 'qd'. The error quaternion is defined as dq = qd * q_conj,
%   where q_conj is the conjugate of q.
%
%   The calculation is performed using matrix multiplication for efficiency,
%   avoiding explicit quaternion multiplication.
%
%   Syntax:
%       q_err = Error_quaternio(qd, q)
%
%   Inputs:
%       qd - 4x1 vector representing the **desired** attitude quaternion
%            in the format [q0; q1; q2; q3], where q0 is the scalar part.
%       q  - 4x1 vector representing the **current** attitude quaternion
%            in the same format [q0; q1; q2; q3].
%
%   Outputs:
%       q_err - 4x1 vector representing the resulting error quaternion.
%
%   Author:
%       Bespi123

% --- Input Validation ---
if ~isvector(qd) || numel(qd) ~= 4 || ~isvector(q) || numel(q) ~= 4
    error('Inputs must be 4-element vectors.');
end

% --- Quaternion Error Calculation ---
% To calculate the error quaternion q_err = qd * conjugate(q), we can
% use a matrix representation of the left-multiplication by qd.
% The conjugate of q is [q0; -q1; -q2; -q3].

% Construct the Xi matrix for left quaternion multiplication by qd.
Xi = [-qd(2), -qd(3), -qd(4);
       qd(1), -qd(4),  qd(3);
       qd(4),  qd(1), -qd(2);
      -qd(3),  qd(2),  qd(1)];

% The error quaternion's vector part (q13_err) is calculated by
% multiplying the transpose of Xi with the vector part of the current
% quaternion. This is an efficient way to perform part of the
% quaternion product. Note that this is equivalent to the vector part
% of the product between qd and the conjugate of q.
q13_err = Xi' * q(2:4);

% The scalar part of the error quaternion (q0_err) is the dot product
% of the desired and current quaternions.
q0_err = qd' * q;

% Assemble the final error quaternion.
q_err = [q0_err; q13_err];

end

