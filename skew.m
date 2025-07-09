function [skewsim] = skew(vector)
%SKEW Generates the skew-symmetric matrix of a 3-element vector.
%
%   This function is essential for cross-product operations, as the
%   cross product of two vectors (a x b) is equivalent to the matrix-vector
%   multiplication of the skew-symmetric matrix of 'a' and the vector 'b'.
%   It is particularly useful for versions of MATLAB prior to R2022a, which
%   introduced a built-in `skew` function.
%
%   Syntax:
%       skewsim = skew(vector)
%
%   Inputs:
%       vector - A 3x1 or 1x3 vector.
%
%   Outputs:
%       skewsim - The corresponding 3x3 skew-symmetric matrix.
%
% AUTHOR: bespi123

% --- Input Validation ---
% Ensure the input is a 3-element vector for the operation to be valid.
if ~isvector(vector) || numel(vector) ~= 3
    error('Input must be a 3-element vector.');
end

% --- Matrix Construction ---
% Construct the skew-symmetric matrix from the vector components.
skewsim = [    0    , -vector(3),  vector(2);
             vector(3),     0    , -vector(1);
            -vector(2),  vector(1),     0    ];
end