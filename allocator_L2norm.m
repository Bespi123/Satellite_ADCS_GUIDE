function T_w_L2norm = allocator_L2norm(W,Tt)
%ALLOCATOR_L2NORM Calculates an allocator based on the L2-norm (Least Squares) method.
%
%   T_w_L2norm = ALLOCATOR_L2NORM(W, Tt) computes the allocation vector
%   T_w_L2norm by solving a least squares problem using the Moore-Penrose
%   pseudoinverse. This method is commonly used in control allocation
%   problems where a desired generalized force/torque (Tt) needs to be
%   distributed among a set of actuators or effectors, represented by the
%   effectiveness matrix W.
%
%   Inputs:
%       W  - Effectiveness matrix (m x n), where 'm' is the dimension of
%            the generalized forces/torques and 'n' is the number of
%            actuators/effectors. Each column of W represents the
%            contribution of an actuator to the total generalized force/torque.
%       Tt - Desired generalized force/torque vector (m x 1). This vector
%            represents the total force or torque that needs to be generated.
%
%   Output:
%       T_w_L2norm - Allocated actuator commands/forces (n x 1). This vector
%                    contains the individual commands or forces for each
%                    actuator/effector that, when combined through W,
%                    minimize the L2-norm of the error between the desired
%                    and achieved generalized forces/torques.
%
%   Method:
%       The function solves the problem W * T_w_L2norm = Tt in a least
%       squares sense by utilizing the Moore-Penrose pseudoinverse (pinv)
%       of W. The solution is given by:
%           T_w_L2norm = pinv(W) * Tt
%
%   See also PINV, LQR, QPALLOCATOR.

    % Calculate the pseudoinverse matrix for W using the Moore-Penrose Pseudoinverse.
    % This handles cases where W is not square or is singular, providing the
    % least-squares solution.
    W_pseudo = pinv(W);

    % Calculate the L2-norm method solution for the actuator allocation.
    % This step effectively determines the actuator commands that minimize
    % the squared difference between the desired force/torque (Tt) and
    % what the actuators can produce (W * T_w_L2norm).
    T_w_L2norm = W_pseudo * Tt;

end
