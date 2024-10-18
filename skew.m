% AUTHOR: 
%  LS2125204: Brayan Espinoza
% DESCRIPTION: 
% Perform the skew symetric. Needed for Matlab 2017. 
function [skewsim] = skew(vector)
skewsim=[0,-vector(3),vector(2);
        vector(3),0,-vector(1);
        -vector(2),vector(1),0];
end

