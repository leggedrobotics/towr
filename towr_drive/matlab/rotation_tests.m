% Tests with the wheels coordinate frame

clc; clear;

% v_c = [0.5 0 0]';
v_w = [0.5 0 0]';

xyz = [0 0 20]*pi/180;
w_R_b = GetRotationMatrixBaseToWorld(xyz);
b_R_w = inv(w_R_b);

% v_b = b_R_w * v_w

% n = [-2 0 1]';
% n = [0 0 1]';
% n = [1 0 1]';
[n, tx, ty] = GetTerrainBasis(1.05, "Flat"); %"Step")
n = n/norm(n)
heading = w_R_b * [1 0 0]';
ey = cross(n,heading); %ey = ey / norm(ey);
ex = cross(ey,n);      %ex = ex / norm(ex);
ez = cross(ex,ey);     %ez = ez / norm(ez);
w_R_c = [ex/norm(ex) ey/norm(ey) ez/norm(ez)]
c_R_w = w_R_c';
quat = rotm2quat(w_R_c);
% v_w = w_R_c * v_c
v_c = c_R_w * v_w

% xyz = [0 -atan(2) 0];
% w_R_c = GetRotationMatrixBaseToWorld(xyz)
% v_w = w_R_c * v_c
% w_R_c' * v_w;


