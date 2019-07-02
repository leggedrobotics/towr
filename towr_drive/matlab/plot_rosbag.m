%% Plot towr rosbags in matlab
% Author: Vivian S. Medeiros
% 11/04/2019

clc;
clear;
close all;

%% Extract the desired 3D vectors from the bag
filePath = '../bags/anymal_wheels_matlab.bag';
% filePath = '../bags/multiple_slopes_bear_matlab.bag';
% filePath = '../bags/step_bear_v2_matlab.bag';
% filePath = '/home/vivian/.ros/towr_trajectory.bag';
bag_all = rosbag(filePath);

t0 = bag_all.StartTime;
T  = bag_all.EndTime;

selectOptions = {'Time', [t0 T] };
bag = select(bag_all, selectOptions{:});

% base motion
bag_base_pose = select(bag, 'Topic', 'base_pose');
ts_base_pos = timeseries(bag_base_pose);

bag_base_vel  = select(bag, 'Topic', 'base_vel_lin');
ts_base_vel_lin = timeseries(bag_base_vel);

bag_base_vel  = select(bag, 'Topic', 'base_vel_ang');
ts_base_vel_ang = timeseries(bag_base_vel);

bag_base_acc  = select(bag, 'Topic', 'base_acc_lin');
ts_base_acc_lin = timeseries(bag_base_acc);

bag_base_acc  = select(bag, 'Topic', 'base_acc_ang');
ts_base_acc_ang = timeseries(bag_base_acc);

% endeffector motion
bag_foot = select(bag, 'Topic', 'foot_pos_0');
ts_foot_LF = timeseries(bag_foot);
bag_foot_vel = select(bag, 'Topic', 'foot_vel_0');
ts_vel_LF = timeseries(bag_foot_vel);
bag_foot_acc = select(bag, 'Topic', 'foot_acc_0');
ts_acc_LF = timeseries(bag_foot_acc);

bag_foot = select(bag, 'Topic', 'foot_pos_1');
ts_foot_RF = timeseries(bag_foot);
bag_foot_vel = select(bag, 'Topic', 'foot_vel_1');
ts_vel_RF = timeseries(bag_foot_vel);
bag_foot_acc = select(bag, 'Topic', 'foot_acc_1');
ts_acc_RF = timeseries(bag_foot_acc);

bag_foot = select(bag, 'Topic', 'foot_pos_2');
ts_foot_LH = timeseries(bag_foot);
bag_foot_vel = select(bag, 'Topic', 'foot_vel_2');
ts_vel_LH = timeseries(bag_foot_vel);
bag_foot_acc = select(bag, 'Topic', 'foot_acc_2');
ts_acc_LH = timeseries(bag_foot_acc);

bag_foot = select(bag, 'Topic', 'foot_pos_3');
ts_foot_RH = timeseries(bag_foot);
bag_foot_vel = select(bag, 'Topic', 'foot_vel_3');
ts_vel_RH = timeseries(bag_foot_vel);
bag_foot_acc = select(bag, 'Topic', 'foot_acc_3');
ts_acc_RH = timeseries(bag_foot_acc);

% endeffector forces
bag_force = select(bag, 'Topic', 'foot_force_0');
ts_force_LF = timeseries(bag_force);

bag_force = select(bag, 'Topic', 'foot_force_1');
ts_force_RF  = timeseries(bag_force);

bag_force = select(bag, 'Topic', 'foot_force_2');
ts_force_LH  = timeseries(bag_force);

bag_force = select(bag, 'Topic', 'foot_force_3');
ts_force_RH  = timeseries(bag_force);


%% define the plotting range and other additional quantities
t = ts_base_pos.Time; 
dyn_con_idx = rem(t,0.1) < 1e-5;
kyn_con_idx = rem(t,0.08) < 1e-5;

% base motion
base_pos  = [ts_base_pos.Data(:,1), ts_base_pos.Data(:,2), ts_base_pos.Data(:,3)];
base_quat = [ts_base_pos.Data(:,4), ts_base_pos.Data(:,5), ...
             ts_base_pos.Data(:,6), ts_base_pos.Data(:,7)];
[base_yaw, base_pitch, base_roll] = quat2angle(base_quat,'XYZ');

% base velocity
base_vel_lin = [ts_base_vel_lin.Data(:,1), ts_base_vel_lin.Data(:,2), ts_base_vel_lin.Data(:,3)];
base_vel_ang = [ts_base_vel_ang.Data(:,1), ts_base_vel_ang.Data(:,2), ts_base_vel_ang.Data(:,3)];

% base acceleration
base_acc_lin = [ts_base_acc_lin.Data(:,1), ts_base_acc_lin.Data(:,2), ts_base_acc_lin.Data(:,3)];
base_acc_ang = [ts_base_acc_ang.Data(:,1), ts_base_acc_ang.Data(:,2), ts_base_acc_ang.Data(:,3)];

% foot motion
pos_LF = [ts_foot_LF.Data(:,1), ts_foot_LF.Data(:,2), ts_foot_LF.Data(:,3)];
pos_RF = [ts_foot_RF.Data(:,1), ts_foot_RF.Data(:,2), ts_foot_RF.Data(:,3)];
pos_LH = [ts_foot_LH.Data(:,1), ts_foot_LH.Data(:,2), ts_foot_LH.Data(:,3)];
pos_RH = [ts_foot_RH.Data(:,1), ts_foot_RH.Data(:,2), ts_foot_RH.Data(:,3)];

% foot velocity
vel_LF = [ts_vel_LF.Data(:,1), ts_vel_LF.Data(:,2), ts_vel_LF.Data(:,3)];
vel_RF = [ts_vel_RF.Data(:,1), ts_vel_RF.Data(:,2), ts_vel_RF.Data(:,3)];
vel_LH = [ts_vel_LH.Data(:,1), ts_vel_LH.Data(:,2), ts_vel_LH.Data(:,3)];
vel_RH = [ts_vel_RH.Data(:,1), ts_vel_RH.Data(:,2), ts_vel_RH.Data(:,3)];

% foot acceleration
acc_LF = [ts_acc_LF.Data(:,1), ts_acc_LF.Data(:,2), ts_acc_LF.Data(:,3)];
acc_RF = [ts_acc_RF.Data(:,1), ts_acc_RF.Data(:,2), ts_acc_RF.Data(:,3)];
acc_LH = [ts_acc_LH.Data(:,1), ts_acc_LH.Data(:,2), ts_acc_LH.Data(:,3)];
acc_RH = [ts_acc_RH.Data(:,1), ts_acc_RH.Data(:,2), ts_acc_RH.Data(:,3)];

% foot force
force_LF = [ts_force_LF.Data(:,1), ts_force_LF.Data(:,2), ts_force_LF.Data(:,3)];
force_RF = [ts_force_RF.Data(:,1), ts_force_RF.Data(:,2), ts_force_RF.Data(:,3)];
force_LH = [ts_force_LH.Data(:,1), ts_force_LH.Data(:,2), ts_force_LH.Data(:,3)];
force_RH = [ts_force_RH.Data(:,1), ts_force_RH.Data(:,2), ts_force_RH.Data(:,3)];

%% Check dynamics

Ixx =  0.248057547486776;
Iyy =  0.650151658461464;
Izz =  0.620944203624185;
Ixy =  0.001097838800893; 
Ixz = -0.003945011648535;
Iyz = -0.002135691054868;

m = 19.642;
J = [Ixx Ixy Ixz; 
     Ixy Iyy Iyz;
     Ixz Iyz Izz];  % inertia
n_ee = 4;
g = [0; 0; -9.81]; % gravity acceleration

% linear
acc_lin_dyn = (1/m)*(force_LF' + force_RF' + force_LH' + force_RH') + g;

figure()
subplot(3,1,1)
plot(t,base_acc_lin(:,1)); %hold on; plot(t,acc_lin_dyn(1,:),'k--');
hold on; plot(t(dyn_con_idx),base_acc_lin(dyn_con_idx,1),'ro');
grid on; xlabel('t [s]'); ylabel('a_x [m/s^2]')
title('Dynamics (linear)')
subplot(3,1,2)
plot(t,base_acc_lin(:,2)); %hold on; plot(t,acc_lin_dyn(2,:),'k--'); 
hold on; plot(t(dyn_con_idx),base_acc_lin(dyn_con_idx,2),'ro');
grid on; xlabel('t [s]'); ylabel('a_y [m/s^2]')
subplot(3,1,3)
plot(t,base_acc_lin(:,3)); %hold on; plot(t,acc_lin_dyn(3,:),'k--'); 
hold on; plot(t(dyn_con_idx),base_acc_lin(dyn_con_idx,3),'ro'); grid on;
% xlabel('t [s]'); ylabel('a_z [m/s^2]'); legend('from opt','from forces','constraint enforced points')

% angular
n = size(base_acc_ang,1);
forces_ee(:,:,1) = force_LF;
forces_ee(:,:,2) = force_RF;
forces_ee(:,:,3) = force_LH;
forces_ee(:,:,4) = force_RH;
pos_ee(:,:,1) = pos_LF;
pos_ee(:,:,2) = pos_RF;
pos_ee(:,:,3) = pos_LH;
pos_ee(:,:,4) = pos_RH;

acc_ang_dyn = zeros(size(base_acc_ang));
for i = 1:n
    sum = zeros(1,3);
    w_R_b = quat2rotm(base_quat(i,:));
    J_b = w_R_b*J*w_R_b;
    for j = 1:n_ee
        sum = sum + cross(forces_ee(i,:,j),base_pos(i,:)-pos_ee(i,:,j));
    end
    dyn = inv(J_b)*(sum - cross(base_vel_ang(i,:),J_b*(base_vel_ang(i,:)')))';
    acc_ang_dyn(i,:) = dyn';
end

figure()
subplot(3,1,1)
plot(t,base_acc_ang(:,1)); hold on; plot(t,acc_ang_dyn(:,1),'k--');
hold on; plot(t(dyn_con_idx),base_acc_ang(dyn_con_idx,1),'ro');
grid on; xlabel('t [s]'); ylabel('$$\dot{\omega}_x [rad/s^2]$$','Interpreter','Latex');
title('Dynamics (angular)')
subplot(3,1,2)
plot(t,base_acc_ang(:,2)); hold on; plot(t,acc_ang_dyn(:,2),'k--'); 
hold on; plot(t(dyn_con_idx),base_acc_ang(dyn_con_idx,2),'ro');
grid on; xlabel('t [s]'); ylabel('$$\dot{\omega}_y [rad/s^2]$$','Interpreter','Latex');
subplot(3,1,3)
plot(t,base_acc_ang(:,3)); hold on; plot(t,acc_ang_dyn(:,3),'k--'); 
hold on; plot(t(dyn_con_idx),base_acc_ang(dyn_con_idx,3),'ro');
grid on; xlabel('t [s]'); ylabel('$$\dot{\omega}_z [rad/s^2]$$','Interpreter','Latex'); 
legend('from opt','from forces','constraint enforced points')

%% plot the values

% base
figure()
subplot(3,2,1)
plot(t,base_pos(:,1)); grid on;
xlabel('t [s]'); ylabel('x [m]')
title('Base Linear Position')
subplot(3,2,3)
plot(t,base_pos(:,2)); grid on;
xlabel('t [s]'); ylabel('y [m]')
subplot(3,2,5)
plot(t,base_pos(:,3)); grid on;
xlabel('t [s]'); ylabel('z [m]')

subplot(3,2,2)
plot(t,base_roll*180/pi); grid on;
xlabel('t [s]'); ylabel('roll [deg]')
title('Base Angular Position')
subplot(3,2,4)
plot(t,base_pitch*180/pi); grid on;
xlabel('t [s]'); ylabel('pitch [deg]')
subplot(3,2,6)
plot(t,base_yaw*180/pi); grid on;
xlabel('t [s]'); ylabel('yaw [deg]')

figure()
subplot(3,1,1)
plot(t,base_vel_lin(:,1)); grid on;
hold on; plot(t(dyn_con_idx),base_vel_lin(dyn_con_idx,1),'ro');
xlabel('t [s]'); ylabel('v_x [m/s]')
title('Base Linear Velocity')
subplot(3,1,2)
plot(t,base_vel_lin(:,2)); grid on;
hold on; plot(t(dyn_con_idx),base_vel_lin(dyn_con_idx,2),'ro');
xlabel('t [s]'); ylabel('v_y [m/s]')
subplot(3,1,3)
plot(t,base_vel_lin(:,3)); grid on;
hold on; plot(t(dyn_con_idx),base_vel_lin(dyn_con_idx,3),'ro');
xlabel('t [s]'); ylabel('v_z [m/s]')

% subplot(3,2,2)
% plot(t,base_vel_ang(:,1)); grid on;
% xlabel('t [s]'); ylabel('\omega_x [rad/s]')
% title('Base Angular Velocity')
% subplot(3,2,4)
% plot(t,base_vel_ang(:,2)); grid on;
% xlabel('t [s]'); ylabel('\omega_y [rad/s]')
% subplot(3,2,6)
% plot(t,base_vel_ang(:,3)); grid on;
% xlabel('t [s]'); ylabel('\omega_z [rad/s]')

% foot motion
h = figure();
set(h, 'Name', 'Wheels position');
subplot(3,4,1); plot(t,pos_LF(:,1)); grid on;
xlabel('t [s]'); ylabel('p_x [m]'); title('LF')
subplot(3,4,5); plot(t,pos_LF(:,2)); grid on;
xlabel('t [s]'); ylabel('p_y [m]'); %ylim([0.18 0.2]);
subplot(3,4,9); plot(t,pos_LF(:,3)); grid on;
xlabel('t [s]'); ylabel('p_z [m]')

subplot(3,4,2); plot(t,pos_RF(:,1)); grid on;
xlabel('t [s]'); ylabel('p_x [m]'); title('RF')
subplot(3,4,6); plot(t,pos_RF(:,2)); grid on;
xlabel('t [s]'); ylabel('p_y [m]'); %ylim([-0.2 -0.18]);
subplot(3,4,10); plot(t,pos_RF(:,3)); grid on;
xlabel('t [s]'); ylabel('p_z [m]');

subplot(3,4,3); plot(t,pos_LH(:,1)); grid on;
xlabel('t [s]'); ylabel('p_x [m]'); title('LH')
subplot(3,4,7); plot(t,pos_LH(:,2)); grid on;
xlabel('t [s]'); ylabel('p_y [m]'); %ylim([0.18 0.2]);
subplot(3,4,11); plot(t,pos_LH(:,3)); grid on;
xlabel('t [s]'); ylabel('p_z [m]')

subplot(3,4,4); plot(t,pos_RH(:,1)); grid on;
xlabel('t [s]'); ylabel('p_x [m]'); title('RH')
subplot(3,4,8); plot(t,pos_RH(:,2)); grid on;
xlabel('t [s]'); ylabel('p_y [m]'); %ylim([-0.2 -0.18]);
subplot(3,4,12); plot(t,pos_RH(:,3)); grid on;
xlabel('t [s]'); ylabel('p_z [m]')

h = figure();
set(h, 'Name', 'Wheels position (X x Z)');
subplot(4,1,1); plot(pos_LF(:,1),pos_LF(:,3)); grid on; axis equal;
xlabel('p_x [m]'); ylabel('p_z [m]'); title('LF')
subplot(4,1,2); plot(pos_RF(:,1),pos_RF(:,3)); grid on; axis equal;
xlabel('p_x [m]'); ylabel('p_z [m]'); title('RF')
subplot(4,1,3); plot(pos_LH(:,1),pos_LH(:,3)); grid on; axis equal;
xlabel('p_x [m]'); ylabel('p_z [m]'); title('LH')
subplot(4,1,4); plot(pos_RH(:,1),pos_RH(:,3)); grid on; axis equal;
xlabel('p_x [m]'); ylabel('p_z [m]'); title('RH')

% foot velocity
h = figure();
set(h, 'Name', 'Wheels velocity');
subplot(3,4,1); plot(t,vel_LF(:,1)); grid on;
xlabel('t [s]'); ylabel('v_x [m/s]'); title('LF')
subplot(3,4,5); plot(t,vel_LF(:,2)); grid on;
xlabel('t [s]'); ylabel('v_y [m/s]')
subplot(3,4,9); plot(t,vel_LF(:,3)); grid on;
xlabel('t [s]'); ylabel('v_z [m/s]')

subplot(3,4,2); plot(t,vel_RF(:,1)); grid on;
xlabel('t [s]'); ylabel('v_x [m/s]'); title('RF')
subplot(3,4,6); plot(t,vel_RF(:,2)); grid on;
xlabel('t [s]'); ylabel('v_y [m/s]'); 
subplot(3,4,10); plot(t,vel_RF(:,3)); grid on;
xlabel('t [s]'); ylabel('v_z [m/s]')

subplot(3,4,3); plot(t,vel_LH(:,1)); grid on;
xlabel('t [s]'); ylabel('v_x [m/s]'); title('LH')
subplot(3,4,7); plot(t,vel_LH(:,2)); grid on;
xlabel('t [s]'); ylabel('v_y [m/s]')
subplot(3,4,11); plot(t,vel_LH(:,3)); grid on;
xlabel('t [s]'); ylabel('v_z [m/s]')

subplot(3,4,4); plot(t,vel_RH(:,1)); grid on;
xlabel('t [s]'); ylabel('v_x [m/s]'); title('RH')
subplot(3,4,8); plot(t,vel_RH(:,2)); grid on;
xlabel('t [s]'); ylabel('v_y [m/s]')
subplot(3,4,12); plot(t,vel_RH(:,3)); grid on;
xlabel('t [s]'); ylabel('v_z [m/s]')

% foot acceleration
h = figure();
set(h, 'Name', 'Wheels acceleration');
subplot(3,4,1); plot(t,acc_LF(:,1)); grid on;
xlabel('t [s]'); ylabel('a_x [m/s^2]'); title('LF')
subplot(3,4,5); plot(t,acc_LF(:,2)); grid on;
xlabel('t [s]'); ylabel('a_y [m/s^2]')
subplot(3,4,9); plot(t,acc_LF(:,3)); grid on;
xlabel('t [s]'); ylabel('a_z [m/s^2]')

subplot(3,4,2); plot(t,acc_RF(:,1)); grid on;
xlabel('t [s]'); ylabel('a_x [m/s^2]'); title('RF')
subplot(3,4,6); plot(t,acc_RF(:,2)); grid on;
xlabel('t [s]'); ylabel('a_y [m/s^2]')
subplot(3,4,10); plot(t,acc_RF(:,3)); grid on;
xlabel('t [s]'); ylabel('a_z [m/s^2]')

subplot(3,4,3); plot(t,acc_LH(:,1)); grid on;
xlabel('t [s]'); ylabel('a_x [m/s^2]'); title('LH')
subplot(3,4,7); plot(t,acc_LH(:,2)); grid on;
xlabel('t [s]'); ylabel('a_y [m/s^2]')
subplot(3,4,11); plot(t,acc_LH(:,3)); grid on;
xlabel('t [s]'); ylabel('a_z [m/s^2]')

subplot(3,4,4); plot(t,acc_RH(:,1)); grid on;
xlabel('t [s]'); ylabel('a_x [m/s^2]'); title('RH')
subplot(3,4,8); plot(t,acc_RH(:,2)); grid on;
xlabel('t [s]'); ylabel('a_y [m/s^2]')
subplot(3,4,12); plot(t,acc_RH(:,3)); grid on;
xlabel('t [s]'); ylabel('a_z [m/s^2]')

% foot forces
h = figure();
set(h, 'Name', 'Wheels force');
subplot(3,4,1); plot(t,force_LF(:,1)); grid on;
xlabel('t [s]'); ylabel('f_x [N]'); title('LF')
subplot(3,4,5); plot(t,force_LF(:,2)); grid on;
xlabel('t [s]'); ylabel('f_y [N]')
subplot(3,4,9); plot(t,force_LF(:,3)); grid on;
xlabel('t [s]'); ylabel('f_z [N]')

subplot(3,4,2); plot(t,force_RF(:,1)); grid on;
xlabel('t [s]'); ylabel('f_x [N]'); title('RF')
subplot(3,4,6); plot(t,force_RF(:,2)); grid on;
xlabel('t [s]'); ylabel('f_y [N]')
subplot(3,4,10); plot(t,force_RF(:,3)); grid on;
xlabel('t [s]'); ylabel('f_z [N]')

subplot(3,4,3); plot(t,force_LH(:,1)); grid on;
xlabel('t [s]'); ylabel('f_x [N]'); title('LH')
subplot(3,4,7); plot(t,force_LH(:,2)); grid on;
xlabel('t [s]'); ylabel('f_y [N]')
subplot(3,4,11); plot(t,force_LH(:,3)); grid on;
xlabel('t [s]'); ylabel('f_z [N]')

subplot(3,4,4); plot(t,force_RH(:,1)); grid on;
xlabel('t [s]'); ylabel('f_x [N]'); title('RH')
subplot(3,4,8); plot(t,force_RH(:,2)); grid on;
xlabel('t [s]'); ylabel('f_y [N]')
subplot(3,4,12); plot(t,force_RH(:,3)); grid on;
xlabel('t [s]'); ylabel('f_z [N]')

base_ang = [base_roll, base_pitch, base_yaw];
plot_stability_measure(base_pos, base_ang, pos_ee, base_acc_lin, base_acc_ang, base_vel_ang, J, m ,g, t);

%% compute wheel's torque

wh_radius = 0.05; % [m]
wh_torque = zeros(size(forces_ee));
for i = 1:n
    for j = 1:n_ee
        x = pos_ee(i,1,j);
        [n, tx, ty] = GetTerrainBasis(x, "Flat");
        wh_torque(i,:,j) = (forces_ee(i,:,j) * tx) * wh_radius * 2;    
    end
end

h = figure();
set(h, 'Name', 'Wheels torque');
subplot(1,4,1); plot(t,wh_torque(:,:,1),'m','LineWidth',2.0); grid on;
ylabel('torque [Nm]'); title('LF'); xlabel('t [s]')
ax = gca;
ax.FontSize = 14;
subplot(1,4,2); plot(t,wh_torque(:,:,2),'m','LineWidth',2.0); grid on;
ylabel('torque [Nm]'); title('RF'); xlabel('t [s]')
ax = gca;
ax.FontSize = 14;
subplot(1,4,3); plot(t,wh_torque(:,:,3),'m','LineWidth',2.0); grid on;
ylabel('torque [Nm]'); title('LH'); xlabel('t [s]')
ax = gca;
ax.FontSize = 14;
subplot(1,4,4); plot(t,wh_torque(:,:,4),'m','LineWidth',2.0); grid on;
ylabel('torque [Nm]'); title('RH'); xlabel('t [s]')
ax = gca;
ax.FontSize = 14;





