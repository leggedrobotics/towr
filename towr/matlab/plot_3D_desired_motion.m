function plot_3D_desired_motion(base_pos, pos_ee)

figure()
plot3(base_pos(:,1), base_pos(:,2), base_pos(:,3),'LineWidth',2); hold on;
plot3(pos_ee(:,1,1), pos_ee(:,2,1), pos_ee(:,3,1),'LineWidth',2); hold on;
plot3(pos_ee(:,1,2), pos_ee(:,2,2), pos_ee(:,3,2),'LineWidth',2); hold on;
plot3(pos_ee(:,1,3), pos_ee(:,2,3), pos_ee(:,3,3),'LineWidth',2); hold on;
plot3(pos_ee(:,1,4), pos_ee(:,2,4), pos_ee(:,3,4),'LineWidth',2); hold on;
grid on;
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');

end