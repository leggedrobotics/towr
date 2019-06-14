function plot_stability_measure (com_pos,com_ang, pos_ee, com_acc, omega_dot, omega, I_b, m, g, t)

ee_pos_cw = zeros(size(pos_ee,1),size(pos_ee,2),size(pos_ee,3)+1);
ee_pos_cw(:,:,1) = pos_ee(:,:,1);
ee_pos_cw(:,:,2) = pos_ee(:,:,2);
ee_pos_cw(:,:,3) = pos_ee(:,:,4);
ee_pos_cw(:,:,4) = pos_ee(:,:,3);
ee_pos_cw(:,:,5) = pos_ee(:,:,1);

theta = zeros(size(pos_ee,1), size(pos_ee,3), 2);
for i = 1:size(pos_ee,1)
    w_R_b = GetRotationMatrixBaseToWorld(com_ang(i,:));
    I_w = w_R_b * I_b * w_R_b';
    for ee = 1:size(pos_ee,3)
        a = (ee_pos_cw(i,:,ee+1) - ee_pos_cw(i,:,ee))';
        a_n = a/norm(a);
        proj_a = eye(3) - a_n*a_n';
        l = proj_a * (ee_pos_cw(i,:,ee+1) - com_pos(i,:))';
        l_n = l/norm(l);
        f = proj_a * (m*com_acc(i,:)' + m*g);
        n = (eye(3) - proj_a) * (I_w*omega_dot(i,:)' + cross(omega(i,:)', I_w*omega(i,:)'));
        ft = f + cross(l_n,n)/norm(l);
        ft_n = ft/norm(ft);
        d = -l + (l'*ft_n)*ft_n;
        tmp = cross(ft_n,l_n)' * a_n;
        theta(i,ee,1) = sign(tmp) * acos(ft_n'*l_n);  % with angular loads
        
        f_n = f/norm(f);
        tmp = cross(f_n,l_n)' * a_n;
        theta(i,ee,2) = sign(tmp) * acos(f_n'*l_n);  % without angular loads
        
        (m*com_acc(i,:)' + m*g)'*l;
    end
end

stab = [theta(:,1,1)' ; theta(:,2,1)'; theta(:,3,1)'; theta(:,4,1)'];
stab_min = min(stab);

h = figure();
set(h, 'Name', 'Stability Angle');
subplot(1,4,1); plot(t,theta(:,1,1)*180/pi,t,theta(:,1,2)*180/pi); 
grid on; xlabel('t [s]'); ylabel('\theta [°]'); title('LF-RF')
subplot(1,4,2); plot(t,theta(:,2,1)*180/pi,t,theta(:,2,2)*180/pi); 
grid on; xlabel('t [s]'); ylabel('\theta [°]'); title('RF-RH')
subplot(1,4,3); plot(t,theta(:,3,1)*180/pi,t,theta(:,3,2)*180/pi); 
grid on; xlabel('t [s]'); ylabel('\theta [°]'); title('RH-LH')
subplot(1,4,4); plot(t,theta(:,4,1)*180/pi,t,theta(:,4,2)*180/pi); 
grid on; xlabel('t [s]'); ylabel('\theta [°]'); title('LH-LF')

figure();
plot(t,stab_min*180/pi,'LineWidth',2); 
grid on; xlabel('t [s]'); ylabel('\theta_{min} [°]'); 
title('Stability angle')
ax = gca;
ax.FontSize = 14;

end