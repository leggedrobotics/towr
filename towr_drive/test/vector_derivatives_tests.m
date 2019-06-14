% Testes da derivada

% clc; clear;
% 
% format long 
% 
% syms x y real
% 
% h_offset = 0.2;
% amp = 0.2;
% freq = 2;
% sine_start = 0.5;
% n_cycles = 2;
% sine_end = n_cycles*2*pi/freq + sine_start;
% 
% %% numeric
% disp('Numeric:')
% 
% h = @(x,y) amp*sin(freq*(x-sine_start)) + h_offset; %+ cos(y); 
% dhdx = @(x,y) amp*freq*cos(freq*(x-sine_start));
% dhdy = @(x,y) 0.0; %-sin(y);
% dhdxx = @(x,y) -freq*amp*freq*sin(freq*(x-sine_start));
% dhdxy = @(x,y) 0.0;
% dhdyy = @(x,y) 0.0; %-cos(y);
% dhdyx = @(x,y) 0.0;
% 
% v = @(x,y) [dhdx(x,y); dhdy(x,y); 1.0];
% 
% % example position
% p = [0.34, -0.21];
% 
% v = v(p(1),p(2));
% n = v/norm(v);
% 
% I = eye(3);
% 
% % dn/dx
% dvdx = [dhdxx(p(1),p(2)); dhdyx(p(1),p(2)); 0.0];
% dn = (norm(v) * I(:,1) - v(1)*n)/norm(v)^2;
% dndx1 = dvdx.*dn
% 
% % dn/dy
% dvdy = [dhdxy(p(1),p(2)); dhdyy(p(1),p(2)); 0.0];
% dn = (norm(v) * I(:,2) - v(2)*n)/norm(v)^2;
% dndy1 = dvdy.*dn
% 
% %% towr (after correction)
% 
% disp('Numeric towr:')
% 
% dn = ((eye(3)*norm(v) - v*v'/norm(v))/norm(v)^2)
% dndx2 = dn * dvdx
% 
% dndy2 = dn * dvdy
% 
% %% symbolic
% 
% disp('Symbolic:')
% 
% h = amp*sin(freq*(x-sine_start))+h_offset;%+cos(y);
% dhdx = amp*freq*cos(freq*(x-sine_start));
% dhdy = 0.0; %-sin(y);
% dhdxx = -freq*amp*freq*sin(freq*(x-sine_start));
% dhdxy = 0.0;
% dhdyy = 0.0; %-cos(y);
% dhdyx = 0.0;
% 
% v = [dhdx; dhdy; 1.0];
% dvdx = [dhdxx; dhdyx; 0.0];
% dvdy = [dhdxy; dhdyy; 0.0];
% 
% n = v/sqrt(v'*v);
% 
% dndx = diff(n,x);
% dndy = diff(n,y);
% 
% dndx = double(subs(dndx,{x,y},{p(1),p(2)}));
% dndy = double(subs(dndy,{x,y},{p(1),p(2)}));
% 
% dndx = (dvdx*norm(v) - v*(v'*dvdx)/norm(v))/norm(v)^2;
% dndx3 = double(subs(dndx,{x,y},{p(1),p(2)}))
% 
% dndy = (dvdy*norm(v) - v*(v'*dvdy)/norm(v))/norm(v)^2;
% dndy3 = double(subs(dndy,{x,y},{p(1),p(2)}))
% 
% %% Checks
% 
% dndx2 - dndx3 < 1e-10
% dndx1 - dndx3 < 1e-10

%% 

clear; clc;

syms x y z real

% derivative of a projection matriz wrt vector index
v = [x; y; z]; 
P = eye(3) - v*v'/(v'*v);
dPdx_sym = diff(P,x);
dPdx_sym = simplify(dPdx_sym);

dvdx = [1; 0; 0];

dPdx = -((dvdx*v'+v*dvdx')*(v'*v) - (v*v')*2*(v'*dvdx))/((v'*v)^2);
dPdx = simplify(dPdx);

% derivative of tipover axis

syms eex0 eey0 eez0 real
syms eex1 eey1 eez1 real
syms rx ry rz real
syms ax ay az real
ee0 = [eex0; eey0; eez0];
ee1 = [eex1; eey1; eez1];
com = [rx; ry; rz];
acc = [ax; ay; az];
m = 20;
a = ee1 - ee0;  % RF - LF
P = eye(3) - a*a'/(a'*a);
f = P*(m*acc);
l = P*(ee1 - com);
ln = l/sqrt(l'*l)
fn = f/sqrt(f'*f)

% checks

% wrt base pos (com_pos_)
dldnorm = (eye(3)*sqrt(l'*l) - ((l*l')/sqrt(l'*l))) / (l'*l);

sym = simplify(diff(ln,rx))

% num = simplify(dldnorm*(P*[-1; 0; 0]))
num = simplify(dldnorm*(P*(-eye(3))))

sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,rx,ry,rz},{2+0.35,0.21,0,2+0.35,-0.21,0,2,0,0.4}))
num = double(subs(num,{eex0,eey0,eez0,eex1,eey1,eez1,rx,ry,rz},{2+0.35,0.21,0,2+0.35,-0.21,0,2,0,0.4}))

fex = [60;0;40];
fex'*num

% wrt base acc (com_acc_)
dfdnorm = (eye(3)*sqrt(f'*f) - ((f*f')/sqrt(f'*f))) / (f'*f);

sym = simplify(diff(fn,ax))

num = simplify(dfdnorm*(m*P*[1; 0; 0]))
% num = simplify(dfdnorm*(m*P*(eye(3))))

sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az},{2+0.35,0.21,0,2+0.35,-0.21,0,3,0,5}))
num = double(subs(num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az},{2+0.35,0.21,0,2+0.35,-0.21,0,3,0,5}))

% wrt ee0
dldnorm = (eye(3)*sqrt(l'*l) - ((l*l')/sqrt(l'*l))) / (l'*l);
dfdnorm = (eye(3)*sqrt(f'*f) - ((f*f')/sqrt(f'*f))) / (f'*f);

% dadeex0 = [-1; 0; 0];
dadeex0 = [0; 0; -1];
dPdeex0 = -((dadeex0*a'+a*dadeex0')*(a'*a) - (a*a')*2*(a'*dadeex0))/((a'*a)^2);

sym = diff(fn,eex0);
df_num = dfdnorm * dPdeex0 * (m*acc);

sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))
df_num = double(subs(df_num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

sym = diff(ln,eex0);
dl_num = dldnorm * dPdeex0 * (ee1 - com);

sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))
dl_num = double(subs(dl_num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

stab = fn'*ln;
sym = diff(stab,eez0);
sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

num = df_num'*ln + fn'*dl_num;
num = double(subs(num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

% wrt ee1
dldnorm = (eye(3)*sqrt(l'*l) - ((l*l')/sqrt(l'*l))) / (l'*l);
dfdnorm = (eye(3)*sqrt(f'*f) - ((f*f')/sqrt(f'*f))) / (f'*f);

dadeex1 = [0; 0; 1];
% dadeex1 = [0; 0; 1];
dPdeex1 = -((dadeex1*a'+a*dadeex1')*(a'*a) - (a*a')*2*(a'*dadeex1))/((a'*a)^2);

sym = diff(fn,eex1);
df_num = dfdnorm * dPdeex1 * (m*acc);

sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))
df_num = double(subs(df_num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

sym = diff(ln,eex1);
dl_num = dldnorm * (dPdeex1*(ee1-com) + P*dadeex1);

sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))
dl_num = double(subs(dl_num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

stab = fn'*ln;
sym = diff(stab,eez1);
sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

num = df_num'*ln + fn'*dl_num;
num = double(subs(num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

% total derivative
stab = acos(fn'*ln);

sym = diff(stab,ax);
sym = double(subs(sym,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))

dldnorm = (eye(3)*sqrt(l'*l) - ((l*l')/sqrt(l'*l))) / (l'*l);
dfdnorm = (eye(3)*sqrt(f'*f) - ((f*f')/sqrt(f'*f))) / (f'*f);
dadeex1 = [-1; 0; 0];
dPdeex1 = -((dadeex1*a'+a*dadeex1')*(a'*a) - (a*a')*2*(a'*dadeex1))/((a'*a)^2);
df_num = dfdnorm * dPdeex1 * (m*acc);
dl_num = dldnorm * (dPdeex1*(ee1-com) + P*dadeex1);
num = (-1.0/(sqrt(1-(fn'*ln)^2)))*(df_num'*ln + fn'*dl_num);
num = double(subs(num,{eex0,eey0,eez0,eex1,eey1,eez1,ax,ay,az,rx,ry,rz},{3+0.35,0.21,0,2+0.35,-0.21,0,3,0,5,2,0,0.4}))
