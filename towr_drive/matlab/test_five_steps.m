close all
clear
clc

% Build FiveSteps terrain

x = 0:0.025:5.0;
h = zeros(size(x));
for i = 1:length(x)
    h(i) = five_steps(x(i));
end

figure()
% hold on
plot(x,h); grid on; axis equal;

function h = five_steps(x)

height_ref  = 0.0;
step_start  = 1.0;
step_width  = 0.1;
step_height = 0.2;
dist_steps	= 0.5; % 0.4;
slope = 2;
num_steps = 5;

h = 0.0;

for i = 0:(num_steps-1) 
    if (x >= (step_start + i*dist_steps))
      h = i*step_height + slope*(x-(step_start + i*dist_steps));
    end
    if (x >= (step_start + i*dist_steps + step_width))
      h = (i+1)*step_height;
    end
end

h = h + height_ref;

end