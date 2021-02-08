function [time, idx] = build_time_vec (vec, t)

time = zeros(size(vec));
idx = ones(size(vec));
for i = 1:length(vec)
    time(i) = sum(vec(1:i));
    tmp = find(abs(t-time(i)) < 0.001);
    id = 1;
    if (isempty(tmp))
        tmp = find(abs(t-time(i)) < 0.002);  % try with a higher tol
        if  (isempty(tmp) == 0)
            [m, id] = min(abs(t(tmp)-time(i)));
            idx(i) = tmp(id);
        else
            idx(i) = length(t);
        end
    else
        idx(i) = tmp(id);
    end
%     disp([i time(i), tmp', t(tmp)', tmp(id)])
    
end

time = [0; time];
idx = [1; idx];

end