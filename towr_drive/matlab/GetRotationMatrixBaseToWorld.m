function M = GetRotationMatrixBaseToWorld (xyz)

x = xyz(1);
y = xyz(2);
z = xyz(3);

M = [ cos(y)*cos(z), cos(z)*sin(x)*sin(y) - cos(x)*sin(z), sin(x)*sin(z) + cos(x)*cos(z)*sin(y);
      cos(y)*sin(z), cos(x)*cos(z) + sin(x)*sin(y)*sin(z), cos(x)*sin(y)*sin(z) - cos(z)*sin(x);
            -sin(y),                        cos(y)*sin(x),                        cos(x)*cos(y)];
        
end