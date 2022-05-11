function [pos] = f_spline_fit(position, z)
    pos(1,:) = spline(z, position(1,:),-1:0.000001:1);
    pos(2,:) = spline(z, position(2,:),-1:0.000001:1);
    pos(3,:) = spline(z, position(3,:),-1:0.000001:1);
end

