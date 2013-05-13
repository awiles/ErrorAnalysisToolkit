function data = convertPathForSimulink(filename)
%% convert paths used in FLE testing to time stamp format for Simulink.
load(filename);

M = length(xfrm);
data = zeros(8,M); % data rows: t, x, y, z, q0, qx, qy, qz
for i = 1:M
    t = i/60;
    q = rm2quat(xfrm{i}.R);
    data(:,i) = [t, xfrm{i}.pos, q]';
end