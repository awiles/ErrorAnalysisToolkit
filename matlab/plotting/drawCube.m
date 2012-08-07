function drawCube ( origin, size )

if(length(size) == 1)
    sz = size*ones(1,3);
elseif(length(size) == 3)
    sz = size;
else
    error('Invalid dimensions for size variable. Length should be 1 (cube) or 3 (rectangular prism).');
end
x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*sz(1)+origin(1);
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*sz(2)+origin(2);
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*sz(3)+origin(3);
for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),'w');
    set(h,'facecolor', [0 .5 .6], 'edgecolor','k', 'FaceAlpha', 0.25)
end 