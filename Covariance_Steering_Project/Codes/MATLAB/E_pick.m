function out = E_pick(m)
global dim_xfull
global dim_x
temp = [];
for i=1:m-1
    temp = [temp,zeros(dim_x)];
end
temp = [temp,eye(dim_x)];
temp = [temp,zeros(dim_x, dim_xfull-length(temp(1,:)))];
out = temp;
end
