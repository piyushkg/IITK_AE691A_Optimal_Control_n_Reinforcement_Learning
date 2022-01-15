function plot_cov_elps(mean,covariance,col, lw)
[evec,eval]=eig(covariance);
min_eig = min(max(eval)); 
max_eig = max(max(eval));
theta = linspace(0,2*pi);
% Ellipse parameter
a = sqrt(min_eig); aa = a*cos(theta);
b = sqrt(max_eig); bb = b*sin(theta);
Ell2 = [aa;bb];
bias = (evec')*Ell2;
Ell(1,:) = bias(1,:)+mean(1); Ell(2,:) = bias(2,:)+mean(2);
plot(Ell(1,:),Ell(2,:),col,'LineWidth',lw);
grid on;