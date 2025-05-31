function [xx_tt,yy_tt] = theta_calculate5(xx, yy, a, b, l)
%theta_calculate5 在盘入、盘出曲线上知道xx、yy点的坐标，求下一个把手的空间坐标。
% 这个函数主要是用来在计算前把手在盘入、盘出曲线，以及确定该把手也在盘入、盘出曲线上是，快速计算该把手坐标的函数
% xx,yy是前把手的横坐标和纵坐标
% a、b为螺旋参数
% l为该把手与前一个把手的固定距离

% xx = (a + b*theta) * cos(theta)
% yy = (a + b*theta) * sin(theta)
% xx^2 + yy^2 = (a + b*theta)^2 (a = 0)
% theta = 1/b * sqrt(xx^2 + yy^2);

% 要求是xx、yy一定要在曲线上
% xx、yy一起确定的曲线的theta
theta = 1/b * sqrt(xx^2 + yy^2); 

% l^2 = (xx - xx_tt)^2 + (yy - yy_tt)^2
% 其中xx_tt = (a + b*theta) * cos(theta)，yy_tt = (a + b*theta) * cos(theta)
% 带入参数方程后，令theta = x
fun = @(x) xx^2 + yy^2 - 2 * xx * (a + b * x) * cos(x) - 2 * yy * (a + b * x) * sin(x) + (a + b * x)^2 - l^2;

% 扫描区间
span = 5;
x_scan = linspace(theta, theta + span, 200);
y_scan = arrayfun(fun, x_scan);


% 查找第一个变号区间
found = false;
for i = 1:length(x_scan)-1
    if y_scan(i) * y_scan(i+1) < 0
        x_left = x_scan(i);
        x_right = x_scan(i+1);
        found = true;
        break;
    end
end

if ~found
    error('未找到变号点，函数可能无解或区间太小。');
end

% 使用 fzero 在区间中求零点
x = fzero(fun, [x_left, x_right]);

% xx_tt和yy_tt的参数方程
xx_tt = (a + b * x) * cos(x);
yy_tt = (a + b * x) * sin(x);

end