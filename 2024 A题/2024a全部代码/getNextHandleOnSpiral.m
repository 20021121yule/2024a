function [xx_tt,yy_tt] = getNextHandleOnSpiral(xx, yy, a, b, l, flag_in, flag_out)
%theta_calculate5 在盘入、盘出曲线上知道xx、yy点的坐标，求下一个把手的空间坐标。
% 这个函数主要是用来在计算前把手在盘入、盘出曲线，以及确定该把手也在盘入、盘出曲线上是，快速计算该把手坐标的函数
% xx,yy是前把手的横坐标和纵坐标
% a、b为螺旋参数
% l为该把手与前一个把手的固定距离
% flag_spiral是前后点在曲线上的种类，如果输入的是flag_in(输入模式为：1, 0), x_san应该是linspace(theta, theta + span, 200);
%                                如果输入的是flag_out(输入模式为：0, 1), x_san应该是linspace(theta - span, theta, 200);

% xx = (a + b*theta) * cos(theta)
% yy = (a + b*theta) * sin(theta)
% xx^2 + yy^2 = (a + b*theta)^2 (a = 0)
% theta = 1/b * sqrt(xx^2 + yy^2);

% 要求是xx、yy一定要在曲线上
% xx、yy一起确定的曲线的theta
theta = 1/b * sqrt(xx^2 + yy^2);

% l^2 = (xx - xx_tt)^2 + (yy - yy_tt)^2
% 对于盘入曲线：其中xx_tt = (a + b*theta) * cos(theta)，yy_tt = (a + b*theta) * sin(theta)
% 对于盘出曲线：其中xx_tt = (a + b*theta) * cos(theta + pi)，yy_tt = (a + b*theta) * sin(theta + pi)
% 带入参数方程后，令theta = x
if  flag_in == 1 && flag_out == 0
fun = @(x) xx^2 + yy^2 - 2 * xx * (a + b * x) * cos(x) - 2 * yy * (a + b * x) * sin(x) + (a + b * x)^2 - l^2;
elseif flag_out == 1 && flag_in == 0
fun = @(x) xx^2 + yy^2 - 2 * xx * (a + b * x) * cos(x + pi) - 2 * yy * (a + b * x) * sin(x + pi) + (a + b * x)^2 - l^2; 
end

% 扫描区间
span = 5;
% 根据 flag_in / flag_out 选择扫描方向
if flag_in == 1 && flag_out == 0
    x_scan = linspace(theta, theta + span, 200);
elseif flag_out == 1 && flag_in == 0
    x_scan = linspace(theta - span, theta, 200);
else
    error('必须且只能指定一个标志：flag_in=1 或 flag_out=1。');
end

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
if flag_in == 1 && flag_out == 0
    xx_tt = (a + b * x) * cos(x);
    yy_tt = (a + b * x) * sin(x);
elseif flag_out == 1 && flag_in == 0
    xx_tt = (a + b * x) * cos(x + pi);
    yy_tt = (a + b * x) * sin(x + pi);
end

end
