function [xx_tt, yy_tt] = getNextHandleOn_In_AC(xx, yy, l, a, b, theta_i)
%getNextHandleOn_In_AC 来计算前把手在AC上，后把手在盘入曲线上的情况
% 目标函数
fun = @(x) xx^2 + yy^2 - 2 * xx * (a + b * x) * cos(x) - 2 * yy * (a + b * x) * sin(x) + (a + b * x)^2 - l^2;

% 扫描区间
span = 5;
x_scan = linspace(theta_i, theta_i + span , 200);
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

xx_tt = (a + b * x) * cos(x);
yy_tt = (a + b * x) * sin(x);

end