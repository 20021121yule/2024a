function x = theta_calculate1(xx_tt, yy_tt, l, a, b, theta_i)
%theta_calculate1 来计算前把手在AC上，后把手在盘入曲线上的情况
% 目标函数
fun = @(x) xx_tt^2 + yy_tt^2 - 2 * xx_tt * (a + b * x) * cos(x) - 2 * yy_tt * (a + b * x) * sin(x) + (a + b * x)^2 - l^2;

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

end