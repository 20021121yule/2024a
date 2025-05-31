function x = theta_calculate3(theta_1, l, a, b)
% 前后把手都在中心对称曲线的情况
% 更高效的版本，使用 fzero + 初步变号区间扫描

% 目标函数
fun = @(x) cosine_therom(x, theta_1, l, a, b);

% 扫描区间
span = 5;
x_scan = linspace(theta_1 - span, theta_1, 200);
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
