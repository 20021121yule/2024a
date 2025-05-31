function [xx, yy] = theta_calculate2(xx_tt, yy_tt, O, O_i, radius, l, theta_O_i, theta_Oi_tt)
%theta_calculate2 来计算前把手在AC上，后把手在也在AC上的情况
% O是调转区域圆心坐标
% O_i中i的取值为1、2，分别是O1的圆心坐标或者O2的圆心坐标
% radius是圆的半径
% 直接返回在圆上的空间坐标，因为算出圆上相对于O_i的极角没有任何意义
% theta_O_i 是theta_O1 = linspace(gamma1, gamma1 - theta1, 1000) 根据这里定义的
% theta_Oi_tt 是前把手的极角

fun = @(x) (O_i(1) - O(1) + radius * cos(x) - xx_tt)^2 + (O_i(2) - O(2) + radius * sin(x) - yy_tt)^2 - l^2;

% 扫描区间
span = 60;
x_scan = linspace(theta_Oi_tt - 0.1, max(theta_O_i) + span, 200);
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

% 用计算出来的theta角x来换算空间坐标
xx = O_i(1) - O(1) + radius * cos(x);
yy = O_i(2) - O(2) + radius * sin(x);

end