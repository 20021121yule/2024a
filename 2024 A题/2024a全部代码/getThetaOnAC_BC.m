function x = getThetaOnAC_BC(xx_tt, yy_tt, O, O_i, radius, l, theta_min, theta_max)
%getThetaOnAC_BC 来计算前把手在AC或BC或盘出曲线上时，该把手也在AC或BC上的情况
% O是调转区域圆心坐标
% O_i中i的取值为1、2，分别是O1的圆心坐标或者O2的圆心坐标
% radius是圆的半径
% theta_min是搜索区间最小值
% theta_max是搜索区间最大值

% theta_min:1.当前把手在AC上时，此把手对应相对于弧AC圆的搜索区间最小值应该是从前把手相对O_i的极角开始找，也即theta_min取值为theta_O1_tt(ii - 1)
%           2.当前把手在BC上时，如果此把手在AC上，那么theta_min取值就应该是gamma1 - theta1，因为此时的把手是在弧AC上距离C点最近的把手了（注意按照运动方向，在O1圆上的极角是随时间逐渐减少的）
%           3.当前把手在BC上时，如果此把手在BC上，那么theta_min取值就应该是gamma2，因为对于O2圆来说，按照运动方向，相对于O2的极角是在增大的，所以这个搜索区间最小值应该是gamma2
%           4.当前把手在盘出曲线上时，此把手只有一种可能，即一定在BC上，那么直接可以在BC对应的区间内搜索，因为只有可能有一个解，所以theta_min取值为gamma2。

% theta_max:1.当前把手在AC上时，此把手对应相对于弧AC圆的搜索区间最大值应该是A点相对于O1的极角，也即gamma1
%           2.当前把手在BC上时，如果此把手在AC上，那么此时最大值可以直接取gamma1 - theta1,因为在弧AC上这个区间当中一定只有一个解（几何关系）
%           3.当前把手在BC上时，如果此把手在BC上，那么搜索区间最大值应该是前把手相对于O2的极角，也就是theta_O2_tt(ii - 1)
%           4.当前把手在盘出曲线上时，此把手只有一种可能，即一定在BC上，那么直接可以在BC对应的区间内搜索，因为只有可能有一个解，所以theta_max取值为gamma2 + theta2。

% 对于O1圆，相对于O1的极角取值范围是[gamma1, gamma1 - theta1](按时间推进，由大到小)
% 对于O2圆，相对于O2的极角取值范围是[gamma2, gamma2 + theta2](按时间推进，由小到大)


fun = @(x) (O_i(1) - O(1) + radius * cos(x) - xx_tt)^2 + (O_i(2) - O(2) + radius * sin(x) - yy_tt)^2 - l^2;

% 扫描区间
eps = 0.1; % 注意这里的eps保证theta_min和theta_max也是这个方程的解
x_scan = linspace(theta_min - eps, theta_max + eps, 200);
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