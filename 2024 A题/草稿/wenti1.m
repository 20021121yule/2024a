clear;
% 等距螺线参数，r = a + b * theta
radius = 0.55;
a = 0;
b = 1 / (2 * pi) * radius;

theta1 = linspace(16*2*pi, 0 ,1000);
r = a + b * theta1;
xx = r .* cos(theta1);
yy = r .* sin(theta1);

% 龙头随时间的运动
delta_t = 1; % 注意这里，精度问题，我们采用delta_t = 1,如果需要更好的精度就再调小即可
start_time = 0; % 启动时间
end_time = 300; % 终止时间
step = (end_time - start_time) / delta_t; % 时间步

% 注意在我们的规定当中，theta是逐渐在减小的（是负数），取值范围是[0, - 16 * 2pi]。
theta_t = 16 * 2 * pi * ones(1,step + 1); % 0-300s所有的theta取值都在这里了，我们再用坐标变换再求其余坐标即可

% 暂时用来计算：由龙头来确定其余部分
head_num = 223;
l1 = 3.41; % 第一节长度
l2 = 2.20; % 第二节长度
l = 0; % 占位，遍历时长度
theta = zeros(1, head_num + 1); % 这个theta储存223节的theta值

% --- Dynamic Plotting Setup ---
figure;
hold on;
% 螺线轨迹：细灰色虚线 (plot once)
plot(xx, yy, '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.0, 'DisplayName', '螺线轨迹');

% Initialize plot handles for the dynamic elements
h_line = plot(NaN, NaN, 'k-', 'LineWidth', 2.0, 'DisplayName', '板凳连接线');
h_points = plot(NaN, NaN, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', '板凳位置');
h_title = title('动态板凳位置'); % Dynamic title to show current time

axis equal;
xlabel('x'); ylabel('y');
legend('Location', 'best');
grid on;
hold off;

% Get the axes handle to control the view
ax = gca;
% Set initial axis limits if you want a fixed view
% xlim([min(xx) max(xx)]);
% ylim([min(yy) max(yy)]);

% --- Animation Loop ---
for tt = 1:step
    time = (tt - 1) * delta_t; % Calculate current time
    theta_t(tt+1) = - (1/sqrt(b^2 + (a+b*theta_t(tt))^2)) * delta_t + theta_t(tt); % 随时间迭代, 这是龙头的theta值
    theta(1) = theta_t(tt+1); % 每次迭代时初始值是tt时刻下的龙头的theta值

    for ii = 1:head_num
        % 先检查是否是第一节，如果是，那么用第一节长度
        if ii == 1
            l = l1;
        else
            l = l2; % 否则，用第二节长度
        end
        theta(ii + 1) = theta_calculate(theta(ii), l, a, b);
    end

    % Update the x and y coordinates for the current time step
    xx_tt = (a + b*theta) .* cos(theta);
    yy_tt = (a + b*theta) .* sin(theta);

    % Update the plot data
    set(h_line, 'XData', xx_tt, 'YData', yy_tt);
    set(h_points, 'XData', xx_tt, 'YData', yy_tt);
    set(h_title, 'String', sprintf('当前时间：%.2f s', time));

    drawnow limitrate; % Force MATLAB to update the plot
    % pause(0.01); % Optional: Add a small pause to control animation speed
end

fprintf('动画播放完毕！\n');