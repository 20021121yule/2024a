clear;
% 等距螺线参数，r = a + b * theta
radius = 0.55;
a = 0;
b = 1 / (2 * pi) * radius;

% 龙头随时间的运动
delta_t = 0.1; % 注意这里，精度问题，我们采用delta_t = 1,如果需要更好的精度就再调小即可
start_time = 0; % 启动时间
end_time = 400; % 终止时间
step = (end_time - start_time) / delta_t; % 时间步

% 注意在我们的规定当中，theta是逐渐在减小的（是负数），取值范围是[0, - 16 * 2pi]。
theta_t = 16 * 2 * pi * ones(1,step + 1); % 0-300s所有的theta取值都在这里了，我们再用坐标变换再求其余坐标即可

% 暂时用来计算：由龙头来确定其余部分
head_num = 223; % head_num为前把手的值，最后一个把手也算一个前把手
l1 = 2.26; % 第一节长度,l1是把手长度，实际长度才是3.40
l11 = 3.41;
l2 = 1.65; % 第二节长度，l2是把手长度，实际长度才是2.20
l22 = 2.20;
w = 0.3; % 宽度
l = l2 * ones(1,2); % 注意这里的l在第二问中需要变成向量
l33 = l22 * ones(1,2);

theta = zeros(1, head_num + 1); % 这个theta储存224节的theta值

flag = 1;

% A B D E 储存每个板凳的4个端点，大小为1 * headnum
A = zeros(2,head_num);
B = zeros(2,head_num);
D = zeros(2,head_num);
E = zeros(2,head_num);

tt = 1;
time = 0.00;

X = zeros(5, head_num);
Y = zeros(5, head_num);

theta1 = linspace(0, - 16*2*pi,1000);
r = a + b * theta1;
xx = r .* cos(theta1);
yy = r .* sin(theta1);
% --- Dynamic Plotting Setup ---
figure;
hold on;
% 螺线轨迹：细灰色虚线 (plot once)
plot(xx, yy, '-', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.0, 'DisplayName', '螺线轨迹');

% Initialize plot handles for the dynamic elements
h_line = plot(NaN, NaN, 'k-', 'LineWidth', 2.0, 'DisplayName', '板凳连接线');
h_border = plot(NaN, NaN, 'b-', 'LineWidth', 0.6, 'DisplayName', '板凳边界');
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

while(flag)
    theta_t(tt+1) = - (1/sqrt(b^2 + (a+b*theta_t(tt))^2)) * delta_t + theta_t(tt);% 随时间迭代, 这是龙头的theta值
    theta(1) = theta_t(tt+1); % 每次迭代时初始值是tt时刻下的龙头的theta值

    for ii = 1:head_num
        % 先检查是否是第一节，如果是，那么用第一节长度
        if ii == 1
            l(1) = l1;
            l33(1) = l11;
        else
            l(1) = l2; % 否则，用第二节长度
            l33(1) = l2;
        end

        theta(ii + 1) = theta_calculate(theta(ii), l(1), a, b); % 用l(1)来计算由ii递推到ii+1角的值

        theta_current1 = zeros(1,2); % 临时储存ii和ii+1
        theta_current1(1) = theta(ii);
        theta_current1(2) = theta(ii + 1);

        [A(:,ii), B(:,ii), D(:,ii), E(:,ii)] = terminal_calculate(theta_current1, l33(1), w, a, b); % 用l(1)来确定这个矩形
    end

    A_current = zeros(2,1);
    B_current = zeros(2,1);
    D_current = zeros(2,1);
    E_current = zeros(2,1);

    % 检测数目
    detected_num = 20; % 注意是，从头的下一个开始的3个，也就是说从第三个板凳开始数3个。

    theta_current4 = zeros(1,2); %theta_current4表示要检测的板凳的2个极角
    % 第1个板凳的角点
    X1 = [A(1, 1), B(1, 1), D(1, 1), E(1, 1), A(1, 1)];
    Y1 = [A(2, 1), B(2, 1), D(2, 1), E(2, 1), A(2, 1)];

    for jj = 1:detected_num
        theta_current4(1) = theta(2 + jj);
        theta_current4(2) = theta(2 + jj + 1);

        [A_current(:,1), B_current(:,1), D_current(:,1), E_current(:,1)] = terminal_calculate(theta_current4, l33(2), w, a, b);
        X_current = [A_current(1, 1), B_current(1, 1), D_current(1, 1), E_current(1, 1), A_current(1, 1)];
        Y_current = [A_current(2, 1), B_current(2, 1), D_current(2, 1), E_current(2, 1), A_current(2, 1)];

        % 情况一：第1个板凳的角点是否落入当前板凳
        s1 = inpolygon(A(1, 1), A(2, 1), X_current, Y_current);
        s2 = inpolygon(B(1, 1), B(2, 1), X_current, Y_current);
        s3 = inpolygon(D(1, 1), D(2, 1), X_current, Y_current);
        s4 = inpolygon(E(1, 1), E(2, 1), X_current, Y_current);
        s = [s1; s2; s3; s4];

        % 情况二：当前板凳的角点是否落入第1个板凳
        t1 = inpolygon(A_current(1,1), A_current(2,1), X1, Y1);
        t2 = inpolygon(B_current(1,1), B_current(2,1), X1, Y1);
        t3 = inpolygon(D_current(1,1), D_current(2,1), X1, Y1);
        t4 = inpolygon(E_current(1,1), E_current(2,1), X1, Y1);
        t = [t1; t2; t3; t4];

        % 任意方向有点落入都算碰撞
        if any(s) || any(t)
            flag = 0;
            fprintf('发生碰撞：第1个板凳 与 第 %d 个板凳相撞\n', jj + 2);
            break;
        end
    end

    tt = tt + 1;
    time = time + delta_t;

    % Update the x and y coordinates for the current time step
    xx_tt = (a + b*theta) .* cos(theta);
    yy_tt = (a + b*theta) .* sin(theta);

    % Update the plot data
    set(h_line, 'XData', xx_tt, 'YData', yy_tt);
    set(h_points, 'XData', xx_tt, 'YData', yy_tt);
    set(h_title, 'String', sprintf('当前时间：%.2f s', time));

    % 连线画出每个板凳矩形轮廓
    for ii = 1:head_num
        % 每个点都是 [x; y]，拼成闭合多边形
        X(:,ii) = [A(1, ii), B(1, ii), D(1, ii), E(1, ii), A(1, ii)]';
        Y(:,ii) = [A(2, ii), B(2, ii), D(2, ii), E(2, ii), A(2, ii)]';
    end

    % 拼接所有矩形为一条长向量，并用 NaN 分隔
    X_border = [];
    Y_border = [];
    for ii = 1:head_num
        X_border = [X_border, X(:,ii)', NaN];
        Y_border = [Y_border, Y(:,ii)', NaN];
    end

    set(h_border, 'XData', X_border, 'YData', Y_border); % 正确：传的是向量

    drawnow limitrate; % Force MATLAB to update the plot
    % pause(0.01); % Optional: Add a small pause to control animation speed

end

figure;
hold on;
% 画所有端点
plot(A(1, :), A(2, :), 'r.', 'DisplayName', 'A 点');
plot(B(1, :), B(2, :), 'g.', 'DisplayName', 'B 点');
plot(D(1, :), D(2, :), 'b.', 'DisplayName', 'D 点');
plot(E(1, :), E(2, :), 'm.', 'DisplayName', 'E 点');

% 连线画出每个板凳矩形轮廓
for ii = 1:head_num
    % 每个点都是 [x; y]，拼成闭合多边形
    X = [A(1, ii), B(1, ii), D(1, ii), E(1, ii), A(1, ii)];
    Y = [A(2, ii), B(2, ii), D(2, ii), E(2, ii), A(2, ii)];

    plot(X, Y, 'k-'); % 黑色轮廓线
end

title('相撞后所有板凳的端点和轮廓');
xlabel('x');
ylabel('y');
axis equal;

hold off;



