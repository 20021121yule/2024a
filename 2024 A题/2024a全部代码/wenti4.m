clear;

radius = 1.7; % 螺距
a = 0;
b = 1 / (2 * pi) * radius;

theta_a = linspace(10*2*pi, 0,1000);
r_a = a + b * theta_a;
% 曲线(正)
xx1 = r_a .* cos(theta_a);
yy1 = r_a .* sin(theta_a);

% 曲线(中心对称)
xx2 = r_a .* cos(theta_a + pi);
yy2 = r_a .* sin(theta_a + pi);

% 调转区域
r = 4.5; % 调转区域半径
theta_r = linspace(0, 2 * pi, 1000);
% 调转区域二维分布
xx3 = r .* cos(theta_r);
yy3 = r .* sin(theta_r);

% 计算交点，交点规定为A点和B点
theta_intersection = r / b; % 交点的theta
x_A = (a + b * theta_intersection) * cos(theta_intersection);
y_A = (a + b * theta_intersection) * sin(theta_intersection);
A = [x_A, y_A]; % A点坐标

x_B = -x_A;
y_B = -y_A;
B = [x_B, y_B]; % B点坐标

% 计算路径，先计算张角theta1
k1 = ((a + b * theta_intersection) * sin(theta_intersection) - b * cos(theta_intersection)) / ((a + b * theta_intersection) * cos(theta_intersection) + b * sin(theta_intersection)); % 这是盘入曲线的法线的斜率
k2 = (y_B - y_A) / (x_B - x_A); % 这是两个交点连线的斜率

% 应用k1和k2一起求出夹角alpha
alpha = atan(abs((k2 - k1) / (1 + k1 * k2)));
theta1 = pi - 2 * alpha; % 张角theta1

x = r / 3; % OC的长度
R = (r + x) / (4 * cos(alpha)); % 小圆半径
R_large = 2 * R; % 大圆半径

% 求O1坐标
% 先求A点切向量
t_A = [((a + b * theta_intersection) * sin(theta_intersection) - b * cos(theta_intersection)), -((a + b * theta_intersection) * cos(theta_intersection) + b * sin(theta_intersection))];
t_A = t_A / norm(t_A); % A点切向量
n_A = [t_A(2), -t_A(1)]; % A点法向量
O1 = A + n_A * R_large; % O1圆心坐标

% 求C和O的坐标
n_AB = [x_B - x_A, y_B - y_A];
n_AB = n_AB / norm(n_AB);
len_AC = r + x;
C = A + n_AB * (r + x); % C的坐标
O = A + n_AB * r; % 圆心O的坐标

% O1坐标系的极角和O坐标系的极角的变换关系
gamma1 = atan2(y_A - O1(2), x_A - O1(1));
theta_O1 = linspace(gamma1, gamma1 - theta1, 1000); % thetaO1的取值范围

O1_O = O1 - O; % 圆心O到圆心O1的向量
G_O1_xx = R_large * cos(theta_O1); % 圆心1到动点G的向量(x分量)
G_O1_yy = R_large * sin(theta_O1); % 圆心1到动点G的向量(y分量)

% 在O1段圆弧的坐标
xx_O1 = O1_O(1) + G_O1_xx;
yy_O1 = O1_O(2) + G_O1_yy;

% 求O2坐标
n_B = -n_A; %中心对称关系
O2 = B + n_B * R; % O2圆心坐标

gamma2 = atan2(C(2) - O2(2), C(1) - O2(1));
theta2 = theta1; % 几何关系，几何图形相似性
theta_O2 = linspace(gamma2, gamma2 + theta2, 1000); % thetaO2的取值范围

O2_O = O2 - O; % 圆心O到圆心O1的向量
I_O2_xx = R * cos(theta_O2); % 圆心1到动点G的向量(x分量)
I_O2_yy = R * sin(theta_O2); % 圆心1到动点G的向量(y分量)

% 在O2段圆弧的坐标
xx_O2 = O2_O(1) + I_O2_xx;
yy_O2 = O2_O(2) + I_O2_yy;

% 龙把手速度求解
% 龙头随时间的运动
delta_t = 0.1; % 注意这里，精度问题，我们采用delta_t = 1,如果需要更好的精度就再调小即可
start_time = 0; % 启动时间
end_time = 500; % 终止时间
step = (end_time - start_time) / delta_t; % 时间步

% 暂时用来计算：由龙头来确定其余部分
head_num = 223;
l1 = 2.86; % 第一节长度，把手间距
l2 = 1.65; % 第二节长度，把手间距
l = 0; % 占位，遍历时长度

theta_inistialize = 3*2*pi; % 在t = 0时的初始角
theta_in_tt = zeros(1, head_num + 1); % 这个theta储存223节把手在盘入曲线上的theta值
theta_O1_tt = zeros(1, head_num + 1); % 这个theta储存223节把手在O1圆上的theta值
theta_O2_tt = zeros(1, head_num + 1); % 这个theta储存223节把手在O2圆上的theta值
theta_out_tt = zeros(1, head_num + 1); % 这个theta储存223节把手在盘出曲线上的theta值

time1 = NaN; % 从盘入曲线进入大圆的时刻
time2 = NaN; % 从大圆盘入小圆的时刻
time3 = NaN; % 从小圆盘出进入中心对称曲线的时刻

flag_in = ones(1, head_num + 1); % 判定所有的把手是否进入in段，初始时所有都为in
flag_AC = zeros(1, head_num + 1); % 判定所有的把手是否进入AC段
flag_BC = zeros(1, head_num + 1); % 判定所有的把手是否进入BC段
flag_out = zeros(1, head_num + 1); % 判定所有的把手是否进入BC段

xx_tt = zeros(1, head_num + 1); % head_num + 1个把手的横坐标
yy_tt = zeros(1, head_num + 1); % head_num + 1个把手的纵坐标

eps = 0.1; % 曲线判定容限

% --- Dynamic Plotting Setup ---
figure;
hold on;

plot(xx1, yy1, 'k-');
plot(xx2, yy2, 'r-');
plot(xx3, yy3, '-','LineWidth', 2.0)

% O1和O2段路径
plot(xx_O1, yy_O1, 'b-');
plot(xx_O2, yy_O2, 'g-');

% 标记 A B C O1 O2点
plot(x_A, y_A, 'k.', 'MarkerSize', 10, 'LineWidth', 2);  % 加粗大黑点
text(x_A + 0.1, y_A + 0.1, 'A', 'FontSize', 10);          % 正常字体，小号标注
plot(x_B, y_B, 'r.', 'MarkerSize', 10, 'LineWidth', 2);  % 加粗大红点
text(x_B + 0.1, y_B + 0.1, 'B', 'FontSize', 10);
plot(C(1), C(2), 'b.', 'MarkerSize', 10, 'LineWidth', 2);  % 加粗大蓝点
text(C(1) + 0.1, C(2) + 0.1, 'C', 'FontSize', 10);

% Initialize plot handles for the dynamic elements
h_line = plot(NaN, NaN, 'k-', 'LineWidth', 2.0, 'DisplayName', '板凳连接线');
h_points = plot(NaN, NaN, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', '板凳位置');
h_title = title('动态板凳位置'); % Dynamic title to show current time

axis equal;
xlabel('x'); ylabel('y');
grid on;
hold off;

% Get the axes handle to control the view
ax = gca;
% Set initial axis limits if you want a fixed view
% xlim([min(xx) max(xx)]);
% ylim([min(yy) max(yy)]);

% --- Animation Loop ---
for tt = 1:step
    time = (tt - 1) * delta_t;

    % === 段1：阿基米德盘入 ===
    if isnan(time1)
        theta_inistialize = theta_inistialize - (1/sqrt(b^2 + (a + b * theta_inistialize)^2)) * delta_t;
        theta_in_tt(1) = theta_inistialize; % 每次迭代时初始值是tt时刻下的龙头的theta值

        % 龙头前把手坐标
        xx_tt(1) = (a + b*theta_in_tt(1)) .* cos(theta_in_tt(1));
        yy_tt(1) = (a + b*theta_in_tt(1)) .* sin(theta_in_tt(1));

        for ii = 2:head_num + 1
            % 先检查是否是第一节，如果是，那么用第一节长度
            if ii == 2
                l = l1;
            else
                l = l2; % 否则，用第二节长度
            end

            [xx_tt(ii), yy_tt(ii)] = getNextHandleOnSpiral(xx_tt(ii - 1), yy_tt(ii - 1), a, b, l, flag_in(ii), flag_out(ii));
        end

        % 判断是否进入大圆弧（到达A点）
        if theta_inistialize <= theta_intersection
            time1 = time;
        end

        % === 段2：大圆弧 A→C ===
    elseif isnan(time2)
        % 龙头把手在O1上相对于O1的角度
        theta_O1_tt(1) = 1 / R_large * (time1 - time) + gamma1;

        % 龙头把手坐标
        xx_tt(1) = O1_O(1) + R_large * cos(theta_O1_tt(1));
        yy_tt(1) = O1_O(2) + R_large * sin(theta_O1_tt(1));

        % 龙头把手的flag
        flag_in(1) = 0;
        flag_AC(1) = 1;

        % 判断是否进入小圆弧（到达C点）
        if theta_O1_tt(1) <= gamma1 - theta1
            time2 = time;
        end

        for ii = 2:head_num + 1
            % 先检查是否是第一节，如果是，那么用第一节长度
            if ii == 2
                l = l1;
            else
                l = l2; % 否则，用第二节长度
            end

            if flag_AC(ii) == 1
                % 如果第ii个把手也在AC上
                theta_O1_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O1, R_large, l, theta_O1_tt(ii - 1) , gamma1);

                xx_tt(ii) = O1(1) - O(1) + R_large * cos(theta_O1_tt(ii));
                yy_tt(ii) = O1(2) - O(2) + R_large * sin(theta_O1_tt(ii));
                continue; % 计算完直接下一个循环
            end

            if flag_in(ii) == 1
                % 如果第ii个把手在盘入曲线，分两种情况，1.前一个把手在AC上、2.前一个把手也在盘入曲线
                if flag_AC(ii - 1) == 1
                    % 前一个把手在AC上
                    [xx_tt(ii), yy_tt(ii)] = getNextHandleOn_In_AC(xx_tt(ii - 1), yy_tt(ii - 1), l, a, b, theta_intersection);

                    dis1 = sqrt((xx_tt(ii) - A(1))^2 + (yy_tt(ii) - A(2))^2); % dis1代表把手与A点距离
                    if dis1 < eps
                        % 满足该条件，那么说明进入盘出曲线
                        flag_AC(ii) = 1;
                        flag_in(ii) = 0;
                    end

                    continue; % 计算完直接下一个循环
                end

                if flag_in(ii - 1) == 1
                    % 前一个把手也在盘入曲线
                    [xx_tt(ii), yy_tt(ii)] = getNextHandleOnSpiral(xx_tt(ii - 1), yy_tt(ii - 1), a, b, l, flag_in(ii), flag_out(ii));
                    continue;
                end

            end

        end

        % === 段3：小圆弧 C→B ===
    elseif isnan(time3)
        theta_O2_tt(1) = 1 / R * (time - time2) + gamma2;

        % 龙头把手坐标
        xx_tt = O2_O(1) + R * cos(theta_O2_tt(1));
        yy_tt = O2_O(2) + R * sin(theta_O2_tt(1));

        flag_AC(1) = 0;
        flag_BC(1) = 1;

        % 判断是否进入中心对称螺线（到达B点）
        if theta_O2_tt(1) >= gamma2 + theta2
            time3 = time;
            theta_out_tt(1) = theta_intersection;  % 初始化中心对称段起始角度
        end


        for ii = 2:head_num + 1
            % 先检查是否是第一节，如果是，那么用第一节长度
            if ii == 2
                l = l1;
            else
                l = l2; % 否则，用第二节长度
            end

            if flag_BC(ii) == 1
                % 如果第ii个把手也在BC上，只有一种情况，前一个把手也一定在BC上
                theta_O2_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O2, R, l, gamma2 , theta_O2_tt(ii - 1));

                xx_tt(ii) = O2(1) - O(1) + R * cos(theta_O2_tt(ii));
                yy_tt(ii) = O2(2) - O(2) + R * sin(theta_O2_tt(ii));

                continue;
            end

            if flag_AC(ii) == 1
                % 如果第ii个把手也在AC上
                if flag_AC(ii - 1) == 1
                    % 如果前一一个把手也在AC上
                    theta_O1_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O1, R_large, l, theta_O1_tt(ii - 1) , gamma1);
                end

                if flag_BC(ii - 1) == 1
                    % 如果前一个把手在BC上
                    theta_O1_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O1, R_large, l, gamma1 - theta1 , gamma1);
                end

                xx_tt(ii) = O1(1) - O(1) + R_large * cos(theta_O1_tt(ii));
                yy_tt(ii) = O1(2) - O(2) + R_large * sin(theta_O1_tt(ii));

                % 这个时候就要判定是否进入BC段了
                dis2 = sqrt((xx_tt(ii) - C(1))^2 + (yy_tt(ii) - C(2))^2); % dis2代表把手与C点的距离
                if dis2 < eps
                    % 进入BC段
                    flag_AC(ii) = 0;
                    flag_BC(ii) = 1;
                end

                continue; % 计算完直接下一个循环
            end

            if flag_in(ii) == 1
                % 如果第ii个把手在盘入曲线，分两种情况，1.前一个把手在AC上、2.前一个把手也在盘入曲线
                if flag_AC(ii - 1) == 1
                    % 前一个把手在AC上
                    [xx_tt(ii), yy_tt(ii)] = getNextHandleOn_In_AC(xx_tt(ii - 1), yy_tt(ii - 1), l, a, b, theta_intersection);

                    dis1 = sqrt((xx_tt(ii) - A(1))^2 + (yy_tt(ii) - A(2))^2); % dis1代表把手与A点距离
                    if dis1 < eps
                        % 满足该条件，那么说明进入盘出曲线
                        flag_AC(ii) = 1;
                        flag_in(ii) = 0;
                    end

                    continue; % 计算完直接下一个循环
                end

                if flag_in(ii - 1) == 1
                    % 前一个把手也在盘入曲线
                    [xx_tt(ii), yy_tt(ii)] = getNextHandleOnSpiral(xx_tt(ii - 1), yy_tt(ii - 1), a, b, l, flag_in(ii), flag_out(ii));
                    continue;
                end
            end
        end


        % === 段4：中心对称阿基米德盘出 ===
        % 分两段，1.前把手在盘出曲线上和后把手在小圆弧 C→B上，2.前后把手都在前把手在盘出曲线上
    else
        theta_out_tt(1) = theta_out_tt(1) + 1 / sqrt(b^2 + (a + b * theta_out_tt(1))^2) * delta_t;

        % 龙头把手坐标
        xx_tt = (a + b * theta_out_tt(1)) * cos(theta_out_tt(1) + pi);
        yy_tt = (a + b * theta_out_tt(1)) * sin(theta_out_tt(1) + pi);

        flag_BC(1) = 0;
        flag_out(1) = 1;

        for ii = 2:head_num + 1
            % 先检查是否是第一节，如果是，那么用第一节长度
            if ii == 2
                l = l1;
            else
                l = l2; % 否则，用第二节长度
            end

            if flag_out(ii) == 1
                % 第ii把手在盘出曲线上
                [xx_tt(ii),yy_tt(ii)] = getNextHandleOnSpiral(xx_tt(ii - 1), yy_tt(ii - 1), a, b, l, flag_in(ii), flag_out(ii));
            end

            if flag_BC(ii) == 1
                % 有2种情况，前把手在BC上和前把手在盘出曲线上
                if flag_BC(ii - 1) == 1
                    % 前把手在BC上
                    theta_O2_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O2, R, l, gamma2 , theta_O2_tt(ii - 1));
                end

                if flag_out(ii - 1) == 1
                    % 前把手在盘出曲线上
                    theta_O2_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O2, R, l, gamma2 , gamma2 + theta2);
                end

                xx_tt(ii) = O2(1) - O(1) + R * cos(theta_O2_tt(ii));
                yy_tt(ii) = O2(2) - O(2) + R * sin(theta_O2_tt(ii));

                dis3 = sqrt((xx_tt(ii) - B(1))^2 + (yy_tt(ii) - B(2))^2); % dis3代表把手与B点距离
                if dis3 < eps
                    % 满足该条件，那么说明进入盘出曲线
                    flag_BC(ii) = 0;
                    flag_out(ii) = 1;
                end

                continue;
            end

            if flag_AC(ii) == 1
                % 如果第ii个把手也在AC上
                if flag_AC(ii - 1) == 1
                    % 如果前一一个把手也在AC上
                    theta_O1_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O1, R_large, l, theta_O1_tt(ii - 1) , gamma1);
                end

                if flag_BC(ii - 1) == 1
                    % 如果前一个把手在BC上
                    theta_O1_tt(ii) = getThetaOnAC_BC(xx_tt(ii - 1), yy_tt(ii - 1), O, O1, R_large, l, gamma1 - theta1 , gamma1);
                end

                xx_tt(ii) = O1(1) - O(1) + R_large * cos(theta_O1_tt(ii));
                yy_tt(ii) = O1(2) - O(2) + R_large * sin(theta_O1_tt(ii));

                % 这个时候就要判定是否进入BC段了
                dis2 = sqrt((xx_tt(ii) - C(1))^2 + (yy_tt(ii) - C(2))^2); % dis2代表把手与C点的距离
                if dis2 < eps
                    % 进入BC段
                    flag_AC(ii) = 0;
                    flag_BC(ii) = 1;
                end

                continue; % 计算完直接下一个循环
            end

            if flag_in(ii) == 1
                % 如果第ii个把手在盘入曲线，分两种情况，1.前一个把手在AC上、2.前一个把手也在盘入曲线
                if flag_AC(ii - 1) == 1
                    % 前一个把手在AC上
                    [xx_tt(ii), yy_tt(ii)] = getNextHandleOn_In_AC(xx_tt(ii - 1), yy_tt(ii - 1), l, a, b, theta_intersection);

                    dis1 = sqrt((xx_tt(ii) - A(1))^2 + (yy_tt(ii) - A(2))^2); % dis1代表把手与A点距离
                    if dis1 < eps
                        % 满足该条件，那么说明进入盘出曲线
                        flag_AC(ii) = 1;
                        flag_in(ii) = 0;
                    end

                    continue; % 计算完直接下一个循环
                end

                if flag_in(ii - 1) == 1
                    % 前一个把手也在盘入曲线
                    [xx_tt(ii), yy_tt(ii)] = getNextHandleOnSpiral(xx_tt(ii - 1), yy_tt(ii - 1), a, b, l, flag_in(ii), flag_out(ii));
                    continue;
                end
            end
        end

    end

    % === 更新绘图 ===
    set(h_line, 'XData', xx_tt, 'YData', yy_tt);
    set(h_points, 'XData', xx_tt, 'YData', yy_tt);
    set(h_title, 'String', sprintf('当前时间：%.2f s', time));
    drawnow limitrate;
end

fprintf('动画播放完毕！\n');



