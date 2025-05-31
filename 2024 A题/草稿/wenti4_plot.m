clear;

radius = 1.7; % 螺距
a = 0;
b = 1 / (2 * pi) * radius;

theta_a = linspace(16*2*pi, 0,1000);
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
theta_i = r / b; % 交点的theta
x_A = (a + b * theta_i) * cos(theta_i);
y_A = (a + b * theta_i) * sin(theta_i);
A = [x_A, y_A]; % A点坐标

x_B = -x_A;
y_B = -y_A;
B = [x_B, y_B]; % B点坐标

% 计算路径，先计算张角theta1
k1 = ((a + b * theta_i) * sin(theta_i) - b * cos(theta_i)) / ((a + b * theta_i) * cos(theta_i) + b * sin(theta_i)); % 这是盘入曲线的法线的斜率
k2 = (y_B - y_A) / (x_B - x_A); % 这是两个交点连线的斜率

% 应用k1和k2一起求出夹角alpha
alpha = atan(abs((k2 - k1) / (1 + k1 * k2)));
theta1 = pi - 2 * alpha; % 张角theta1

x = r / 3; % OC的长度
R = (r + x) / (4 * cos(alpha)); % 小圆半径
R_large = 2 * R; % 大圆半径

% 求O1坐标
% 先求A点切向量
t_A = [((a + b * theta_i) * sin(theta_i) - b * cos(theta_i)), -((a + b * theta_i) * cos(theta_i) + b * sin(theta_i))];
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

axis([-5 5 -5 5]);
axis equal
hold off;

