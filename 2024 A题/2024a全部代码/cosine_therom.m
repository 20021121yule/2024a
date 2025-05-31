function y = cosine_therom(x, theta_1, l, a, b)
% 几何关系，对应余弦定理
y = (a + b * x) .^2 + (a + b * theta_1) .^2 - l^2 - 2 * (a + b * x) .* (a + b * theta_1) .* cos(x - theta_1);
end