R = 1; S = -0.1; T = 1.8; P = 0; G = 0.8; W = 0.1; theta = 1;
alpha = 0.7;
beta = 0.7;
a = 2*(G-W);
b = T-S-1;
c = a+b;
delta_1 = (b+S)/beta-1;
delta_2 = (b+sqrt(-4*S*(b+S)))/c;
x_c1 = zeros(1,100);
x_c2 = zeros(1,100);
x_c3 = zeros(1,100);
x_p = zeros(1,100);
k = 1;
for delta = 0.00:0.0001:1
    x_c1(k) = 1-(b+S)/((1+delta)*beta);
    x_c2(k) = ((c*delta-b-2*S)+sqrt((c*delta-b)^2+4*S*(b+S)))/(2*c*delta);
    x_c3(k) = ((c*delta-b-2*S)-sqrt((c*delta-b)^2+4*S*(b+S)))/(2*c*delta);
    x_p(k) = (c*delta-b-2*S+(1+delta)*(alpha+beta)-sqrt((c*delta-b-2*S+(1+delta)*(alpha+beta))^2+4*c*delta*(S-(1+delta)*alpha)))/(2*c*delta);
    k = k+1;
end


delta_X1 = 0.1430:0.0001:1;
X1 = 1431:10001;
delta_X2 = 0.6372:0.0001:1;
X2 = 6373:10001;
delta_X3 = 0.6372:0.0001:1;
X3 = 6373:10001;
y1 = x_c1(X1);
plot(delta_X1,x_c1(X1),'k-','LineWidth', 2)
hold on 
fill([delta_X1, fliplr(delta_X1)], [y1, zeros(size(y1))], ...
     [0.9 0.9 0.9], 'EdgeColor', 'none')   
hold off
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])
box on


delta_X1 = 0.1430:0.0001:1;
X1 = 1431:10001;
delta_X2 = 0.6372:0.0001:1;
X2 = 6373:10001;
delta_X3 = 0.6372:0.0001:1;
X3 = 6373:10001;
plot(delta_X2,x_c2(X2),'k-','LineWidth', 2)
hold on 
plot(delta_X3,x_c3(X3),'k--','LineWidth', 2)
plot([0 1],[0 0],'k-','LineWidth', 2)
plot([0.6372 0.6372],[0 0.2649],'k:','LineWidth', 1.2)
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])
hold off

delta_Xp = 0.1430:0.0001:1;
Xp = 1431:10001;
plot(delta_Xp,x_p(Xp),'k--','LineWidth', 2)
hold on
plot([0 1],[1 1],'k-','LineWidth', 2)
plot([0.1430 1],[0 0],'k-','LineWidth', 2)
plot([0.1430 0.1430],[0 1],'k:','LineWidth', 1.2)
hold off
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])
