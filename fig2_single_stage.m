R = 1; S = -0.1; T = 1.8; P = 0; G = 0.8; W = 0.1; theta = 1.0;
alpha = 0.49;
beta = 0.44;
a = 2*(G-W);
b = T-S-1;
c = a+b;
delta_1 = (b+S)/beta-1;
delta_2 = (b+sqrt(-4*S*(b+S)))/c;
x_c1 = zeros(1,100);
x_c2 = zeros(1,100);
x_c3 = zeros(1,100);
x_p = zeros(1,100);
x_pz = zeros(1,100);
k = 1;
for delta = 0.00:0.0001:1
    x_c1(k) = 1-(b+S)/((1+delta)*beta);
    x_c2(k) = ((c*delta-b-2*S)+sqrt((c*delta-b)^2+4*S*(b+S)))/(2*c*delta);
    x_c3(k) = ((c*delta-b-2*S)-sqrt((c*delta-b)^2+4*S*(b+S)))/(2*c*delta);
    x_p(k) = (c*delta-b-2*S+(1+delta)*(alpha+beta)-sqrt((c*delta-b-2*S+(1+delta)*(alpha+beta))^2+4*c*delta*(S-(1+delta)*alpha)))/(2*c*delta);
    x_pz(k) = (c*delta-b-2*S+(1+delta)*(alpha+beta)+sqrt((c*delta-b-2*S+(1+delta)*(alpha+beta))^2+4*c*delta*(S-(1+delta)*alpha)))/(2*c*delta);
    k = k+1;
end


delta_X1 = 0.8182:0.0001:1;
X1 = 8183:10001;
delta_X2 = 0.6372:0.0001:1;
X2 = 6373:10001;
delta_X3 = 0.6372:0.0001:1;
X3 = 6373:10001;
y1 = x_c1(X1);
hold on
plot(delta_X1,x_c1(X1),'k-','LineWidth', 2)
fill([delta_X1, fliplr(delta_X1)], [y1, zeros(size(y1))], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none')
set(gca,'Layer','top')
hold off
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])
box on


delta_X1 = 0.8182:0.0001:1;
X1 = 8183:10001;
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

delta_Xp = 0.6518:0.0001:1;
Xp = 6519:10001;
delta_Xpz = 0.6518:0.0001:0.8182;
Xpz = 6519:8183;
plot(delta_Xp,x_p(Xp),'k-','LineWidth', 2)
hold on
plot(delta_Xpz,x_pz(Xpz),'k--','LineWidth', 2)
plot([0 0.8182],[1 1],'k-','LineWidth', 2)
plot([0.6518 0.6518],[0 0.7748],'k:','LineWidth', 1)
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])