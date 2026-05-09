e = [];
for kk = 0.01:0.01:0.99
    for jj = 0.01:0.01:0.99
        if 1-kk-jj>=0.01
            e = [e;kk,jj,1-kk-jj];
        end
    end
end
% 初始化参数 
R = 1;
S = -0.1;
T = 1.8;
P = 0;
G = 0.8;
W = 0.1;
theta = 1;


beta = 0;

intial = [e(:,1),e(:,2)];

tmax = 4000;
initial_number = size(e,1);
CE_domain = zeros(101,101);

m1 = 0;
for alpha = 0
    m1 = m1+1;
    for delta = 0:0.01:1
        tic
        m2 = round(100*delta)+1;
        for i =1:initial_number
            y = zeros(tmax,2);
            y(1,:) = intial(i,:);
            for t = 1:tmax
                dt = 0.1;
                y1 = y(t,1);
                y2 = y(t,2);
                y3 = 1-y(t,1)-y(t,2);
                pi_c = (1-delta)*(y1*(R/2*(1+theta))+y2*(S)+y3*(R/2*(1+theta)))+...
                    delta*(y1^2*(R/3*(1+theta+theta^2))+2*y1*y2*(G/2*(1+theta))+2*y1*y3*(R/3*(1+theta+theta^2))+y2^2*(S)+2*y2*y3*(G/2*(1+theta))+y3^2*(R/3*(1+theta+theta^2)));
                pi_d = (1-delta)*(y1*(T)+y2*(P)+y3*(T-beta))+...
                    delta*(y1^2*(T/2*(1+theta))+2*y1*y2*(W)+2*y1*y3*(T/2*(1+theta)-beta)+y2^2*(P)+2*y2*y3*(W-beta)+y3^2*(T/2*(1+theta)-2*beta));
                pi_p = (1-delta)*(y1*(R/2*(1+theta))+y2*(S-alpha)+y3*(R/2*(1+theta)))+...
                    delta*(y1^2*(R/3*(1+theta+theta^2))+2*y1*y2*(G/2*(1+theta)-alpha)+2*y1*y3*(R/3*(1+theta+theta^2))+y2^2*(S-2*alpha)+2*y2*y3*(G/2*(1+theta)-alpha)+y3^2*(R/3*(1+theta+theta^2)));
                pi_bar = y1*pi_c + y2*pi_d + y3*pi_p;
                y(t+1,1) = y1 + dt*y1*(pi_c - pi_bar);
                y(t+1,2) = y2 + dt*y2*(pi_d - pi_bar);
                y(t+1,1) = max(0,min(1,y(t+1,1)));
                y(t+1,2) = max(0,min(1,y(t+1,2)));
            end
            
            X = [y 1-y(:,1)-y(:,2)];
            kk = size(X,1);
            CE_domain(m1,m2) = CE_domain(m1,m2) + X(kk,1)+X(kk,3);
        end
        CE_domain(m1,m2) = CE_domain(m1,m2)/initial_number;
        toc
    end
    
end

plot(CE_domain(1,:),'b-', 'LineWidth', 2)
hold on 
plot(CE0109(1,:),'r-', 'LineWidth', 2)
plot(CE0109(2,:),'g-', 'LineWidth', 2)
set(gca, ...
    'FontSize', 16, ...
    'LineWidth', 1.8, ...
    'Box', 'on', ...
    'Layer', 'top');
xlabel('\delta', 'FontSize', 20);
hold off
pbaspect([1 1 1])

