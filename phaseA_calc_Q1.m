clear all;
close all;
clc;

%% Q1: obtaining u

% ws = warm start
% use this uk as the start for the 2nd time: warm start (a better initial guess)

T = 5;
h = 0.01;
N = T/h;
uk_ini = zeros(N+1,1);

[uk_ws,xk_ws,t_ws] = NR_Dual_Pendulum_Swingup_revised(T,uk_ini);

figure(2)
sz = 12;
plot(t_ws,xk_ws(1,:),"LineWidth",1)
hold on
plot(t_ws,xk_ws(2,:),"LineWidth",1)
plot(t_ws,xk_ws(3,:),"LineWidth",1)
plot(t_ws,xk_ws(4,:),"LineWidth",1)
plot(t_ws,xk_ws(5,:),"LineWidth",1)
plot(t_ws,xk_ws(6,:),"LineWidth",1)
set(gca, 'FontName', 'times')
set(gca,'FontSize',sz)
grid on
xlabel('Time [seconds]')
ylabel('Evolution of variables [-]')
legend(["x","\theta_1","\theta_2","dx/dt"...
    ,"d\theta_1/dt","d\theta_2/dt"],'Location','southeast','Orientation','vertical')

figure(3)
plot(t_ws,uk_ws,"LineWidth",1)
set(gca, 'FontName', 'times')
set(gca,'FontSize',sz)
grid on
xlabel('Time [seconds]')
ylabel('Evolution of u_k')

% 2nd trial
[uk,xk,t] = NR_Dual_Pendulum_Swingup_QTrevised(T,uk_ws);

figure(4)
sz = 12;
plot(t,xk(1,:),"LineWidth",1)
hold on
plot(t,xk(2,:),"LineWidth",1)
plot(t,xk(3,:),"LineWidth",1)
plot(t,xk(4,:),"LineWidth",1)
plot(t,xk(5,:),"LineWidth",1)
plot(t,xk(6,:),"LineWidth",1)
set(gca, 'FontName', 'times')
set(gca,'FontSize',sz)
grid on
xlabel('Time [seconds]')
ylabel('Evolution of variables [-]')
legend(["x","\theta_1","\theta_2","dx/dt"...
    ,"d\theta_1/dt","d\theta_2/dt"],'Location','southeast','Orientation','vertical')

figure(5)
plot(t,uk,"LineWidth",1)
set(gca, 'FontName', 'times')
set(gca,'FontSize',sz)
grid on
xlabel('Time [seconds]')
ylabel('Evolution of u_k')

xk(1:6,end)
%% Q2: computing K(t)

%% Q3: Computing L(t)