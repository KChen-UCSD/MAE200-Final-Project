clear all;
close all;
clc;

%% Inputting parameters
% Cart
mc = 10; % kg
% Longer pendulum
m1 = 1; % kg
l1 = 1; % m
% CHECK -------------------------------------------------------------------
% Use half of L1 or L1 ?
I1 = (1/12)*m1*(l1^2);
% Shorter pendulum
m2 = 0.5; % kg
l2 = 0.5; % m
I2 = (1/12)*m2*(l2^2);
g = 9.8; % m/s^2

% Definition of x
% x = [x;theta1;theta2;x_dot;theta1_dot;theta2_dot]

% In the form of E * x_dot = A_star * x + B_star * u
E = [1,0,0,0,0,0;
     0,1,0,0,0,0;
     0,0,1,0,0,0;
     0,0,0,(mc+m1+m2),-m1*l1,-m2*l2;
     0,0,0,-m1*l1,(I1+m1*(l1^2)),0;
     0,0,0,-m2*l2,0,(I2+m2*(l2^2))];

A_star = [0,0,0,1,0,0;
          0,0,0,0,1,0;
          0,0,0,0,0,1;
          0,0,0,0,0,0;
          0,m1*g*l1,0,0,0,0;
          0,0,m2*g*l2,0,0,0];

B_star = [0;0;0;1;0;0];

% In the form of x_dot = A * x + B * u
A = inv(E)*A_star;
B = inv(E)*B_star;

% As given
C = [1,0,0,0,0,0;
     0,1,0,0,0,0;
     0,0,1,0,0,0];

% Size of the system
size_matrix = length(B_star);

%% Q4. Feedback control: determining an optimal K
% Assume a Q
Q_K = eye(size_matrix);
% Assume a R
R_K = 1;

% After inversing, E = I on the LHS
E1 = eye(size_matrix);

% Using icare to solve for X
X = icare(A,B,Q_K,R_K,[],E1,[]);

% K can be computed as X is calculated
K = -inv(R_K)*B'*X;

% Checking that K works
x0 = [0.0114
    0.7622
    0.8823
   -0.6256
   -0.8477
   -1.4435];
% x0 = [.2; pi/12; pi/12; .1; .2; .2]; 
dt = .01; x = x0;

steps = 3000;
num = zeros(steps,1);
time = zeros(steps,1);
x_1 = zeros(steps,1);
x_2 = zeros(steps,1);
x_3 = zeros(steps,1);
x_4 = zeros(steps,1);
x_5 = zeros(steps,1);
x_6 = zeros(steps,1);

for i = 1 : steps

    x_1(i) = x(1); x_2(i) = x(2); x_3(i) = x(3);
    x_4(i) = x(4); x_5(i) = x(5); x_6(i) = x(6);

    num(i) = i;
    time(i) = dt*i;

    u = K * x;
    f1 = A * x + B * u;

    x1 = x + dt/2 * f1;
    u1 = K * x1;
    f2 = A * x1 + B * u1;

    x2 = x + dt/2 * f2;
    u2 = K * x2;
    f3 = A * x2 + B * u2;

    x3 = x + dt * f3;
    u3 = K * x3;
    f4 = A * x3 + B * u3;

    x = x + dt * (f1/6 + (f2+f3)/3 + f4/6);

end

figure(1)
sz = 12;
plot(time,x_1,"LineWidth",1)
hold on
plot(time,x_2,"LineWidth",1)
plot(time,x_3,"LineWidth",1)
plot(time,x_4,"LineWidth",1)
plot(time,x_5,"LineWidth",1)
plot(time,x_6,"LineWidth",1)
set(gca, 'FontName', 'times')
set(gca,'FontSize',sz)
grid on
xlabel('Time [seconds]')
ylabel('Evolution of variables [-]')
legend(["x","\theta_1","\theta_2","dx/dt"...
    ,"d\theta_1/dt","d\theta_2/dt"],'Location','southeast','Orientation','vertical')

figure(3)
plot(time,x_2,"LineWidth",1)
hold on
plot(time,x_3,"LineWidth",1)
set(gca, 'FontName', 'times')
set(gca,'FontSize',sz)
grid on
xlabel('Time [seconds]')
ylabel('Evolution of variables [-]')
legend(["\theta_1","\theta_2"],'Location','northeast','Orientation','vertical')

%% Fine-tuning of K: Demonstration

% Fine-tuning R
R_K = 0.1;

% Using icare to solve for X
X = icare(A,B,Q_K,R_K,[],E1,[]);

% K can be computed as X is calculated
K_tuned = -inv(R_K)*B'*X;

% Checking that K works
x0 = [.2; pi/12; pi/12; .1; .2; .2]; 
dt = .01; x = x0;

steps = 3000;
num = zeros(steps,1);
time = zeros(steps,1);
x_1 = zeros(steps,1);
x_2 = zeros(steps,1);
x_3 = zeros(steps,1);
x_4 = zeros(steps,1);
x_5 = zeros(steps,1);
x_6 = zeros(steps,1);

for i = 1 : steps

    x_1(i) = x(1); x_2(i) = x(2); x_3(i) = x(3);
    x_4(i) = x(4); x_5(i) = x(5); x_6(i) = x(6);

    num(i) = i;
    time(i) = dt*i;

    u = K_tuned * x;
    f1 = A * x + B * u;

    x1 = x + dt/2 * f1;
    u1 = K_tuned * x1;
    f2 = A * x1 + B * u1;

    x2 = x + dt/2 * f2;
    u2 = K_tuned * x2;
    f3 = A * x2 + B * u2;

    x3 = x + dt * f3;
    u3 = K_tuned * x3;
    f4 = A * x3 + B * u3;

    x = x + dt * (f1/6 + (f2+f3)/3 + f4/6);

end

figure(2)
sz = 12;
plot(time,x_1,"LineWidth",1)
hold on
plot(time,x_2,"LineWidth",1)
plot(time,x_3,"LineWidth",1)
plot(time,x_4,"LineWidth",1)
plot(time,x_5,"LineWidth",1)
plot(time,x_6,"LineWidth",1)
set(gca, 'FontName', 'times')
set(gca,'FontSize',sz)
grid on
xlabel('Time [seconds]')
ylabel('Evolution of variables [-]')
legend(["x","\theta_1","\theta_2","dx/dt"...
    ,"d\theta_1/dt","d\theta_2/dt"],'Location','southeast','Orientation','vertical')

%% Q5. State Estimation: determining an optimal K
% Assume a Q
Q_L = eye(size_matrix);

% Assume a R
R_L = 1;

A_use = transpose(A);
B_use = transpose(C);

% Using icare to solve for X
P = icare(A_use,B_use,Q_L,R_L,[],E1,[]);

% L can be computed as X is calculated
L = -P*transpose(C)*inv(R_L);

% norm(eig(A+L*C))
% norm(eig(A+B*K))
% eig(A+L*C)
% eig(A+B*K)

%% Fine-tuning of K: Demonstration

% Assume a Q
Q_L = [100,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        0,0,0,1,0,0;
        0,0,0,0,1,0;
        0,0,0,0,0,1];

% Assume a R
R_L = 0.2;

% Using icare to solve for X
P = icare(A_use,B_use,Q_L,R_L,[],E1,[]);

% L can be computed as X is calculated
L_tuned = -P*transpose(C)*inv(R_L);