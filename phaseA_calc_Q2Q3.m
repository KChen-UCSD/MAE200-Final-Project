clear all;
close all;
clc;

% Load the results from Q.1
load("xk.mat")

steps = size(xk,2);

% Timestep
dT = 0.01;

% Total time
T = (steps-1)*dT;

%% System
% Cart
s.mc = 10; % kg
% Longer pendulum
s.m1 = 1; % kg
s.l1 = 1; % m
% CHECK -------------------------------------------------------------------
% Use half of L1 or L1 ?
s.I1 = (1/12)*s.m1*(s.l1^2);
% Shorter pendulum
s.m2 = 0.5; % kg
s.l2 = 0.5; % m
s.I2 = (1/12)*s.m2*(s.l2^2);

%% Q.2: Computing K(t)
% Constructing E
for i = 1 : steps
    E(:,:,i) = Compute_E(xk(:,i),s);
end

% Constructing A
for i = 1 : steps
    A(:,:,i) = Compute_A(xk(:,i),s);
end

% Constructing B
for i = 1:steps
    B(:,:,i) = [0;0;0;1;0;0];
end

% Constructing C
for i = 1:steps
    C(:,:,i) = [1,0,0,0,0,0;
                0,1,0,0,0,0;
                0,0,1,0,0,0];
end

% Controller: parameters and IC
X_T_con = diag([0 0 0 0 0 0]);
R_con = 1;
Q_con = eye(6);

% Estimator: parameters and IC
X_T_est = diag([0 0 0 0 0 0]); % same as P0
R_est = eye(3);
Q_est = eye(6);


X_ric = DRE(E,A,B,X_T_con,T,R_con,Q_con,1);
P_ric = DRE(E,A,C,X_T_est,T,R_est,Q_est,0);

for i = 1:steps
    K_t(:,:,i) = -inv(R_con)*B(:,:,i)'*X_ric(:,:,i);
end

for i = 1:steps
    L_t(:,:,i) = -P_ric(:,:,i)*C(:,:,i)'*R_est;
end

function E=Compute_E(x,s); I=eye(3); Z=zeros(3);
E=[I Z; Z [s.mc+s.m1+s.m2         -s.m1*s.l1*cos(x(2)) -s.m2*s.l2*cos(x(3));
           -s.m1*s.l1*cos(x(2))  s.I1+s.m1*s.l1^2             0            ;
           -s.m2*s.l2*cos(x(3))          0              s.I2+s.m2*s.l2^2   ]];
end % function Compute_E

function A=Compute_A(x,s); g=9.8;
a42=s.m1*s.l1*(x(8)*sin(x(2))+x(5)^2*cos(x(2))); a45=2*s.m1*s.l1*x(5)*sin(x(2));
a43=s.m2*s.l2*(x(9)*sin(x(3))+x(6)^2*cos(x(3))); a46=2*s.m2*s.l2*x(6)*sin(x(3));
a52=s.m1*s.l1*(g*cos(x(2))-x(7)*sin(x(2))); a63=s.m2*s.l2*(g*cos(x(3))-x(7)*sin(x(3)));
A=[zeros(3) eye(3); 0 -a42 -a43 0 -a45 -a46; 0 a52 0 0 0 0; 0 0 a63 0 0 0];
end % function Compute_A