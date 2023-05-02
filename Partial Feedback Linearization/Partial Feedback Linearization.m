%This is the main function
clear;clc;

%RRBot system 2DOF 
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot F t g  l1 l2 r1 r2 I1 I2 tau1 tau2  'real'
syms m1 m2 'real'
%Declaring Symbolic Representation Variables
syms X1 X2 X3 X4 u lambda 'real'
syms t0 tf q0 q0_dot qf qf_dot 'real'

%Calculate the Langrangian Equation L=K-P
K = (1/2)*m1*(r1^2)*(theta1_dot^2) + (1/2)*I1*(theta1_dot^2) + (1/2)*m2*(l1^2)*(theta1_dot^2) + (1/2)*m2*(r2^2)*((theta2_dot+theta1_dot)^2) + m2*l1*r2*theta1_dot*(theta1_dot+theta2_dot)*cos(theta2) +(1/2)*I2*((theta2_dot+theta1_dot)^2);
%display(K)

P = m1*g*r1*cos(theta1) + m2*g*l1*cos(theta1) + m2*g*r2*cos(theta1 + theta2);
%display(P)

%Lagrangian Equation 
L=K-P;
%display(L)

%Solving for Step4 
dL_dtheta1 = jacobian(L,theta1);
dL_dtheta2 = jacobian(L,theta2);

dL_dtheta1_dot = jacobian(L, theta1_dot);
dL_dtheta2_dot = jacobian(L, theta2_dot);

ddt_dL_dtheta1_dot = jacobian(dL_dtheta1_dot,[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];
ddt_dL_dtheta2_dot = jacobian(dL_dtheta2_dot,[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];
%display(ddt_dL_dtheta1_dot);
%display(ddt_dL_dtheta2_dot);

%Equations of Motion 
eqn1=ddt_dL_dtheta1_dot - dL_dtheta1 -tau1;
eqn2=ddt_dL_dtheta2_dot - dL_dtheta2 -tau2;

%display(eqn1)
%display(eqn2)

%A]Generate cubic polynomial Trajectory 
%Joint1(theta1) and Joint2(theta2)

T = [1 t0 t0^2 t0^3; 0 1 2*t0 3*t0^2; 1 tf tf^2 tf^3; 0 1 2*tf 3*tf^2];
q = [q0;q0_dot;qf;qf_dot];

% a0,a1,a2,a3 
A = inv(T)*q;
A_inv = T\q; 

%Declaring the initial and final positions and velocity conditions
% Joint 1 -> [to,tf] =[0,10] [q0;q0_dot;qf;qf_dot] = [180,0,0,0]
joint_1 = subs(A,[t0,tf,q0,q0_dot,qf,qf_dot],[0,10,deg2rad(180),0,0,0]);
% Joint 1 -> [to,tf] =[0,10] [q0;q0_dot;qf;qf_dot] = [180,0,0,0]
joint_2 = subs(A,[t0,tf,q0,q0_dot,qf,qf_dot],[0,10,deg2rad(90),0,0,0]);

display(joint_1);
display(joint_2);

%joint1 q1_desired(theta1_desired) trajectory values
q1_d = joint_1(1) + joint_1(2)*t + joint_1(3)*t^2 + joint_1(4)*t^3;
q1_dot_d = jacobian(q1_d,t);
q1_ddot_d = jacobian(q1_dot_d,t);

display(q1_d);
display(q1_dot_d);
display(q1_ddot_d);

%joint2 q2_desired(theta2_desired) trajectory values
q2_d = joint_2(1) + joint_2(2)*t + joint_2(3)*t^2 + joint_2(4)*t^3;
q2_dot_d = jacobian(q2_d,t);
q2_ddot_d = jacobian(q2_dot_d,t);

display(q2_d);
display(q2_dot_d);
display(q2_ddot_d);

%q_d values and v_d values 
q_d = [q1_d; q2_d; q1_dot_d;q2_dot_d];
v_d = [q1_ddot_d;q2_ddot_d];

%B] Standard Manipulator Equation Form 
eq1 = eqn1 +tau1;   %eq1 -> M(q)q_ddot + C(q,q_dot)q_dot + G(q)
eq2  = eqn2 +tau2;  %eq2 -> M(q)q_ddot + C(q,q_dot)q_dot + G(q)

%Solve for G(q) matrix by subs q and q_ddot values = 0
g1 = subs(eq1,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot],[0,0,0,0]);
g2 = subs(eq2,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot],[0,0,0,0]);
%G(q) Gravity Vector 
G = [g1;g2];
display(g1)
display(g2)

%Solving for M(q) matrix 
m11 = subs(eq1 - g1, [theta1_dot,theta2_dot,theta2_ddot],[0,0,0])/ theta1_ddot;
m11 = simplify(m11);
m12 = subs(eq1 - g1, [theta1_dot,theta2_dot,theta1_ddot],[0,0,0])/ theta2_ddot;
m12 = simplify(m12);
m21 = subs(eq2 - g2, [theta1_dot,theta2_dot,theta2_ddot],[0,0,0])/ theta1_ddot;
m21 = simplify(m21);
m22 = subs(eq2 - g2, [theta1_dot,theta2_dot,theta1_ddot],[0,0,0])/ theta2_ddot;
m22 = simplify(m22);

%M(q) Mass Matrix
M = [m11 m12; m21 m22;];
display(M)


%Solving for C(q,q_dot) matrix

c1 = (eq1 -m11*theta1_ddot - m12*theta2_ddot - g1 );
c1 = simplify(c1);
c2 = (eq2 -m21*theta1_ddot - m22*theta2_ddot - g2 );
c2 = simplify(c2);

%C(q,q_dot) Coriolis Matrix 
C = [c1;c2];
display(C)


%Substitute the parameters with the actual value 
% [m1,m2,l1,l2,r1,r2,I1,I2,g] -> [1,1,1,1,0.45,0.45,0.084,0.084,9.81]

M = subs(M,[m1,m2,l1,l2,r1,r2,I1,I2,g],[1,1,1,1,0.45,0.45,0.084,0.084,9.81]);
C = subs(C, [m1,m2,l1,l2,r1,r2,I1,I2,g],[1,1,1,1,0.45,0.45,0.084,0.084,9.81]);
G = subs(G, [m1,m2,l1,l2,r1,r2,I1,I2,g],[1,1,1,1,0.45,0.45,0.084,0.084,9.81]);

display(M);
display(C);
display(G);

% Controller Design for our System 
% C] Symbolic Feedback Linearization 
% v = q_ddot -> v1 = q1_ddot ; v2 = q2_ddot
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];

%Step1 :Check the controllability of the system 
Rank_Matrix = rank(ctrb(A,B)); 
%Rank = 4 which is greater than n thus Controllable System 

% State Feedback Design for Eigen Values of my choice 

%lambda = [-1,-2,-1+1i, -1-1i];
lambda = [-10,-5,-10,-5];
Kn = place(A,B,lambda);

sz_K = size(Kn); % kn = 2X4

%Gain Matrix;
display(Kn);



%Step3 of Designing Controller 
%Controller Design
% e = x-x_d -> [q-q_d; q_dot- q_dotd]
% v = -K*e + vd
% tau = M(q)*v + C(q,q_dot) + G(q)

%E] Simulate it forward in time using ODE 45

%Simulation time set for 10 seconds
T=10;

%Setting initial conditions of the system 
y0 = [deg2rad(200), deg2rad(125),0,0];   % [theta_1 , theta_2, theta1_dot, thteta2_dot]

%Calling the ode45 function 
%Simulate Using ode45
[t,y] = ode45(@ode_rrbot,[0,T],y0);

%Reconstruct Control Input 
%This we need to store the tau values after simulating it forward in Time
tau1_list=[];
tau2_list=[];
theta1_d= [];
theta2_d=[];
theta1_dot_d =[];
theta2_dot_d =[];

 for i = 1:size(t)    
    [xi,tau1,tau2,q_d]= ode_rrbot(t(i),y(i));
    tau1_list = [tau1_list,tau1];
    tau2_list = [tau2_list,tau2];
    %trajectories 
    q1_d = q_d(1,1);
    q2_d = q_d(2,1);
    q1_dot_d = q_d(3,1);
    q2_dot_d = q_d(4,1);
    theta1_d = [theta1_d,q1_d];
    theta2_d = [theta2_d,q2_d];
    theta1_dot_d = [theta1_dot_d,q1_dot_d];
    theta2_dot_d = [theta2_dot_d,q2_dot_d];

 end

 tau_list = [];
for i = 1:size(t)
    M_sub = double(subs(M,[theta1,theta2,theta1_dot,theta2_dot],[y(i,1),y(i,2),y(i,3),y(i,4)]));
    C_sub = double(subs(C,[theta1,theta2,theta1_dot,theta2_dot],[y(i,1),y(i,2),y(i,3),y(i,4)]));
    G_sub = double(subs(G,[theta1,theta2,theta1_dot,theta2_dot],[y(i,1),y(i,2),y(i,3),y(i,4)]));

    q1_d = (pi*t(i)^3)/500 - (3*pi*t(i)^2)/100 + pi;
    q2_d = (pi*t(i)^3)/1000 - (3*pi*t(i)^2)/200 + pi/2;
    q1_dot_d = (3*pi*t(i)^2)/500 - (3*pi*t(i))/50;
    q2_dot_d = (3*pi*t(i)^2)/1000 - (3*pi*t(i))/100;
    q1_ddot_d = (3*pi*t(i))/250 - (3*pi)/50;
    q2_ddot_d = (3*pi*t(i))/500 - (3*pi)/100;
    xd = [q1_d;q2_d;q1_dot_d;q2_dot_d];
    vd = [q1_ddot_d;q2_ddot_d];
    
    x = [y(i,1);y(i,2);y(i,3);y(i,4)];
    v = -Kn*(x - xd) + vd;

    Tau = M_sub*v + C_sub+ G_sub;
    Tau = double(Tau);
    tau_list = [tau_list Tau];
end


%Plot the Data %visualize the output 
figure;
text(0.5, 100, 'Blue: State Trajectory and Red: Desired State Trajectory', 'HorizontalAlignment', 'center', 'FontSize', 12);
subplot(2,2,1);
plot(t,rad2deg(y(:,1)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1 in deg','FontSize',10);
hold 'on';
plot(t,rad2deg(theta1_d(1,:)),'r','linewidth',2);

subplot(2,2,2);
plot(t,rad2deg(y(:,2)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2 in deg','FontSize',10);
hold 'on';
plot(t,rad2deg(theta2_d(1,:)),'r','linewidth',2);

subplot(2,2,3);
plot(t,rad2deg(y(:,3)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1dot in degrees','FontSize',10);
hold 'on';
plot(t,rad2deg(theta1_dot_d(1,:)),'r','linewidth',2);

subplot(2,2,4);
plot(t,rad2deg(y(:,4)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2dot in degrees','FontSize',10);
hold 'on';
plot(t,rad2deg(theta2_dot_d(1,:)),'r','linewidth',2);

figure;
subplot(2,1,1);
plot(t,tau_list(1,:),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau1','FontSize',10);
subplot(2,1,2);
plot(t,tau_list(2,:),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau2','FontSize',10);














