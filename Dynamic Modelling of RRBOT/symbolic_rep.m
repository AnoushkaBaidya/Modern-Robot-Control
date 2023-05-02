clear;clc;

%RRBot system 2DOF 
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot F t g  l1 l2 r1 r2 I1 I2 tau1 tau2 'real'
syms m1 m2 'real'

X = sym('X', [4,1]);
X(1)= theta1;
X(2)= theta2;
X(3)= theta1_dot;
X(4)= theta2_dot;

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
%display(eq1)
%display(eq2)

sol = solve([eqn1==0, eqn2==0], [theta1_ddot, theta2_ddot]);

display(sol.theta1_ddot)
display(sol.theta2_ddot)



