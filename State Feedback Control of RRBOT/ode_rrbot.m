function dX = ode_rrbot(t,X)
    
   %initiating values of parameters 
    m1 = 1; m2 = 1; l1 = 1; l2 = 1;r1 = 0.45;r2 = 0.45;I1 = 0.084; I2 = 0.084;g = 9.81; T = 10;

    dX=zeros(4,1);  % create empty vector for [theta1_dot; theta1_ddot; theta2_dot; theta2_ddot]
    X=num2cell(X);
    [theta1,theta2,theta1_dot,theta2_dot]=deal(X{:});

    %tau1=0;tau2=0;
    %initialising torque inputs 
    %Adding Gain Matrix
    k11 = 23.5850;
    k12 = 5.8875;
    k13 = 5.1470;
    k14 = 2.6108;
    k21 = 5.8875;
    k22 = 4.9875;
    k23 = 1.5443;
    k24 = 0.9770;
    
    
    K =[k11,k12,k13,k14; k21,k22,k23,k24];

    X = zeros(4,1);
    X(1,1)= theta1;
    X(2,1)= theta2;
    X(3,1)= theta1_dot;
    X(4,1)= theta2_dot;

    %tau1 =  - (k11 *theta1 +k12*theta2 + k13 *theta1_dot + k14*theta2_dot);
    %tau2 =  - (k21 *theta1 +k22*theta2 + k23 *theta1_dot + k24*theta2_dot);

    u= -K*X;
    
    tau1 = u(1,1);
    tau2 = u(2,1);
 
    %State space representation
    %pasting values from the solve function
    dX(1)= theta1_dot;
    dX(2)= theta2_dot;
    dX(3)= (I2*tau1 - I2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
    dX(4)= -(I2*tau1 - I1*tau2 - I2*tau2 - l1^2*m2*tau2 - m1*r1^2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*tau1*cos(theta2) - 2*l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
 
 

