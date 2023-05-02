function [dX,Tau] = ode_rrbot(t, X, k, P, gamma, B)
    
    %Initiating Acual value parameters of the Robot
    m1 = 1; m2 = 1; l1 = 1; l2 = 1;r1 = 0.45;r2 = 0.45;I1 = 0.084; I2 = 0.084;g = 9.81;
    
    
    %Initialising Nominal Values 
    %m1_hat = 0.75; m2_hat = 0.75; I1_hat = 0.063; I2_hat = 0.063;

    dX=zeros(9,1);  % create empty vector for [theta1_dot; theta1_ddot; theta2_dot; theta2_ddot]
    X=num2cell(X);
    [theta1,theta2,theta1_dot,theta2_dot,alpha1, alpha2, alpha3, alpha4, alpha5]=deal(X{:});
    
    %Initiating Nominal values 
    alpha_hat = [ alpha1; alpha2; alpha3; alpha4; alpha5];
    


     %Initialising desired Trajectories 
    q1_d = (pi*t^3)/500 - (3*pi*t^2)/100 + pi;
    q1_dot_d = (3*pi*t^2)/500 - (3*pi*t)/50;
    q1_ddot_d = (3*pi*t)/250 - (3*pi)/50;
    q2_d = (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
    q2_dot_d = (3*pi*t^2)/1000 - (3*pi*t)/100;
    q2_ddot_d = (3*pi*t)/500 - (3*pi)/100;

    %Desired State Trajectories
    x_d = [q1_d; q2_d; q1_dot_d;q2_dot_d];
    v_d = [q1_ddot_d;q2_ddot_d];
    
    %Current State Trajectories
    x = [theta1;theta2;theta1_dot;theta2_dot];
   
    e = x - x_d;

    %Initialising torque inputs for Robust Control 
    %Adding Gain Matrix
    %k = [2 0 3 0; 0 2 0 3];
    %k = [12 0 7 0; 0 12 0 7];

    %Virtual Input Term 
    v = -k*(e) + v_d;

    %Y_t defined by substituting v (virtual Control Input) 

    Y_t = [v(1), ...
    cos(theta2)*(2*v(1) + v(2)) - 2*sin(theta2)*theta1_dot*theta2_dot - sin(theta2)*theta2_dot^2, ...
    v(2), ...
    -sin(theta1)*g, ...
    -sin(theta1 + theta2)*g; ...
    0, ...
    sin(theta2)*theta1_dot^2 + cos(theta2)*v(1), ...
    v(1) + v(2), ...
    0, ...
    -sin(theta1+theta2)*g];

    Tau = Y_t*alpha_hat;

    %Final Control Input Term
    tau1 = Tau(1);
    tau2 = Tau(2);
    

    %State space representation augmented with alpha values 
    dX(1)= theta1_dot;
    dX(2)= theta2_dot;
    dX(3)= (I2*tau1 - I2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
    dX(4)= -(I2*tau1 - I1*tau2 - I2*tau2 - l1^2*m2*tau2 - m1*r1^2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*tau1*cos(theta2) - 2*l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
    theta1_ddot = dX(3);
    theta2_ddot = dX(4);

    %Derived form Y 
    m_1 = [1, 2*cos(theta2), 0, 0, 0;
      0, cos(theta2), 1, 0, 0];
    m_2 = [0, cos(theta2), 1, 0, 0;
      0, 0, 1, 0, 0];

    Y = [theta1_ddot, ...
    cos(theta2)*(2*theta1_ddot + theta2_ddot) - 2*sin(theta2)*theta1_dot*theta2_dot - sin(theta2)*theta2_dot^2, ...
    theta2_ddot, ...
    -sin(theta1)*g, ...
    -sin(theta1 + theta2)*g; ...
    0, ...
    sin(theta2)*theta1_dot^2 + cos(theta2)*theta1_ddot, ...
    theta1_ddot + theta1_ddot, ...
    0, ...
    -sin(theta1+theta2)*g];

    M_hat = [m_1*alpha_hat, m_2*alpha_hat];

    phi = M_hat\Y;

    alpha_hat_dot = -gamma\(phi'*B'*P*e);

    dX(5) = alpha_hat_dot(1);
    dX(6) = alpha_hat_dot(2);
    dX(7) = alpha_hat_dot(3);
    dX(8) = alpha_hat_dot(4);
    dX(9) = alpha_hat_dot(5);

    display(t)
 
end
