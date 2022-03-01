clc
clear
clf

r2=1;
r1t2=[];
r2t2=[];
MaxP1t2=[];
MaxP2t2=[];
MaxP3t2=[];

% these are the variables and their initial conditions

%number of cycles to simulate
 numberofcycles=2;

%Iterations
plate_value=1;
ratio=1;
norm_Vb=plate_value/ratio;
norm_Va=norm_Vb;

% center of frame to leg in frame 1
Vb1= -norm_Vb.*[-0.165; 0.17775; 0];
Vb2= -norm_Vb.*[0.165; 0.17775; 0];
Vb3= -norm_Vb.*[0; -0.17775; 0];

% paramters and properties of the frame
height= ratio*2*(abs(Vb3(2)) + abs(Vb1(2)));
base= ratio*(abs(Vb1(1)) + abs(Vb2(2)));

% input the thickness and density
thickness=0.05;
density= 2700;

%total mass of the frame plus the person on it
mass =base*height*thickness*2700/2;

%initialisation of variables
x=0.01;

splines=100000;

omega_Zcenter=4.8/2;

omega_theta=5.2/2;

omega_phi=5/2;

iteration_calculation=min([omega_Zcenter, omega_theta, omega_phi]);

Number_of_iterations=(2*pi()/iteration_calculation)/x;

iterationrounding=round(Number_of_iterations);

if iterationrounding<Number_of_iterations
    Number_of_iterations=iterationrounding+1;
else
    Number_of_iterations=iterationrounding;
end

 Number_of_iterations= Number_of_iterations*numberofcycles;

while r2<=3
r1=1;
r1t=[];
r2t=[];
MaxP1t=[];
MaxP2t=[];
MaxP3t=[];
while r1<=3
    MAX_P1=0;
    MAX_P2=0;
    MAX_P3=0;

% center of frame to leg in frame 0
Va1= -norm_Va.*[r1*0.165; r2*-0.17775; 0]; 
Va2= -norm_Va.*[r1*-0.165; r2*-0.17775; 0];
Va3= -norm_Va.*[0; r2*0.17775; 0];

%matrix initialization
force1=[];
force2=[];
force3=[];
add_on=[];
MAX_Matrix=[];
Length1=[];
Length2=[];
Length3=[];

%time initialization
t=0;
iteration=0;
while iteration<Number_of_iterations
% height, roll and pitch of the frame

Zcenter= 0.05 * cos (t * omega_Zcenter);
Zcenter_dot= 0.05*omega_Zcenter*(-1)*sin(t*omega_Zcenter);
Zcenter_double_dot=  0.05*(omega_Zcenter^2)*(-1)*cos(t*omega_Zcenter);

theta= (pi/6)*cos (t * omega_theta);
theta_dot= (pi/6)*omega_theta*(-1)*sin(t*omega_theta);
theta_double_dot=  (pi/6)*(omega_theta^2)*(-1)*cos(t*omega_theta);

phi= (pi/4)*sin (t * omega_phi);
phi_dot= (pi/4)*omega_phi*(-1)*sin(t*omega_phi);
phi_double_dot= (pi/4)*(omega_phi^2)*(-1)*cos(t*omega_phi);

% center of the frame coordinate
Vp=[0; 0; Zcenter+0.5];

% rotational vector
platform_pitch_matrix= [ 1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta) ];
platform_roll_matrix= [ cos(phi), 0, -sin(phi); 0, 1, 0; sin(phi), 0, cos(phi) ];
Mall= platform_pitch_matrix*platform_roll_matrix;

% leg at frame 0 to leg at frame 1
Displacement_of_link_1= Vp - Mall * Vb1;
Displacement_of_link_2= Vp - Mall * Vb2;
Displacement_of_link_3= Vp - platform_pitch_matrix * Vb3;

% angles of the links
cos_alpha_link_1= -(Displacement_of_link_1(1)-Va1(1))/(sum((Displacement_of_link_1-Va1).^2));
cos_beta_link_1= -(Displacement_of_link_1(2)-Va1(2))/(sum((Displacement_of_link_1-Va1).^2));
cos_delta_link_1= -(Displacement_of_link_1(3)-Va1(3))/(sum((Displacement_of_link_1-Va1).^2));

cos_alpha_link_2= (Displacement_of_link_2(1)-Va2(1))/(sum((Displacement_of_link_2-Va2).^2));
cos_beta_link_2= (Displacement_of_link_2(2)-Va2(2))/(sum((Displacement_of_link_2-Va2).^2));
cos_delta_link_2= (Displacement_of_link_2(3)-Va2(3))/(sum((Displacement_of_link_2-Va2).^2));

cos_alpha_link_3= (Displacement_of_link_3(1)-Va3(1))/(sum((Displacement_of_link_3-Va3).^2));
cos_beta_link_3= (Displacement_of_link_3(2)-Va3(2))/(sum((Displacement_of_link_3-Va3).^2));
cos_delta_link_3= (Displacement_of_link_3(3)-Va3(3))/(sum((Displacement_of_link_3-Va3).^2));

%resultant of forces
force_resultants= [(density*thickness*height*(base^3)/48)*phi_double_dot ; 
                   (density*thickness*(height^3)*base/9)*theta_double_dot; 
                    mass*(Zcenter_double_dot + 9.81)                    ];

%coefficients
coefficients= [cos_delta_link_1*(Displacement_of_link_1(1))+cos_alpha_link_1*(Displacement_of_link_1(3)-Zcenter),cos_delta_link_2*(Displacement_of_link_2(1))+cos_alpha_link_2*(Displacement_of_link_2(3)-Zcenter),0;
               cos_delta_link_1*(Displacement_of_link_1(2))-cos_beta_link_1*(Displacement_of_link_1(3)-Zcenter),cos_delta_link_2*(Displacement_of_link_2(2))-cos_beta_link_2*(Displacement_of_link_2(3)-Zcenter),cos_delta_link_3*(Displacement_of_link_3(2))+cos_beta_link_3*(Displacement_of_link_3(3)-Zcenter);
               cos_delta_link_1,cos_delta_link_2,cos_delta_link_3];
               
%force exerted by leg
forces= coefficients*force_resultants;

%force matrices
force1=[force1, forces(1)];
force2=[force2, forces(2)];
force3=[force3, forces(3)];

% leg speed
L1=sqrt(sum(Displacement_of_link_1.^2));
L2=sqrt(sum(Displacement_of_link_2.^2));
L3=sqrt(sum(Displacement_of_link_3.^2));

Length1=[Length1,L1];
Length2=[Length2,L2];
Length3=[Length3,L3];

MAX_Matrix=[MAX_Matrix,t];
t=t+x;
add_on=[add_on,1];
iteration=iteration+1;
end

Length1=[0, Length1,0];
Length2=[0, Length2,0];
Length3=[0, Length3,0];

MAX_Matrix=MAX_Matrix+4.*add_on;
MAX_Matrix=[0, MAX_Matrix, max(MAX_Matrix)+4];

%Displacement

spline1=spline(MAX_Matrix, [0 Length1 0]);
spline2=spline(MAX_Matrix, [0 Length2 0]);
spline3=spline(MAX_Matrix, [0 Length3 0]);
xq=linspace(min(MAX_Matrix), max(MAX_Matrix), splines);

displacement1matrix= ppval(spline1, xq);
displacement2matrix= ppval(spline2, xq);
displacement3matrix= ppval(spline3, xq);

%velocity

Velocity1=mkpp(MAX_Matrix, [3.*spline1.coefs(:, 1), 2.*spline1.coefs(:, 2), spline1.coefs(:, 3)]);
Velocity2=mkpp(MAX_Matrix, [3.*spline2.coefs(:, 1), 2.*spline2.coefs(:, 2), spline2.coefs(:, 3)]);
Velocity3=mkpp(MAX_Matrix, [3.*spline3.coefs(:, 1), 2.*spline3.coefs(:, 2), spline3.coefs(:, 3)]);

velocity1matrix= ppval(Velocity1, xq);
velocity2matrix= ppval(Velocity2, xq);
velocity3matrix= ppval(Velocity3, xq);

%acceleration

Acceleration1=mkpp(MAX_Matrix, [6.*spline1.coefs(:, 1), 2.*spline1.coefs(:, 2)]);
Acceleration2=mkpp(MAX_Matrix, [6.*spline2.coefs(:, 1), 2.*spline2.coefs(:, 2)]);
Acceleration3=mkpp(MAX_Matrix, [6.*spline3.coefs(:, 1), 2.*spline3.coefs(:, 2)]);

acceleration1matrix= ppval(Acceleration1, xq);
acceleration2matrix= ppval(Acceleration2, xq);
acceleration3matrix= ppval(Acceleration3, xq);

%forces

forcesforspline1=[0, force1,0];
forcesforspline2=[0, force2,0];
forcesforspline3=[0, force3,0];

forcespline1=spline(MAX_Matrix, [0, forcesforspline1, 0]);
forcespline2=spline(MAX_Matrix, [0, forcesforspline2, 0]);
forcespline3=spline(MAX_Matrix, [0, forcesforspline3, 0]);

force1matrix= ppval(forcespline1, xq);
force2matrix= ppval(forcespline2, xq);
force3matrix= ppval(forcespline3, xq);

%power

power1matrix= abs(force1matrix.*velocity1matrix);
power2matrix= abs(force2matrix.*velocity2matrix);
power3matrix= abs(force3matrix.*velocity3matrix);

MAX_P1=max(power1matrix);
MAX_P2=max(power2matrix);
MAX_P3=max(power3matrix);

r1t=[r1t;r1];

r2t=[r2t;r2];

MaxP1t=[MaxP1t; MAX_P1];
MaxP2t=[MaxP2t; MAX_P2];
MaxP3t=[MaxP3t; MAX_P3];

r1=r1+0.03;
end
r1t2=[r1t2;r1t];

r2t2=[r2t2;r2t];

MaxP1t2=[MaxP1t2; MaxP1t];
MaxP2t2=[MaxP2t2; MaxP2t];
MaxP3t2=[MaxP3t2; MaxP3t];

r2=r2+0.03;
end
Powermatrix= MaxP1t2+MaxP2t2+MaxP3t2;
figure(1)
plot3(r1t2,r2t2,Powermatrix)
minPower=min(Powermatrix);

minPowerindex= find(Powermatrix==minPower);

optimumr1value= r1t2(minPowerindex)

optimumr2value= r2t2(minPowerindex)

MATRIXPOWER= [r1t2, Powermatrix, r2t2];

r1= optimumr1value;

r2= optimumr2value;

%step 2:simulating using the optimal values

% center of frame to leg in frame 0
Va1= -norm_Va.*[r1*0.165; r2*-0.17775; 0]; 
Va2= -norm_Va.*[r1*-0.165; r2*-0.17775; 0];
Va3= -norm_Va.*[0; r2*0.17775; 0];

%matrix initialization
force1=[];
force2=[];
force3=[];
add_on=[];
MAX_Matrix=[];
Length1=[];
Length2=[];
Length3=[];

%time initialization
t=0;

Number_of_iterations=Number_of_iterations*100;

x=x/100;

iteration=0;

while iteration<Number_of_iterations
% height, roll and pitch of the frame


Zcenter= 0.05 * sin (t * omega_Zcenter);
Zcenter_dot= 0.05*omega_Zcenter*cos(t*omega_Zcenter);
Zcenter_double_dot=  0.05*(omega_Zcenter^2)*(-1)*sin(t*omega_Zcenter);

theta= (pi/6)*sin (t * omega_theta);
theta_dot= (pi/6)*omega_theta*cos(t*omega_theta);
theta_double_dot=  (pi/6)*(omega_theta^2)*(-1)*sin(t*omega_theta);

phi= (pi/4)*sin (t * omega_phi);
phi_dot= (pi/4)*omega_phi*cos(t*omega_phi);
phi_double_dot= (pi/4)*(omega_phi^2)*(-1)*sin(t*omega_phi);

% center of the frame coordinate
Vp=[0; 0; Zcenter+0.5775];

% rotational vector
platform_pitch_matrix= [ 1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta) ];
platform_roll_matrix= [ cos(phi), 0, -sin(phi); 0, 1, 0; sin(phi), 0, cos(phi) ];
Mall= platform_pitch_matrix*platform_roll_matrix;

% leg at frame 0 to leg at frame 1
Displacement_of_link_1= Vp - Mall * Vb1;
Displacement_of_link_2= Vp - Mall * Vb2;
Displacement_of_link_3= Vp - platform_pitch_matrix * Vb3;

% angles of the links
cos_alpha_link_1= -(Displacement_of_link_1(1)-Va1(1))/(sum((Displacement_of_link_1-Va1).^2));
cos_beta_link_1= -(Displacement_of_link_1(2)-Va1(2))/(sum((Displacement_of_link_1-Va1).^2));
cos_delta_link_1= -(Displacement_of_link_1(3)-Va1(3))/(sum((Displacement_of_link_1-Va1).^2));

cos_alpha_link_2= (Displacement_of_link_2(1)-Va2(1))/(sum((Displacement_of_link_2-Va2).^2));
cos_beta_link_2= (Displacement_of_link_2(2)-Va2(2))/(sum((Displacement_of_link_2-Va2).^2));
cos_delta_link_2= (Displacement_of_link_2(3)-Va2(3))/(sum((Displacement_of_link_2-Va2).^2));

cos_alpha_link_3= (Displacement_of_link_3(1)-Va3(1))/(sum((Displacement_of_link_3-Va3).^2));
cos_beta_link_3= (Displacement_of_link_3(2)-Va3(2))/(sum((Displacement_of_link_3-Va3).^2));
cos_delta_link_3= (Displacement_of_link_3(3)-Va3(3))/(sum((Displacement_of_link_3-Va3).^2));

%resultant of forces
force_resultants= [(density*thickness*height*(base^3)/48)*phi_double_dot ; 
                   (density*thickness*(height^3)*base/9)*theta_double_dot; 
                    mass*(Zcenter_double_dot + 9.81)                    ];

%coefficients
coefficients= [cos_delta_link_1*(Displacement_of_link_1(1))+cos_alpha_link_1*(Displacement_of_link_1(3)-Zcenter),cos_delta_link_2*(Displacement_of_link_2(1))+cos_alpha_link_2*(Displacement_of_link_2(3)-Zcenter),0;
               cos_delta_link_1*(Displacement_of_link_1(2))-cos_beta_link_1*(Displacement_of_link_1(3)-Zcenter),cos_delta_link_2*(Displacement_of_link_2(2))-cos_beta_link_2*(Displacement_of_link_2(3)-Zcenter),cos_delta_link_3*(Displacement_of_link_3(2))+cos_beta_link_3*(Displacement_of_link_3(3)-Zcenter);
               cos_delta_link_1,cos_delta_link_2,cos_delta_link_3];
               
%force exerted by leg
forces= coefficients*force_resultants;

%force matrices
force1=[force1, forces(1)];
force2=[force2, forces(2)];
force3=[force3, forces(3)];

% leg speed
L1=sqrt(sum(Displacement_of_link_1.^2));
L2=sqrt(sum(Displacement_of_link_2.^2));
L3=sqrt(sum(Displacement_of_link_3.^2));

Length1=[Length1,L1];
Length2=[Length2,L2];
Length3=[Length3,L3];

MAX_Matrix=[MAX_Matrix,t];
t=t+x;
add_on=[add_on,1];
iteration=iteration+1;
end

Length1=[1, Length1,1];
Length2=[1, Length2,1];
Length3=[1, Length3,1];

MAX_Matrix=MAX_Matrix+1.*add_on;
MAX_Matrix=[0, MAX_Matrix, max(MAX_Matrix)+1];

%Displacement

spline1=spline(MAX_Matrix, [0 Length1 0]);
spline2=spline(MAX_Matrix, [0 Length2 0]);
spline3=spline(MAX_Matrix, [0 Length3 0]);
xq=linspace(min(MAX_Matrix), max(MAX_Matrix), splines);

displacement1matrix= ppval(spline1, xq);
displacement2matrix= ppval(spline2, xq);
displacement3matrix= ppval(spline3, xq);

figure (2)
plot(xq, displacement1matrix)
figure (3)
plot(xq, displacement2matrix)
figure (4)
plot(xq, displacement3matrix)

%velocity

Velocity1=mkpp(MAX_Matrix, [3.*spline1.coefs(:, 1), 2.*spline1.coefs(:, 2), spline1.coefs(:, 3)]);
Velocity2=mkpp(MAX_Matrix, [3.*spline2.coefs(:, 1), 2.*spline2.coefs(:, 2), spline2.coefs(:, 3)]);
Velocity3=mkpp(MAX_Matrix, [3.*spline3.coefs(:, 1), 2.*spline3.coefs(:, 2), spline3.coefs(:, 3)]);

velocity1matrix= ppval(Velocity1, xq);
velocity2matrix= ppval(Velocity2, xq);
velocity3matrix= ppval(Velocity3, xq);

figure (5)
plot(xq, velocity1matrix)
figure (6)
plot(xq, velocity2matrix)
figure (7)
plot(xq, velocity3matrix)

%acceleration

Acceleration1=mkpp(MAX_Matrix, [6.*spline1.coefs(:, 1), 2.*spline1.coefs(:, 2)]);
Acceleration2=mkpp(MAX_Matrix, [6.*spline2.coefs(:, 1), 2.*spline2.coefs(:, 2)]);
Acceleration3=mkpp(MAX_Matrix, [6.*spline3.coefs(:, 1), 2.*spline3.coefs(:, 2)]);

acceleration1matrix= ppval(Acceleration1, xq);
acceleration2matrix= ppval(Acceleration2, xq);
acceleration3matrix= ppval(Acceleration3, xq);

figure (8)
plot(xq, acceleration1matrix)
figure (9)
plot(xq, acceleration2matrix)
figure (10)
plot(xq, acceleration3matrix)

%forces

forcesforspline1=[0, force1,0];
forcesforspline2=[0, force2,0];
forcesforspline3=[0, force3,0];

forcespline1=spline(MAX_Matrix, [0, forcesforspline1, 0]);
forcespline2=spline(MAX_Matrix, [0, forcesforspline2, 0]);
forcespline3=spline(MAX_Matrix, [0, forcesforspline3, 0]);

force1matrix= ppval(forcespline1, xq);
force2matrix= ppval(forcespline2, xq);
force3matrix= ppval(forcespline3, xq);

figure (11)
plot(xq, force1matrix)
figure (12)
plot(xq, force2matrix)
figure (13)
plot(xq, force3matrix)

%power

power1matrix= abs(force1matrix.*velocity1matrix);
power2matrix= abs(force2matrix.*velocity2matrix);
power3matrix= abs(force3matrix.*velocity3matrix);

figure (14)
plot(xq, power1matrix)
figure (15)
plot(xq, power2matrix)
figure (16)
plot(xq, power3matrix)

rotationalvelocity1=(100/2.31).*velocity1matrix;
rotationalvelocity2=(100/2.31).*velocity2matrix;
rotationalvelocity3=(100/2.31).*velocity3matrix;

mazen=[rotationalvelocity1; rotationalvelocity2; rotationalvelocity3; xq];