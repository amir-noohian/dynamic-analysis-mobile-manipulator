clc
clear

syms t
syms q1 q2 q3 q4 q5 q6 q7 q8 q9 % q1 = X, q2 = Y, q3 = Z, q4 = psi, q5 = theta, q6 = phi, q7 = theta1, q8 = theta2, q9 = theta3
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 dq9 % dq1 = dX/dt, ...
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 ddq8 ddq9 % ddq1 = d^2X/dt^2, ...
syms ms m1 m2 m3 a b c L1 L2 L3 L4 Radius tau1 tau2 tau3

q = [q1; q2; q3; q4; q5; q6; q7; q8; q9]; % The vector of generalized coordinates
dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7; dq8; dq9]; % The vector of generalized velocities
ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6; ddq7; ddq8; ddq9]; % The vector of generalized accelerations

%% Evaluating the angular momentum and the kinetic energy of the satellite (without the manipulator)

VG_XYZ = [dq1; dq2; dq3]; % The velocity of the center of mass of the satellite in the inertia coordinates.

% Psi, theta, and phi are the rotation angles from the inertia coordinates to the body coordinates
RZ = [cos(q4) sin(q4) 0;-sin(q4) cos(q4) 0; 0 0 1];
RY = [cos(q5) 0 -sin(q5);0 1 0;sin(q5) 0 cos(q5)];
RX = [1 0 0;0 cos(q6) sin(q6);0 -sin(q6) cos(q6)];
R = RX*RY*RZ; % It converts coordinates from the inertia to the body one.

ws_xyz = RX*RY*[0; 0; dq4] + RX*[0; dq5; 0] +  [dq6; 0; 0]; % The angular velocity of the satellite.
Ixx = 1/12*ms*(b^2+c^2);
Iyy = 1/12*ms*(c^2+a^2);
Izz = 1/12*ms*(a^2+b^2);
Is_xyz = [Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % The inertia matrix of the satellite in the body coordinates.
Hs_xyz = Is_xyz*ws_xyz; % Angular momentum of the satellite in body coordinates.

% Kinetic energy of the satellite
Ts = 1/2*ms*(VG_XYZ.'*VG_XYZ) + 1/2*ws_xyz.'*Hs_xyz;

%% Evaluating the angular momentum and the kinetic energy of the arms of the manipulator.

% These rotation matrices transfer coordinates from the new one to old one.
R1 = [cos(q7) -sin(q7) 0; sin(q7) cos(q7) 0; 0 0 1]; % From the first arm coordinates to the body one.
R2 = [cos(q8) -sin(q8) 0; 0 0 1; -sin(q8) -cos(q8) 0]; % From the second arm coordinates to the first one.
R3 = [cos(q9) -sin(q9) 0; sin(q9) cos(q9) 0; 0 0 1]; % From the thrid arm coodinate to the second one.

% The coorinates of each joint relative to the previous joint in the previous joint coordinates.
rp1_p0_xyz = [0; 0; 0];
rp2_p1_1 = [0; 0; 0];
rp3_p2_2 = [L3; 0; 0];

% The angular velocity of the each arm of the manipulator in its own coordinates.
w1_1 = R1.' * ws_xyz + [0; 0; dq7];
w2_2 = R2.' * w1_1 + [0; 0; dq8];
w3_3 = R3.' * w2_2 + [0; 0; dq9];

% The location and the velocity of the base of the manipulator (reference point).
rp0_G_xyz = [0; 0; L1+L2];
Vp0_xyz = R * VG_XYZ + cross(ws_xyz,rp0_G_xyz);

% The velocity of each joint in its coordinates.
Vp1_1 = R1.' * (Vp0_xyz + cross(ws_xyz,rp1_p0_xyz));
Vp2_2 = R2.' * (Vp1_1 + cross(w1_1,rp2_p1_1));
Vp3_3 = R3.' * (Vp2_2 + cross(w2_2,rp3_p2_2));

% The inertia matrix of each arms of the manipulator (It is assumed that the arms are all cylinder and have the same radius).
IG1_1 = [1/12*m1*(3*Radius^2+L2^2) 0 0; 0 1/12*m1*(3*Radius^2+L2^2) 0; 0 0 1/2*m1*Radius^2];
IG2_2 = [1/2*m2*Radius^2 0 0; 0 1/12*m2*(3*Radius^2+L3^2) 0; 0 0 1/12*m2*(3*Radius^2+L3^2)];
IG3_3 = [1/12*m3*(3*Radius^2+L4^2) 0 0; 0 1/2*m3*Radius^2 0; 0 0 1/12*m3*(3*Radius^2+L4^2)];

% The location of the center of mass of each arm relative to its joint.
rG1_p1_1 = [0; 0; -L2/2];
rG2_p2_2 = [L3/2; 0; 0];
rG3_p3_3 = [0; L4/2; 0];

% The veloctiy of the center of mass of each arm relative to its joint.
VG1_1 = Vp1_1 + cross(w1_1,rG1_p1_1);
VG2_2 = Vp2_2 + cross(w2_2,rG2_p2_2);
VG3_3 = Vp3_3 + cross(w3_3,rG3_p3_3);

% The velocity of the center of mass of arms in the inertia coordinates.
VG1_XYZ = R.'*R1*VG1_1;
VG2_XYZ = R.'*R1*R2*VG2_2;
VG3_XYZ = R.'*R1*R2*R3*VG3_3;

% Angular momentum of each arm.
HG1_1 = IG1_1 * w1_1; 
HG2_2 = IG2_2 * w2_2; 
HG3_3 = IG3_3 * w3_3; 

% The kinetic energy of each arm.
T1 = 1/2*m1*(VG1_1.'*VG1_1) + 1/2*w1_1.'*HG1_1;
T2 = 1/2*m2*(VG2_2.'*VG2_2) + 1/2*w2_2.'*HG2_2;
T3 = 1/2*m3*(VG3_3.'*VG3_3) + 1/2*w3_3.'*HG3_3;

%% Evaluating the acceleration of the center of mass of the satellite and the each arm of the manipulator in the inertia coordinates

s = [q, dq];
ds = [dq, ddq];

aG_XYZ = [0; 0; 0];
for i=1:18
aG_XYZ = aG_XYZ + diff(VG_XYZ , s(i)) * ds(i); % Evaluating the time derivation of the components of the veloctiy vector
end
aG_xyz = R * aG_XYZ; % Converting the components from the inetria coordinates to the body one.

aG1_1 = [0; 0; 0];
for i=1:18
aG1_1 = aG1_1 + diff(VG1_1 , s(i)) * ds(i); % Evaluating the time derivation of the components of the veloctiy vector
end
aG1_1 = aG1_1 + cross(w1_1,VG1_1); % Adding the time derivation of the unit vectors of the veloctiy vector
aG1_xyz = R1* aG1_1; % Converting the components from the first arm coordinates to the body one.


% The comments for the following commands are similiar to the previous comments (The arms are just replaced).
aG2_2 = [0; 0; 0];
for i=1:18
aG2_2 = aG2_2 + diff(VG2_2 , s(i)) * ds(i);
end
aG2_2 = aG2_2 + cross(w2_2,VG2_2);
aG2_xyz = R1*R2* aG2_2;

aG3_3 = [0; 0; 0];
for i=1:18
aG3_3 = aG3_3 + diff(VG3_3 , s(i)) * ds(i);
end
aG3_3 = aG3_3 + cross(w3_3,VG3_3);
aG3_xyz = R1*R2*R3* aG3_3;

%% The first three-equation of motion derived from three force equations for the whole system (Both the satellite and the manipulator)

Eq1_3 = ms * aG_xyz + m1 * aG1_xyz + m2 * aG2_xyz + m3 * aG3_xyz;

%% Evaluating the rate of angular momentum of the satellite and the arms of the manipulator in the body coordinates

% Converting the components of the angular momentum from each member coordinates to the body coordinates
HG1_xyz = R1*HG1_1;
HG2_xyz = R1*R2*HG2_2;
HG3_xyz = R1*R2*R3*HG3_3;

% The location of the center of mass of each arm of the manipulator relative to the center of mass of the satellite in the body coordinates
rG1_G_xyz = R1*rG1_p1_1 + (rp1_p0_xyz + rp0_G_xyz);
rG2_G_xyz = R1*R2*rG2_p2_2 + R1*rp2_p1_1 + (rp1_p0_xyz + rp0_G_xyz);
rG3_G_xyz = R1*R2*R3*rG3_p3_3 + R1*R2*rp3_p2_2 + R1*rp2_p1_1 + (rp1_p0_xyz + rp0_G_xyz);

% The time derivation of the angular momentum of the satellite.
dHs_xyz = [0; 0; 0];
for i=1:18
dHs_xyz = dHs_xyz + diff(Hs_xyz , s(i)) * ds(i); % Evaluating the time derivation of the components of the angular momentum vector
end
dHs_xyz = dHs_xyz + cross(ws_xyz,Hs_xyz); % Adding the time derivation of the unit vectors of the angular momentum vector

% The comments of the following commands are similar to the previous comments.
dHG1_xyz = [0; 0; 0];
for i=1:18
dHG1_xyz = dHG1_xyz + diff(HG1_xyz , s(i)) * ds(i);
end
dHG1_xyz = dHG1_xyz + cross(ws_xyz,HG1_xyz);

dHG2_xyz = [0; 0; 0];
for i=1:18
dHG2_xyz = dHG2_xyz + diff(HG2_xyz , s(i)) * ds(i);
end
dHG2_xyz = dHG2_xyz + cross(ws_xyz,HG2_xyz);

dHG3_xyz = [0; 0; 0];
for i=1:18
dHG3_xyz = dHG3_xyz + diff(HG3_xyz , s(i)) * ds(i);
end
dHG3_xyz = dHG3_xyz + cross(ws_xyz,HG3_xyz);

%% The second three-equation of motion through momentum equations for the whole system (Both the satellite and the manipulator)

Eq4_6 = dHs_xyz + dHG1_xyz + cross(rG1_G_xyz, m1*aG1_xyz) + dHG2_xyz + cross(rG2_G_xyz, m2*aG2_xyz) + dHG3_xyz + cross(rG3_G_xyz, m3*aG3_xyz);

%% The equation number 9 through momentum for the third link of the manipulator (In oder to be similar to the lagrange equations, this equation is named "the quation number 9")

dHG3_3 = R3.'*R2.'*R1.'*dHG3_xyz;
aG3_3 = R3.'*R2.'*R1.'*aG3_xyz;
Eq9 = (dHG3_3 + cross(rG3_p3_3, m3*aG3_3)).'*[0; 0; 1] - tau3;

%% The equation number 8 through momentum for the second and thrid link of the manipulator

dHG2_2 = R2.'*R1.'*dHG2_xyz;
aG2_2 = R2.'*R1.'*aG2_xyz;
Eq8 = (R3*dHG3_3 + cross(R3*rG3_p3_3 + rp3_p2_2, R3*m3*aG3_3) + dHG2_2 + cross(rG2_p2_2, m2*aG2_2)).'*[0; 0; 1] - tau2;

%% The equation number 7 through momentum for the first, second and thrid link of the manipulator

dHG1_1 = R1.'*dHG1_xyz;
aG1_1 = R1.'*aG1_xyz;
Eq7 = (R2*R3*dHG3_3 + cross(R2*R3*rG3_p3_3 + R2*rp3_p2_2 + rp2_p1_1, R2*R3*m3*aG3_3) + R2*dHG2_2 + cross(R2*rG2_p2_2 + rp2_p1_1, R2*m2*aG2_2) + dHG1_1 + cross(rG1_p1_1, m1*aG1_1)).'*[0; 0; 1] - tau1;

%% Evaulating the M and B matrices

Eq = [Eq1_3; Eq4_6; Eq7; Eq8; Eq9];

M = simplify(jacobian(Eq, ddq)); % The mass matrix
B = simplify(Eq - M*ddq); % The B matrix (not included the ddqs)

%% Evaluating the total energy, the linear momentum, and the angular momentum

V = 0; % The gravitational potential energy is ignored.
T = Ts + T1 + T2 + T3;
E = simplify(T + V); % The total energy of system.

%%
P = simplify(ms*VG_XYZ + m1*VG1_XYZ + m2*VG2_XYZ + m3*VG3_XYZ); % The total linear momentum of the system

rG_O = [q1; q2; q3];
HSat_O = R.'*Hs_xyz + cross(rG_O,ms*VG_XYZ);

rG1_O = R.'*R1*rG1_p1_1 + R.'*(rp1_p0_xyz + rp0_G_xyz) + rG_O;
HG1_O = R.'*R1*HG1_1 + cross(rG1_O,m1*VG1_XYZ);

rG2_O = R.'*R1*R2*rG2_p2_2 + R.'*R1*rp2_p1_1 + R.'*(rp1_p0_xyz + rp0_G_xyz) + rG_O;
HG2_O = R.'*R1*R2*HG2_2 + cross(rG2_O,m2*VG2_XYZ);

rG3_O = R.'*R1*R2*R3*rG3_p3_3 + R.'*R1*R2*rp3_p2_2 + R.'*R1*rp2_p1_1 + R.'*(rp1_p0_xyz + rp0_G_xyz) + rG_O;
HG3_O = R.'*R1*R2*R3*HG3_3 + cross(rG3_O,m3*VG3_XYZ);

HTot_O = HSat_O + HG1_O + HG2_O + HG3_O;
