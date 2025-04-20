function T = transformation(theta)
T11 = (cos(theta))^2;
T12 = (sin(theta))^2;
T13 = 2 * sin(theta) * cos(theta);
T21 = (sin(theta))^2;
T22 = (cos(theta))^2;
T23 = -2 * sin(theta) * cos(theta);
T31 = - sin(theta) * cos(theta);
T32 = sin(theta) * cos(theta);
T33 = (cos(theta))^2 - (sin(theta))^2;
T = [T11, T12, T13; T21, T22, T23; T31, T32, T33];
end