function Qbar = Qbar(P,theta)
V21 = ( P(2) * P(3) ) / P(1);
Q11 = P(1) / (1 - P(3) * V21);
Q12 = ( ( P(3) * P(2) ) / ( 1 - P(3) * V21));
Q22 = P(2) / ( 1 - P(3) * V21);
Q66 = P(4);
Q16 = 0;
Q26 = 0;
Q_11 = Q11 * (cos(theta)) ^ 4  + 2 * (Q12 + 2 * Q66) * ( (sin(theta))^2 ) * ( (cos(theta))^2 ) + Q22 *  (sin(theta))^4 ;
Q_12 = (Q11 + Q22 - 4 * Q66) *( (cos(theta) ) ^ 2 * (sin(theta)) ^ 2 ) + Q12 * ( (sin(theta))^4 + (cos(theta))^4 );
Q_22 = Q11 *  ( sin(theta) ) ^ 4  +  2 * (Q12 + 2*Q66) * ( (cos(theta) )^2 * (sin(theta))^2 ) + Q22*( cos(theta) ) ^ 4;
Q_16 = ( Q11 - Q12 - 2 * Q66) * ( sin(theta) ) * ( cos(theta) ) ^ 3 - (Q22 - Q12 - 2 * Q66) * ( sin(theta) ) ^3 * ( cos(theta) );
Q_26 = ( Q11 - Q12 - 2 * Q66)* ( sin(theta) ) ^3 * ( cos(theta) );
Q_66 = ( Q11 + Q22 - 2 * Q12 - 2 * Q66) * ( (cos(theta))^2 * (sin(theta))^2 )+ Q66 * ( (sin(theta))^4 + ( cos(theta) ) ^ 4 );
Qbar = [Q_11,Q_12,Q_16;Q_12,Q_22,Q_26;Q_16,Q_26,Q_66] * 10 ^ 9;
end