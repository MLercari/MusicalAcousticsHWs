syms rho c S1 L1 L2 x1 Z0 k;
Zc = (i*rho*c)/S1*(cot(k*L1) + (k*x1)^(-1))^(-1);
Zin = Z0*( (Zc*cos(k*L2) + i*Z0*sin(k*L2)) / (i*Zc*sin(k*L2) + Z0*cos(k*L2)) );
