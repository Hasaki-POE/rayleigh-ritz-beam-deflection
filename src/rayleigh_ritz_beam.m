% Rayleigh-Ritz Beam Deflection Approximation
% 
% This script approximates the transverse displacement field u(x)
% of a fixed-free Euler-Bernoulli beam using the Rayleigh-Ritz method.
%
% A sixth-order polynomial trial function is assumed. Displacement
% constraints are enforced at the fixed end, midspan, and free end. The
% remaining Ritz coefficients are found by minimizing bending strain energy.
syms x l E I a2 a3 a1 a0 A P a4 a5 a6 u1 u2

u = a0 + a1 * x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 + a6*x^6;

du = diff(u,x);
eq1 = subs(u,  x, 0)   == 0;      % u(0) = 0
eq2 = subs(du, x, 0)   == 0;      % u'(0) = 0
eq3 = subs(u,  x, l/2) == u1;     % u(l/2) = u1
eq4 = subs(u,  x, l)   == u2;     % u(l) = u2

%eliminating coefficients
sol_bc = solve([eq1,eq2,eq3,eq4],[a0,a1,a2,a3])
%sol_bc tells that a0 and a1 is 0
%substitute
u_bc = subs(u,[a0,a1,a2,a3],[sol_bc.a0, sol_bc.a1, sol_bc.a2, sol_bc.a3]);

U = (1/2)*int(E*I*(diff(u_bc,x,2))^2, x, 0, l);
U = simplify(U);


W = 0;
Pi = U-W;

dp4 = diff(Pi,a4);
dp5 = diff(Pi,a5);
dp6 = diff(Pi,a6);



sol = solve([dp4, dp5, dp6],[a4, a5, a6])
a2_final = simplify(subs(sol_bc.a2, [a4, a5, a6], [sol.a4, sol.a5, sol.a6]))

a3_final = simplify(subs(sol_bc.a3, [a4, a5, a6], [sol.a4, sol.a5, sol.a6]))
u_final = simplify(subs(u_bc,[a4,a5,a6], [sol.a4, sol.a5, sol.a6]))
