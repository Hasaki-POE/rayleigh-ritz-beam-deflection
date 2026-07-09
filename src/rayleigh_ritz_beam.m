% Rayleigh-Ritz Beam Deflection Approximation
% 
% This script approximates the transverse displacement field u(x)
% of a fixed-free Euler-Bernoulli beam using the Rayleigh-Ritz method.
%
% A sixth-order polynomial trial function is assumed. Displacement
% constraints are enforced at the fixed end, midspan, and free end. The
% remaining Ritz coefficients are found by minimizing bending strain energy.
syms x L E I a2 a3 a1 a0 A P a4 a5 a6 u1 u2

% Define sixth-order polynomial function
u = a0 + a1 * x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 + a6*x^6;

% first derivative
du = diff(u,x);
% Constraints
eq1 = subs(u,  x, 0)   == 0;      % u(0) = 0
eq2 = subs(du, x, 0)   == 0;      % u'(0) = 0
eq3 = subs(u,  x, L/2) == u1;     % u(L/2) = u1
eq4 = subs(u,  x, L)   == u2;     % u(L) = u2

% eliminating coefficients
sol_bc = solve([eq1,eq2,eq3,eq4],[a0,a1,a2,a3])
%sol_bc tells that a0 and a1 is 0

% substitute
u_bc = subs(u,[a0,a1,a2,a3],[sol_bc.a0, sol_bc.a1, sol_bc.a2, sol_bc.a3]);

% conpute bending strain energy
U = (1/2)*int(E*I*(diff(u_bc,x,2))^2, x, 0, L);
U = simplify(U);

% No external work due to force is included since u1 and u2 are prescribed displacements
W = 0;
Pi = U-W;

% minimize potential energy to remaining coefficients
dp4 = diff(Pi,a4);
dp5 = diff(Pi,a5);
dp6 = diff(Pi,a6);

sol = solve([dp4, dp5, dp6],[a4, a5, a6])

% final coefficients
a2_final = simplify(subs(sol_bc.a2, [a4, a5, a6], [sol.a4, sol.a5, sol.a6]))

a3_final = simplify(subs(sol_bc.a3, [a4, a5, a6], [sol.a4, sol.a5, sol.a6]))
u_final = simplify(subs(u_bc,[a4,a5,a6], [sol.a4, sol.a5, sol.a6]))

%%
% Example numerical values for plotting
L_val = 1;
E_val = 1;
I_val = 1;
u1_val = -0.02;
u2_val = -0.05;

% Substitute numerical values into final displacement field
u_plot = subs(u_final, [L, E, I, u1, u2], [L_val, E_val, I_val, u1_val, u2_val]);

% Generate x values
x_vals = linspace(0, L_val, 200);

% Evaluate displacement
u_vals = double(subs(u_plot, x, x_vals));

% Plot displacement curve
figure;
plot(x_vals / L_val, u_vals, 'LineWidth', 2);
grid on;
xlabel('Normalized Position, x/L');
ylabel('Transverse Displacement, u(x)');
title('Rayleigh-Ritz Beam Deflection Shape');
