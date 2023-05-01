function [t,y] = neuron_simple(p_time, p_ampl, dt, t_end, bool,I_b)

% p_time: time which perturbation is applied (in a delta shape)
% p_ampl: amplitude of the perturbation

% Set initial conditions. I ensure to depart really close to the limit
% cycle.
v0 = -50.33067; % mV
h0 = 0.0985757;
r0 = 0.00171423;

% Set simulation parameters
tspan = [0 t_end]; % ms. Simulació llarga per a que poder comparar fases a l'infinit
t = tspan(1):dt:tspan(2);

% Set model parameters
C = 1; % microF/cm^2
g_L = 0.05; % mS/cm^2
E_L = -70; % mV
g_Na = 3; % mS/cm^2
E_Na = 50; % mV
g_K = 5; % mS/cm^2
E_K = -90; % mV
g_T = 5; % mS/cm^2
E_T = 0; % mV
% I_b = 5; %mA/cm^2

% Define the auxiliary functions
h_inf = @(v) 1/(1+exp((v+41)/4));
r_inf = @(v) 1/(1+exp((v+84)/4));
alpha_h = @(v) 0.128*exp(-(v+46)/18);
beta_h = @(v) 4/(1+exp(-(v+23)/5));
tau_h = @(v) 1/(alpha_h(v) + beta_h(v));
tau_r = @(v) (28+exp(-(v+25)/10.5));
m_inf = @(v) 1/(1+exp(-(v+37)/7));
p_inf = @(v) 1/(1+exp(-(v+60)/6.2));

% Build the perturbation in form of delta
pert = zeros(1,length(t));
if p_ampl ~= 0
    pert(int16(p_time/dt)) = p_ampl;
end

% Define the thalamic neuron model
y = zeros(length(t),3);
y(1,:) = [v0 h0 r0]; % Initialize the first row with the initial conditions
for i = 1:length(t)-1
    v = y(i,1);
    h = y(i,2);
    r = y(i,3);
    y(i+1,:) = [v h r] + [pert(i);0;0]' + dt * [(1/C)*(-g_L*(v-E_L) - g_Na*m_inf(v)^3*h*(v-E_Na) - g_K*((0.75*(1-h))^4)*(v-E_K) - g_T*p_inf(v)^2*r*(v-E_T) + I_b);
                                 (h_inf(v) - h)/tau_h(v);
                                 (r_inf(v) - r)/tau_r(v);]';
end

if bool
    % Plot the results
    plot(t,y(:,1));
    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');
    title('Perturbation of Thalamic neuron model');
    hold on  % Este hold on va genial per a fer un parell de crides a la funció i veure la Figura 1A.
end
end