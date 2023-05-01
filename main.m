close all
clear all 
clc

alpha = 1;   % Regulates the energy constraint. Defaul is set to 1.
             % Preguntar Jose tot el formalisme que vam fer a classe
             % Alpha no juega ninbun papel

% Define main parameters
dt = 10^-3;
tol = 10^-8;
it_max = 14;
T_original =  8.3995;    % igual que la del article 8.3995; 
w = 2*pi/T_original;
T1 = 1.2 * T_original;
t_end = 100;

%% PRC and simple drawis


% Per a; fer el dibuix Fig. 1A. Treure comentaris dels plots, baxar el temps de simulació.
% Potser canviar la funció neurona simple incloent el temps d'integració
% com a input i també l'opció de treure plot o no per poder mostrar aquestes dues amb 15 ms i les demès sense
% plot a 1000.
figure()
neuron_simple(0,0,dt,20,1,5);
neuron_simple(6,-5,dt,20,1,5);

% Si vull obtenir la seua corva he de simular per t molt gran. Dins de la
% funció fico t_end 1000.
[phase_diff] = sweeping_phase_method(20,dt);
% Hi ha vegades que no surt el plot d'aqui! PRC

%% Control in phase space

% Interpolo (precisió) i obtinc la derivada de Z. Necessito mateixa
% ressolució de t que de theta.
tspan = [0 T1]; % ms. Simulació llarga per a que poder comparar fases a l'infinit
t = tspan(1):dt:tspan(2);

theta = linspace(0,2*pi,20);
theta_int = linspace(0,2*pi,length(t)+2);  %interpolated domain
Z = spline(theta,phase_diff,theta_int);

dtheta = theta_int(2)-theta_int(1);
Z_dot = diff(Z)/dtheta;
Z_dot_dot = diff(Z_dot)/dtheta;

theta_int = theta_int(1:end-2);
Z = Z(1:end-2);
Z_dot = Z_dot(1:end-1);


% Resolc el sistema acoplat + Implemento Newton per fer shooting.

lambda0 = -15;
fun = tol + 1;  % per inicialitzar
it = 0;

while it <= it_max && abs(fun) > tol
    y = zeros(length(t),2);     % Sistema
    x = zeros(length(t),2);     % Variational system
    y(1,:) = [0 lambda0]; % Initialize the first row with the initial conditions
    x(1,:) = [0 1];
    
    dydt = @(t, theta_ix, lambda) [w + (lambda * Z(theta_ix)^2)/ (2*alpha);
                 -(lambda^2*Z(theta_ix)*Z_dot(theta_ix))/(2*alpha)];

    for i = 1:length(t)-1

% XAPO EL EULER PER ERRORS D'INTEGRACIÓ
%         %theta_ev serpa el theta que utilitze per evaluear les funcions Z
%         %aquest theta ha d'estar el més prop possible a theta real
%         [~, theta_ev_ind] = min(abs(theta_int-mod(theta,2*pi)));
%         %aquest index serà el mateix que per Z
% 
%         y(i+1,:) = [theta lambda] + dt * [w + (lambda * Z(theta_ev_ind)^2)/ 2;
%                                      -(lambda^2*Z(theta_ev_ind)*Z_dot(theta_ev_ind))/2;]';
% 
%         x(i+1,:) = [v1 v2] + dt * [lambda * Z(theta_ev_ind) * Z_dot(theta_ev_ind) * v1 + Z(theta_ev_ind)^2 / 2 * v2;
%                                  -(lambda^2 * v1)/2 * (Z_dot(theta_ev_ind)^2 + Z_dot_dot(theta_ev_ind) * Z_dot(theta_ev_ind)) - lambda * Z(theta_ev_ind) * Z_dot(theta_ev_ind) * v2;]';
%    

% FICO EL RK4

        theta = y(i,1);
        lambda = y(i,2);
        [~, theta_ev_ind] = min(abs(theta_int-mod(theta,2*pi)));
        k1 = dt*dydt(t(i), theta_ev_ind, lambda);
        [~, theta_ev_ind_k1] = min(abs(theta_int-mod(theta + k1(1)/2,2*pi)));
        k2 = dt*dydt(t(i)+dt/2, theta_ev_ind_k1, lambda + k1(2)/2);
        [~, theta_ev_ind_k2] = min(abs(theta_int-mod(theta + k2(1)/2,2*pi)));
        k3 = dt*dydt(t(i)+dt/2, theta_ev_ind_k2, lambda + + k2(2)/2);
        [~, theta_ev_ind_k3] = min(abs(theta_int-mod(theta + k3(1),2*pi)));
        k4 = dt*dydt(t(i)+dt, theta_ev_ind_k3, lambda + k3(2));
        y(i+1,:) = y(i,:) + (1/6)*(k1 + 2*k2 + 2*k3 + k4)';
        
        % COM LES VARIACIONALS SON PER CONVERGIR JA EM VAN BÉ
        v1 = x(i,1);
        v2 = x(i,2);
        x(i+1,:) = [v1 v2] + dt * 1/alpha * [lambda * Z(theta_ev_ind) * Z_dot(theta_ev_ind) * v1 + Z(theta_ev_ind)^2 / 2 * v2;
                         -(lambda^2 * v1)/2 * (Z_dot(theta_ev_ind)^2 + Z_dot_dot(theta_ev_ind) * Z_dot(theta_ev_ind)) - lambda * Z(theta_ev_ind) * Z_dot(theta_ev_ind) * v2;]';
    end

    df = x(int16(T1/dt),1);
    fun = y(int16(T1/dt),1) - 2*pi;
    
    lambda0 = lambda0 - fun/df;    
    it = it + 1;
end

% Trobar Z(theta(t)). Tot ho tinc en funció de theta_ind (fine grid 0,2pi)
theta_t_ind = zeros(1,length(y(:,1)));
for i = 1:length(y(:,1))
    [~, theta_t_ind(i)] = min(abs(theta_int-y(i,1)));
end

theta_t = theta_int(theta_t_ind);
Z_theta_t = Z(theta_t_ind);

u_opt = (y(:,2).*Z_theta_t')/(2*alpha);

% Create a 2x2 grid of subplots
subplot(2,2,1);
plot(t, y(:,1));
xlim([0 T1])
% ylim([0 10])
ylabel('$\theta(t)$','Interpreter','latex');
xlabel('t');

subplot(2,2,2);
plot(theta_int, Z);
xlim([0 theta_int(end)])
% ylim([-0.2 0.2])
ylabel('$Z(\theta(t))$','Interpreter','latex');
xlabel('$\theta$','Interpreter','latex');

subplot(2,2,3);
plot(t, y(:,2));
% ylim([-18 -14])
xlim([0 T1])
ylabel('$\lambda(t)$','Interpreter','latex');
xlabel('t');

subplot(2,2,4);
plot(t, u_opt);
% ylim([-2 0.5])
xlim([0 T1])
ylabel('$u_{opt}(t)$','Interpreter','latex');
xlabel('t');

% LA RAÒ PER LA QUE NO ENGANXA BÉ. 
% Z(theta(t)) és periòdica. He vist que el PRC ho és (extrems) i també imposo que theta sigui periòdic.
% lambda no ho acaba de ser. i no precisament per error d'integració (rk4) 
% POSTPROCESSING DE LA LAMBDA!
nova_y2 = y(1:length(y)-100,2);
dnova_y2 = diff(nova_y2)./diff(t(1:end-100));
x2 = t(end-99:end);
nova_y2_2 = spline([t(end-99) t(end) t(end-100)], [nova_y2(end) lambda0 nova_y2(end-1)],x2);
nova_y2_def = [nova_y2' nova_y2_2];

u_opt_2 = (nova_y2_def'.*Z_theta_t')/2;

% Plotejo fins a 10. A partir de u tendeix a un valor constant
%% Control in original space
t_end = 500;
tspan = [0 t_end]; % ms. Simulació llarga per a que poder comparar fases a l'infinit
t = tspan(1):dt:tspan(2);

% u_opt_T1 = u_opt(1:(T1/dt));   % Temps on he caluclat la u
% Com el sistema és periòdic amb freqüència natural w he de tornar a
% aplicar el mateix per a no tornar a w si vull mantindre'm a T1!

G_opt = trapz(t(1:int16(T1/dt)),u_opt(1:int16(T1/dt)).^2);

% L'EMPALME ENTRE FUNCIONS CONSEQUTIVES ESTÀ MALAMENT FET, HAURIA DE SER
% SMOOTH. SI NO ENGANXA BÉ ES PERQUÈ NO HE FICAT PROU O CORRECTAMENT ELS
% PUNTS PER CALCULAR PRC

% distances = pdist2(u_opt(1), u_opt);
% [sortedDistances, indices] = sort(distances);
% selected = indices(1:4);
% sel_sel = selected>9500;
% selected_ind = selected(sel_sel);
% 
% 
% fun_per_enganxar = u_opt(1:selected_ind);
% n_concat = ceil(t_end/t(selected_ind));   % El denominador hauria de ser T1
% 
% u_opt_conc = repmat(fun_per_enganxar, n_concat, 1);
% u_opt_conc = u_opt_conc(1:length(t));

n_concat = ceil(t_end/T1);
u_opt_conc = repmat(u_opt, n_concat, 1);
u_opt_conc = u_opt_conc(1:length(t));

figure()
plot(t,u_opt_conc)
ylabel('Optimal control');
xlabel('time (s)');
xlim([0 100])



figure()
[~, controlled] = neuron_controlled(u_opt_conc, t_end, dt);
[~, uncontrolled] = neuron_controlled(zeros(1,length(t)), t_end, dt);
legend('controlled trajectory', 'uncontrolled trajectory')
xlim([0 100])

figure()
[~, controlled] = neuron_controlled(u_opt_conc, t_end, dt);
[~, uncontrolled] = neuron_controlled(zeros(1,length(t)), t_end, dt);
legend('controlled trajectory', 'uncontrolled trajectory')
xlim([400 500])


[spikes_controlleed] = detect_spikes(controlled,dt);
[spikes_uncontrolled] = detect_spikes(uncontrolled,dt);

figure()
plot(diff(spikes_controlleed))
yline(T1);
ylabel('T controlled');
xlabel('Spike events');
legend('Controlled neuron inter-spike interval','Desired periodicity','Location','best')

T_controlled = mean(diff(spikes_controlleed(int16(end/2):end)));    %Aquí a la neurona controlada no ha trobat el seu lloc encara. He de donar-li temps. 
T_uncontrolled = mean(diff(spikes_uncontrolled));
disp(['The original period has changed a ' num2str((T_controlled-T_uncontrolled)/T_uncontrolled*100) '% in the controlled system'])
% Els del paper asseguren haver conseguit un 1.185, jo en això tinc 1.195957%
% el que demano que augmenti el periode: T1 = 1.2 * T_original.


figure()
plot3(controlled(:,1),controlled(:,2),controlled(:,3))
hold on
plot3(uncontrolled(:,1),uncontrolled(:,2),uncontrolled(:,3))
hold on
plot3(uncontrolled(1,1), uncontrolled(1,2), uncontrolled(1,3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
hold on
plot3(controlled(1,1)+u_opt(1), controlled(1,2), controlled(1,3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r');
xlim([min(controlled(:,1)) max(controlled(:,1))]);
ylim([min(controlled(:,2)) max(controlled(:,2))]);
zlim([min(controlled(:,3))  max(controlled(:,3))]);
xlabel('V');
zlabel('r');
legend('Controlled', 'Uncontrolled','Location','best');
ylabel('h');
title('Trajectories in original space')
%% Comparar G_opt amb u constant.
% L'efecte d'aquest control (res de closed loops) és equivalent a augmentar
% la intensitat base. Anem a fer un plot T vs Ib i amb això podrem
% identificar el valor de u equivalent.

Ib_vec = linspace(2,12,50);
T_vec = zeros(1,length(Ib_vec));
for i = 1:length(Ib_vec)
    [~,y] = neuron_simple(0, 0, dt, 100, 0, Ib_vec(i));
    spikes = detect_spikes(y, dt);
    T_vec(i) = mean(diff(spikes));
end

figure()
plot(Ib_vec,T_vec,'-o')
yline(T1)
xlabel('Ib');
ylabel('T');

% Es veu que clava prou bé, no em cal interpolar per trobar el valor exacte
[~,index] = min(abs(T1-T_vec));
Ib_T1 = Ib_vec(index);

% El valor d'I_b està fixat a 5, per tant tot el que no sigui per això
% haurà d'anar a u_cte. Factor de C_m multiplicant, amb valor 1

u_cte = (Ib_T1 - 5) * 1;
t_end = 500;
u_cte = u_cte * ones(1,t_end/dt);
[t,y_cont_cte] = neuron_controlled(u_cte, t_end, dt);
[t,y_cont_opt] = neuron_controlled(u_opt_conc, t_end, dt);
title('Controlled models')
legend('Constant control','Optimal control')
xlim([0 100])

[t,y_cont_cte] = neuron_controlled(u_cte, t_end, dt);
[t,y_cont_opt] = neuron_controlled(u_opt_conc, t_end, dt);
title('Controlled models')
legend('Constant control','Optimal control')
xlim([400 500])
% Els dos controls fan el que han de fer, però cadascun a una fase
% diferent.

G_cte = trapz(t(1:int16(T1/dt)),u_cte(1:int16(T1/dt)).^2);

G = [G_cte, G_opt];
disp(['Total control energy consumed:',  newline  ,'Constant control: ' num2str(G(1)),  newline ,' Optimal control: ' num2str(G(2))])

% Comparison between both controls
figure()
plot(t,u_opt_conc)
hold on
plot(t,[u_cte u_cte(1)])
xlim([0 T1])
legend('Optimal control','Constant control')
%% Check de les roots de lambda0

tspan = [0 T1]; % ms. Simulació llarga per a que poder comparar fases a l'infinit
t = tspan(1):dt:tspan(2);

theta = linspace(0,2*pi,20);
theta_int = linspace(0,2*pi,length(t));  %interpolated domain
Z = spline(theta,phase_diff,theta_int);

dtheta = theta_int(2)-theta_int(1);
Z_dot = diff(Z)/dtheta;
Z_dot_dot = diff(Z_dot)/dtheta;

theta_int = theta_int(1:end-2);
Z = Z(1:end-2);
Z_dot = Z_dot(1:end-1);


% Resolc el sistema acoplat + Implemento Newton per fer shooting.

lambda0s = linspace(-25,-10,10);
thetas_end = zeros(1,length(lambda0s));

for kk = 1:length(lambda0s)
      lambda0 = lambda0s(kk);

    dydt = @(t, theta_ix, lambda) [w + (lambda * Z(theta_ix)^2)/ 2;
                     -(lambda^2*Z(theta_ix)*Z_dot(theta_ix))/2];

    % Use the 4th-order Runge-Kutta method to integrate the system
    y = zeros(length(t), 2);
    y(1,:) = [0 lambda0]; % Initialize the first row with the initial conditions
    for i = 1:length(t)-1
        theta = y(i,1);
        lambda = y(i,2);
        [~, theta_ev_ind] = min(abs(theta_int-mod(theta,2*pi)));
        k1 = dt*dydt(t(i), theta_ev_ind, lambda);
        [~, theta_ev_ind_k1] = min(abs(theta_int-mod(theta + k1(1)/2,2*pi)));
        k2 = dt*dydt(t(i)+dt/2, theta_ev_ind_k1, lambda + k1(2)/2);
        [~, theta_ev_ind_k2] = min(abs(theta_int-mod(theta + k2(1)/2,2*pi)));
        k3 = dt*dydt(t(i)+dt/2, theta_ev_ind_k2, lambda + + k2(2)/2);
        [~, theta_ev_ind_k3] = min(abs(theta_int-mod(theta + k3(1),2*pi)));
        k4 = dt*dydt(t(i)+dt, theta_ev_ind_k3, lambda + k3(2));
        y(i+1,:) = y(i,:) + (1/6)*(k1 + 2*k2 + 2*k3 + k4)';
    end
     thetas_end(kk) = y(end,1);
end


figure()
scatter(lambda0s,thetas_end)
hold on
yline(2*pi)
ylabel('$\theta(t_{end})$','Interpreter','latex');
xlabel('$\lambda_0$','Interpreter','latex');

% specify folder path
foldername = 'C:\Users\raulc\Documents\MATLAB\Mamme\Semi\pics';

% save all plots to folder
for i = 1:length(findall(0,'Type','figure'))
    fig = figure(i);
    filename = fullfile(foldername, sprintf('plot_%d.png', i));
    saveas(fig, filename);
end