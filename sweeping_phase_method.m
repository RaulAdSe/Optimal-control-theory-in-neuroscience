function [phase_diff] = sweeping_phase_method(num_pulses, dt)

% simulate unperturbed neuron
[~, unperturbed] = neuron_simple(0, 0, dt, 500, 0, 5);

% detect unperturbed spikes
unperturbed_spikes = detect_spikes(unperturbed, dt);

% set times of perturbations and preallocate phase difference array.
% Perturbation takes place at the beggining!
ts_per = linspace(unperturbed_spikes(2),unperturbed_spikes(3),num_pulses);
% ts_per = [unperturbed_spikes(2), ts_per, unperturbed_spikes(3)];    % Include extrems to ensure continuity

phase_diff = zeros(1,num_pulses);

% loop over perturbation phases
for j = 1:num_pulses
    
    % simulate perturbed neuron
    [~, perturbed] = neuron_simple(ts_per(j), -0.3, dt, 500, 0, 5);
    
    % detect perturbed spikes
    perturbed_spikes = detect_spikes(perturbed, dt);
    
    % compute phase difference: in the infinite-time limit
    % andinfinitesimally shifted from base points on gamma
    phase_diff(j) = (perturbed_spikes(end) - unperturbed_spikes(end));    %Phase difference is measured at the end
end

% The spike time changes are scaled to obtain the corresponding phase changes. This is done by dividing the time difference between the perturbed and unperturbed spikes by the period of the limit cycle.
phase_diff = phase_diff./8.3995*2*pi;
% Scale PRC to the magnitude of perturvation
phase_diff = phase_diff./0.3;

theta = linspace(0,2*pi,num_pulses);

figure()
plot(theta,phase_diff,'-o');
xlim([0 2*pi])
xlabel('$\theta$','Interpreter','latex');
ylabel('$\frac{\partial \theta}{\partial v}$','Interpreter','latex','Rotation',0);
title('PRC');
end