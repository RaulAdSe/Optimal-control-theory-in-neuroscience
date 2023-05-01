function [spikes] = detect_spikes(signal,dt)
% detects spikes in a signal based on threshold crossing
spike_threshold = -10; % in mV. Ho he copiat del paper.
detection_window = 5; % in ms

spikes = [];

for i = 2:length(signal)-1
    if signal(i) > spike_threshold && signal(i) > signal(i-1) && signal(i) > signal(i+1)
        if isempty(spikes) || (i - spikes(end) > round(detection_window/dt))
            spikes = [spikes i];
        end
    end
end

spikes = spikes*dt; % convert to time in ms
end