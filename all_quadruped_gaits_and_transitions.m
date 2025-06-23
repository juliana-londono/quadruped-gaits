% Simulates a sequence of gait transitions in the 5-gait network using 
% targeted theta pulses, and also simulates the network starting from
% different initial conditions.
% Generates panels B and D of Figure 4 in the manuscript:
% "Attractor-based models for sequences and pattern generation in 
% neural circuits"
%
% Panel B: activity from different initial conditions showing 5 gaits.
% Panel D: pulse-driven transitions between the gaits.
%
% Written by Juliana Londono

close all
clear all
addpath('functions')

gait_idx = [13,15,17,19,23]; %neuron that needs to be stimulated for each gait
sA = [
    0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1;
    1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1;
    0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1;
    1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1;
    1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
    0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0;
    1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0;
    0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0;
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
];
n = size(sA,1);

gait_colors = [153,103,50;...%1
    203,33,40;...%2
    18,128,178;...%3
    42,153,71;...%4
    1,1,1;...%5
    1,1,1;...%6
    1,1,1;...%7
    1,1,1;...%8
    1,1,1;...%9
    1,1,1;...%10
    1,1,1;...%11
    1,1,1;...%12
    241,96,106;...%13
    241,96,106;...%14
    127,197,237;...%15
    127,197,237;...%16
    105,193,134;...%17
    105,193,134;...%18
    255,207,75;...%19
    255,207,75;...%20
    255,207,75;...%21
    255,207,75;...%22
    239,87,37;...%23
    239,87,37]/255;%24

e = 0.25;
d = 0.5;

%% plot different initial conditions

T = 50;
theta = 1;

for i = 1:numel(gait_idx)

    %solve ODE
    X0 = zeros(n,1);
    X0(gait_idx(i)) = 0.3;
    soln = sA2soln(sA,T,X0,e,d,theta);

    %plot solution
    figure

    %greyscale
    subplot(2,1,1)
    plot_grayscale(soln.X);
    title(['\epsilon = ', num2str(e), ', \delta = ', num2str(d)])

    %ratecurves
    subplot(2,1,2)
    plot_ratecurves(soln.X(:,1:4),soln.time,gait_colors(1:4,:));
    xlim([0,sum(T)]);

end

%% transitions
pulsed_neurons = {17,19,15,13,23,17};
w = 2; h = 2; b = 1;
no_pulses = size(pulsed_neurons,2);
T = [88, 2, 98, 2, 57, 2, 60, 2, 76, 2, 74, 2, 83];

%pulse times
for i = 2:2:2*no_pulses
    T(i) = w;
end
   b = b*ones(n,length(T));

for j = 1:no_pulses %for each pulse:
    sigma = pulsed_neurons{j};
    for i = 1:length(sigma) %for each neuron in the j-th group
        b(sigma(i),2*j) = h;
    end
end

%this t,pulse is useful to plot the pulses plot(t,pulse,'-k');
t(1) = 0;
for i=1:size(b,2)
    t(2*i) = sum(T(1:i));
    t(2*i+1) = sum(T(1:i));
    pulse(:,2*i-1) = b(:,i); %doubled theta pulse? yes, to plot that little rectangle
    pulse(:,2*i) = b(:,i);
end
t = t(1:size(pulse,2)); % cut off last value

%solve ODE
X0 = zeros(n,1);
X0(19) = 0.3;
soln = sA2soln(sA,T,X0,e,d,b);

%plot solution
figure

%greyscale
subplot(2,1,1)
plot_grayscale(soln.X);
title(['\epsilon = ', num2str(e), ', \delta = ', num2str(d)])

%pulse
subplot(7,1,4)
colororder(gait_colors)
plot(t,pulse);
ylim([min(min(pulse)-0.2),max(max(pulse)+0.2)]);
xlim([0,sum(T)]);
yticks(unique(b))
ylabel('\theta')

%ratecurves
subplot(2,1,2)
plot_ratecurves(soln.X(:,1:4),soln.time,gait_colors(1:4,:));
xlim([0,sum(T)]);

