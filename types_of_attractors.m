% Generates panels A to C of Figure 1 of the paper:
% "Attractor-based models for sequences and pattern generation in neural circuits"
%
% sA1: fixed point attractors
% sA2: sequential attractors
% sA3: fusion attractors
%
% Written by Juliana Londono, July 4, 2024

close all
clear all
addpath('functions')

% --- Graph 1: Static attractor (fixed point) ---
sA{1} = [0 1 0 0;
         1 0 0 1;
         1 1 0 1;
         0 0 1 0];
X0{1} = 0.1 * [rand(1,2), 0, 0];  % Initial condition
proj{1} = [...                   % 2D projection basis
    0.4484, 0.7720;
    0.3658, 0.9329;
    0.7635, 0.9727;
    0.6279, 0.1920];

% --- Graph 2: Sequential attractor ---
sA{2} = [0 1 0 0 0 1;
         1 0 0 0 0 1;
         1 1 0 0 0 0;
         0 0 1 0 1 0;
         0 0 1 1 0 0;
         0 0 0 1 1 0];
X0{2} = 0.1*[rand(1,2),0,0,0,0];
proj{2} = [0.1389, 0.4849;
    0.6963, 0.3935;
    0.0938, 0.6714;
    0.5254, 0.7413;
    0.5303, 0.5201;
    0.8611, 0.3477];

% --- Graph 3: Fusion attractor ---
sA{3} = [0 1 1 1 1;
         1 0 1 1 1;
         1 1 0 0 1;
         1 1 1 0 0;
         1 1 0 1 0];
X0{3} = 0.1*[rand(1,2) 0 0 0] + 0.1*rand(1,5);
proj{3} = [0.2428, 0.3947;
    0.4424, 0.6834;
    0.6878, 0.7040;
    0.3592, 0.4423;
    0.7363, 0.0196];

T = 25; % total simulation time
e = 0.25; % epsilon parameter for CTLNs
d = 0.5; % delta parameter for CTLNs
theta = 1; %theta parameter for CTLNs (external input)

for i=1:size(sA,2)
    figure
    n = size(sA{i},1);
    colors = lines(n);

    % solve ODE
    soln = sA2soln(sA{i},T,X0{i},e,d,theta);

    % plot graph G
    subplot(2,3,1)
    plot_graph(sA{i},colors);
    title(['graph ' int2str(i)]);

    % plot projection
    subplot(2,3,4)
    plot_projection(soln.X,proj{i});
    hold on;
    plot_projection(soln.X,proj{i},[.5,1],[.7 0 0]);
    hold off;
    title('projection (red is last half)')

    % rate curves
    subplot(2,3,[2:3])
    plot_ratecurves(soln.X,soln.time,colors);
    xlabel('time'); 
    hold off;
    subplot(2,3,[5:6])

    % grayscale
    plot_grayscale(soln.X);
    xlabel('time'); 
    title(['X0 = [' num2str(round(X0{i},4)) ']']);
     clear colors
end