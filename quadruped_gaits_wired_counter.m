% Builds and simulates a 3-layer CTLN network for quadruped gait sequences.
% Layer 1 (L1): counter network driving the sequence.
% Layer 2 (L2): cycle of auxiliary nodes controlling access to Layer 3.
% Layer 3 (L3): gait network with 5 gait attractors (bound, pace, trot, walk, pronk).
%
% Generates panels D and E of Figure 5 in the manuscript:
% "Attractor-based models for sequences and pattern generation in neural circuits"
%
% Panel D: activity in all 3 layers under theta pulses to L1 and L2, grayscale plot,
% pulse plot, and rate curves for limb (L3), L2, and L1 neurons.
%
% Panel E: full connectivity matrix (W) for the 3-layer network.
%
% Written by Juliana Londono

% 13 bound
% 15 pace
% 17 trot
% 19 walk
% 23 pronk

close all
clear all
addpath('functions')


%sequence must be a sequence of auxiliary nodes to stimulate!
sequence = [15,13,23,17,19,23,17]; %this is the one from the paper!
aux_nodes = sort(unique(sequence));
indep_sz = size(aux_nodes,2);
m = size(sequence,2); %number of pulses
e = .25; d = .5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define your layers
%first layer is counter
sAcell{1,1} = [
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0;
]; %counter
%second layer is independent set
sAcell{2,2} = zeros(indep_sz,indep_sz);
%third layer is gaits
sAcell{3,3} = [0 1 0 1 1 0 1 0 1 0 0 0 1 0 0 1 1 0 0 0 0 1 0 1;
    1 0 1 0 0 1 0 1 0 1 0 0 0 1 0 1 0 1 1 0 0 0 0 1;
    0 1 0 1 1 0 1 0 0 0 1 0 0 1 1 0 1 0 0 0 1 0 0 1;
    1 0 1 0 0 1 0 1 0 0 0 1 1 0 1 0 0 1 0 1 0 0 0 1;
    1 0 1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 0 0 0 1 0 0;
    0 1 0 1 1 0 1 0 0 1 0 0 0 1 0 1 0 1 1 0 0 0 0 0;
    1 0 1 0 0 1 0 1 0 0 1 0 0 1 1 0 1 0 0 0 1 0 0 0;
    0 1 0 1 1 0 1 0 0 0 0 1 1 0 1 0 0 1 0 1 0 0 0 0;
    1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
    0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
    0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];

%number of layers
N = size(sAcell,2);

% convert to W and get size of layers
for i = 1:N
    Wcell{i,i} = graph2net(sAcell{i,i},e,d);
    n_array(i) = size(Wcell{i,i},1);
end

%network size
n = sum(n_array);

%connections from indep. set to gaits:
temp = zeros(n_array(3),n_array(2)); %zeros(size of receiver,size of sender)
cols = 1:n_array(2); %sender: indep. set
rows = aux_nodes; %receiver: aux nodes of gaits
for i = 1:size(cols,2)
    temp(rows(i),cols(i)) = 1;
end
sAcell{3,2} = temp;

%translate gait sequence to L2 indices:
for i=1:size(sequence,2)
    idx = find(aux_nodes==sequence(i));
    translated_sequence(i) = idx;
end

%connection from counter to indep set:
% temp = zeros(n_array(2),n_array(1)); %zeros(size of receiver,size of sender)
% cols = [1:2:2*m,2:2:2*m+1]; %sender: counter
% rows = [translated_sequence,translated_sequence]; %receiver: indep.set
% for i = 1:size(cols,2)
%     temp(rows(i),cols(i)) = 1;
% end
sAcell{2,1} =  [0,0,0,0,1,1,0,0,0,0,0,0,0,0;
                1,1,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,1,1,0,0,0,0;
                0,0,1,1,0,0,0,0,0,0,1,1,0,0;
                0,0,0,0,0,0,1,1,0,0,0,0,0,0];

% order of gaits:
% pace walk bound pronk trot walk
sequence_names = {'rand','pace','walk','bound','pronk','trot','walk'};
sequence_names = string(strjoin(sequence_names));

% 1 bound 15
% 2 pace 16
% 3 trot 17
% 4 walk 18
% 5 pronk 19

%build the W with feed forward layered structure
for i = 1:N
    for j = 1:N
        if j+1==i %feedforward connections are as specified above
            A = sAcell{i,j};
            Wcell{i,j} = -ones(size(A))+A*e + (A-ones(size(A)))*d;
        elseif i>j | j>i  %everyone else is 0
            Wcell{i,j} = zeros(size(Wcell{i,i},1),size(Wcell{j,j},1));
        end
    end
end

%build the big W from the blocks
for i = 1:size(Wcell,1)
    for j = 1:size(Wcell,1)
        %record the partition indices taus
        if i==j
            tau{i} = (sum(n_array(1:i-1))+1):sum(n_array(1:i));
        end
        W((sum(n_array(1:i-1))+1):sum(n_array(1:i)),(sum(n_array(1:j-1))+1):sum(n_array(1:j))) = Wcell{i,j};
    end
end

figure(1)
imagesc(W)
hold on
for j=1:size(n_array,2)-1
    xline(sum(n_array(1:j))+0.5,'m')
    yline(sum(n_array(1:j))+0.5,'m')
end
hold off
colormap(gray)
colorbar('Ticks',unique(W),...
    'TickLabels',{'-1-\delta','-1+\epsilon','0'})
xticks(1:n)
yticks(1:n)
title(['epsilon = ',num2str(e),', delta = ', num2str(d),...
    ', sequence = ', num2str(sequence)])

w = 2; %width of pulses

T = [25 w 50 w 35 w 30 w 50 w 40 w 25 w 31];

theta = ones(n,1);
theta(tau{2}) = 0; %silence layer 2 external inputs
b = repmat(theta,1,length(T));

h = 12; %pulse height
for j=1:7 
    b(tau{1},2*j) = h/2;
    b(tau{2},2*j) = h;
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

% intial condition coincides with clique, no gait
X0 = zeros(1,n); %all nodes off...
X0([1,2]) = 0.6;
% X0([20:43]) = 0.2*rand(size(X0([20:43])));
X0([20:43]) = [0.0452, 0.0341, ...
        0.0455, 0.0871, 0.0622, ...
        0.1847, 0.0860, 0.0370, ...
        0.1810, 0.1959, 0.0878, ...
        0.0222, 0.0516, 0.0817, ...
        0.1190, 0.0524, 0.1206, ...
        0.1422, 0.0443, 0.0235, ...
        0.0593, 0.0638, 0.0848, ...
        0.1016];


%solve ODE
soln = threshlin_ode(W,b,T,X0);

gait_colors = [202,140,43;...%1
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

counter_nodes = 1:n_array(1);
indep_nodes = (n_array(1)+1):(n_array(1)+n_array(2));
leg_nodes = (n_array(1)+n_array(2))+1:(n_array(1)+n_array(2))+4;
%gait_specific_aux_nodes = aux_nodes + n_array(1)+n_array(2);
gait_specific_aux_nodes = (n_array(1)+n_array(2))+13:(n_array(1)+n_array(2))+24;


%color of cycle nodes matches the color of CPG neuron
cycle_nodes_colors = gait_colors(aux_nodes,:);

counter_nodes_colors = [166,206,227;...
31,120,180;...
178,223,138;...
51,160,44;...
251,154,153;...
227,26,28;...
253,191,111;...
255,127,0;...
202,178,214;...
106,61,154;...
210,180,140;...
177,89,40;...
201,148,199;...
221,28,119]/255;

colors = [counter_nodes_colors;cycle_nodes_colors;gait_colors];

% plot solution and pulses
figure(3)
subplot(11,1,1:4)
%plot_grayscale(soln.X(:,n_array(1)+n_array(2)+4:-1:1));
plot_grayscale(soln.X(:,n:-1:1));
%plot_grayscale(soln.X);
hold on
for j=1:size(n_array,2)-1
    %these divide the layers:
    yline(-sum(n_array(1:j))+n+0.5,'m')
    %these divide the gaits:
    %yline(-sum(n_array(1:j))+n+0.5,'m')
end
hold off
title(['Sequence: ', sequence_names])

subplot(11,1,5)
plot(t,pulse,'-k');
title(['theta pulses with width = ',num2str(w)])
xlim([0,max(t)])
ylim([0.5,h+0.5]);
yticks(unique(pulse))

subplot(11,1,6:7)
plot_ratecurves(soln.X(:,leg_nodes),soln.time,colors(leg_nodes,:));
xlim([0,max(t)])
title(['leg nodes (',num2str(leg_nodes),') of L3'])
legend(strsplit(num2str(leg_nodes)))

subplot(11,1,8:9)
plot_ratecurves(soln.X(:,indep_nodes),soln.time,colors(indep_nodes,:));
xlim([0,max(t)])
title('L2 nodes')
legend(strsplit(num2str(indep_nodes)))

subplot(11,1,10:11)
plot_ratecurves(soln.X(:,[counter_nodes]),...
    soln.time,colors([counter_nodes],:));
xlim([0,max(t)])
title('L1 nodes')

legend(strsplit(num2str([counter_nodes])))
