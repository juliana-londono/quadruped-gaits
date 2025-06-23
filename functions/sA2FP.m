function FPstruct = sA2FP(sA,e,d,theta,print_flag)

% function FPstruct = sA2FP(sA,e,d,theta,print_flag)
%
% FPstruct = structure with the following fields:
% FPstruct.FP = {sigma1, sigma2,...} list of fixed point supports
% FPstruct.coreFP = subset of FP that are core motifs
% FPstruct.index = [1 -1 -1 1 ...] vector of with the index for each element of FP
% FPstruct.fixpts = mxn matrix, each row is a fixed point, ordered as in FP
% FPstruct.corefixpts = lxn matrix, fixed points corresponding to coreFP
% FPstruct.core_idx = vector with indices of the core fixpts
% FPstruct.stability = vector of booleans checking stability of fixed pt for elements of FP
%
% finds FP(G) = FP(G,e,d), and the associated fixed points, for a CTLN
% with graph G give by adj mtx sA
%
% this is a wrapper for "check_fixedpt"; also calls 
% graph2net.m, check_core_motif.m, and check_stability.m
%
% created july 19, 2018, as a fancier version of sA2fixpts.m
% updated may 30, 2020 to use "core motif" language instead of "min permitted"
% updated on June 11, 2020 to check fixed points grouped by the
%   size of their supports
% updated on June 16, 2020 to print out fixed point supports
% last modified on June 17, 2020 to include a print_flag of whether or not
%   to print out the fixed point supports
% last modified on June 21, 2020 to include core_idx and stability in the FPstruct 



if nargin<2 || isempty(e)
    e = []; % use default of graph2net
end;

if nargin<3 || isempty(d)
    d = []; % use default of graph2net
end;

if nargin<4 || isempty(theta)
    theta = 1;
end;

if nargin<5 || isempty(print_flag)
    print_flag=1;
end;

%...................................
n = size(sA,1);
W = graph2net(sA,e,d);
b = theta*ones(n,1);

FP = {}; % cell array for fixed point supports
fixpts = []; % matrix for fixed points
stability=[];
index=[];
TF_core_motif = []; % binary vector with 1 for each min permitted elt of FP
j = 0;

for k=1:n
    sets=nchoosek(1:n,k); % create all the subsets of 1..n of size k
    for i=1:size(sets,1)
        sig = sets(i,:);
        [TF x_fp] =  check_fixedpt(W,sig,b);
        if TF
            j = j+1;
            FP(j) = {sig};
            fixpts = [fixpts; x_fp'];
            index(j) = sign(det(eye(length(sig))-W(sig,sig)));
            stability(j)=check_stability(W,sig);
            TF_core_motif(j) = check_core_motif(sA(sig,sig),e,d);
            if print_flag
                fprintf(['sgn = ' int2str(index(j)) ', core = ' int2str(TF_core_motif(j)) ', supp = ' int2str(sig) '\n']);
            end
        end;
    end;
end;


core_idx = find(TF_core_motif); % get indices for core motifs

if print_flag
    fprintf(['num of fixed points = ' int2str(size(fixpts,1)) '\n']);
    fprintf(['num of core motifs = ' int2str(length(core_idx)) '\n']);
end

% package everything into FPstruct
FPstruct.FP = FP;
FPstruct.coreFP = FP(core_idx);
FPstruct.index = index;
FPstruct.fixpts = fixpts;
FPstruct.corefixpts = fixpts(core_idx,:);
FPstruct.core_idx = core_idx;
FPstruct.stability=stability;