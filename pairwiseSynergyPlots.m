% Pairwise model-predicted vs actual 3-drug combinations
% T-ALL cell lines
%  3: CCRF-CEM
% 10: DND-41
%  9: Jurkat
%  6: MOLT-4

% B-ALL cell lines
%  4: 697
%  2: KOPN-8
%  1: NALM-6
% 11: RCH-ACV
%  7: REH
%  8: RES4;11
%  5: UOC-B1

% ETP-ALL cell lines
% 13: LOUCY
% 12: PEER

clear;clc;close all


lineage = 'T+ETP-ALL';

switch lineage
case 'T-ALL'
    cellLineIdx = [3, 10, 9, 6];
case 'B-ALL'
    cellLineIdx = [4, 2, 1, 11, 7, 8, 5];
case 'ETP-ALL'
    cellLineIdx = [13, 12];
case 'T+ETP-ALL'
    cellLineIdx = [3, 10, 9, 6, 13, 12];
otherwise
    disp('Invalid lineage name')
end


% to override the lineage switch case, uncomment the line below
% and type in a number between 1-13 for specific cell line
%cellLineIdx = 13;

% calculate predicted viability and synergy scores
for i = 1:numel(cellLineIdx)
    results{i} = pair_vs_3drug(cellLineIdx(i));
    close all
end

% vertically concatenate tables
data = results{1};
for i = 2:numel(cellLineIdx)
    data = [data; results{i}];
end


% plotting
% figure 1: only the pairwise interactions
figure
hold on
plot(data.g123_actual,data.g123_wood,'b.', 'MarkerSize', 10)
%plot(data.g123_actual, data.g123_bliss, 'r.', 'MarkerSize', 10)
hold off

xlabel('Actual Viability')
ylabel('Predicted Viability')

axis([0 1.1 0 1.1])

rline = refline([1,0]);
rline.Color = [0.1, 0.1, 0.1];

legend({'Pairwise Interaction Model', 'Reference Line'}, 'Location', 'Southeast')


% figure 2: pairwise and total interaction
figure
hold on
plot(data.g123_actual,data.g123_wood,'b.', 'MarkerSize', 10)
plot(data.g123_actual, data.g123_bliss, 'r.', 'MarkerSize', 10)
hold off

xlabel('Actual Viability')
ylabel('Predicted Viability')

axis([0 1.1 0 1.1])

rline = refline([1,0]);
rline.Color = [0.1, 0.1, 0.1];

legend({'Pairwise Interaction Model', 'Total Interaction', 'Reference Line'}, 'Location', 'Southeast')