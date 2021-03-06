function results = pair_vs_3drug(j)
tiledlayout(1,2)
cellLineCol = 5:17;

% import datafile
data = importData; % imports file james_CL_invert.csv
conc1 = data.c1;
conc2 = data.c2;
conc3 = data.c3;

idx3drug = find(conc1~=0 & conc2~=0 & conc3~=0);

% This data is arranged such that no concentrations are 0
% looking at all possible 3 drug combination sets


%cellLineNum = 13; % NEED TO DO THIS ACROSS ALL CELL LINES
cellLineNum = j;
E = data{:,cellLineCol(cellLineNum)};
%E = data{:, cellLineCol(j)};
for i = 1:length(idx3drug) % for every 3 drug combination
    g1idx(i,1) = find(conc1==conc1(idx3drug(i)) & conc2==0 & conc3==0);
    g2idx(i,1) = find(conc1==0 & conc2==conc2(idx3drug(i)) & conc3==0);
    g3idx(i,1) = find(conc1==0 & conc2==0 & conc3==conc3(idx3drug(i)));
    g12idx(i,1) = find(conc1==conc1(idx3drug(i)) & conc2==conc2(idx3drug(i)) & conc3==0);
    g13idx(i,1) = find(conc1==conc1(idx3drug(i)) & conc2==0 & conc3==conc3(idx3drug(i)));
    g23idx(i,1) = find(conc1==0 & conc2==conc2(idx3drug(i)) & conc3==conc3(idx3drug(i)));
end

% since this is from inverted dataset, don't need to subtract
g1 = E(g1idx); % subtract from 1 to get effect (survival?)
g2 = E(g2idx);
g3 = E(g3idx);
g12 = E(g12idx);
g13 = E(g13idx);
g23 = E(g23idx);
g123_actual = E(idx3drug);

% Application of models
g123_wood = wood(g1,g2,g3,g12,g13,g23);
g123_bliss = bliss(g1,g2,g3);

E3 = zeros(height(g1), 1); % initialize variable
DA = E3;
for i = 1:numel(E3)
    [DA(i), E3(i)] = RBI3(g1(i), g2(i), g3(i), g12(i), g13(i), g23(i), g123_actual(i));
end



% for each pair, also display synergy score and categorization
DA12 = zeros(height(g1), 1);
DA13 = DA12;
DA23 = DA12;

for i = 1:height(g1)
    [DA12(i), syn12{i}] = RBI(g1(i), g2(i), g12(i));
    [DA13(i), syn13{i}] = RBI(g1(i), g3(i), g13(i));
    [DA23(i), syn23{i}] = RBI(g2(i), g3(i), g23(i));
end
syn12 = syn12';
syn13 = syn13';
syn23 = syn23';

for i = 1:height(g1)
    D_crit = 0.5; % threshold for criticality
    if abs(E3(i)) < D_crit
        syn123{i} = 'additivity';
    elseif E3(i) < 0
        syn123{i} = 'synergy';
    elseif E3(i) > 0
        syn123{i} = 'antagonistic';
    else
        syn123{i} = 'Unknown';
    end
end
syn123 = syn123';

% Plotting
nexttile
plot(g123_actual,g123_wood,'b.', 'MarkerSize', 10)
hold on
plot(g123_actual, g123_bliss, 'r.', 'MarkerSize', 10)
xlabel('Experimental')
ylabel('Model Predicted')
axis([0 1.1 0 1.1])
refline([1,0])
hold off
legend({'Pairwise Model', 'Net Interaction'}, 'Location', 'Southeast')

nexttile
plot(g123_actual, g123_bliss, 'r.')
xlabel('Experimental')
ylabel('Model Predicted')
axis([0 1.1 0 1.1])
refline(1,0)


% R2 results
mdl_emergent = fitlm(g123_actual, g123_bliss);
mdl_pairwise = fitlm(g123_actual, g123_wood);

N = numel(g123_actual); % number of samples
K_emergent = 3; % number of parameters in emergent model, g1, g2, g3
K_pairwise = 6; % number of parameters in pairwise model, + g12, g23, g13

R2_emergent = mdl_emergent.Rsquared.Ordinary;
R2_emergent_adj = 1 - ((1-R2_emergent)*(N-1)/(N-K_emergent-1));

R2_pairwise = mdl_pairwise.Rsquared.Ordinary;
R2_pairwise_adj = 1 - ((1-R2_pairwise)*(N-1)/(N-K_pairwise-1));


% custom R2 calculations
residuals = g123_actual - g123_wood;

meany = mean(g123_actual);
error2 = (g123_actual - g123_bliss).^2;
distmean2 = (g123_actual - meany).^2;

sumsqerr = sum(error2);
sumsqvar = sum(distmean2);
R2_emer = 1 - sumsqerr/sumsqvar;
R2_ermer_adj = 1 - ((1-R2_emer)*(N-1)/(N-3-1));

error3 = (g123_actual - g123_wood).^2;
distmean3 = (g123_actual - meany).^2;
sumsqerr = sum(error3);
sumsqvar = sum(distmean3);
R2_pair = 1 - sumsqerr/sumsqvar;
R2_pair_adj = 1 - ((1-R2_pair)*(N-1)/(N-6-1));



% results table
results = table(repmat(cellLineNum, length(idx3drug), 1), idx3drug, ...
    conc1(idx3drug), conc2(idx3drug), conc3(idx3drug), ...
    g1, g2, g3, g12, g13, g23, g123_actual, g123_bliss, g123_wood, ...
    DA, E3, DA12, DA13, DA23, syn12, syn13, syn23, syn123, residuals);

results = renamevars(results,["Var1","Var3","Var4","Var5"], ...
["cellLine","conc1","conc2","conc3"]);

% CUSTOM FUNCTIONS ============================================================

% Wood et al paper
function g123 = wood(g1, g2, g3, g12, g13, g23)
g123 = g1.*g23 + g2.*g13 + g3.*g12 - 2.*g1.*g2.*g3;
end

% Zimmer et al paper
function g123 = zimmer(g1, g2, g3, g12, g13, g23)
g123 = g12.*g13.*g23./(g1.*g2.*g3);
end

function g123 = bliss(g1, g2, g3)
g123 = g1.*g2.*g3;
end

% Beppler et al paper
function [DA, E3, DAscaled, E3scaled] = beppler(g1, g2, g3, g12, g13, g23, g123)
DA = g123 - g1.*g2.*g3; % difference from Bliss
E3 = g123 - wood(g1, g2, g3, g12, g13, g23);

% rescaling
% when DA or E3 are negative, replace g122 with 0
for i = 1:length(DA)
    if DA(i) < 0 % synergistic
        DAscaled(i,1) = DA(i)./(abs(0-bliss(g1(i),g2(i),g3(i))));
    elseif DA(i) > 0
        DAscaled(i,1) = DA(i)./(abs(min([g1(i), g2(i), g3(i)])-bliss(g1(i),g2(i),g3(i))));
    elseif DA(i) == 0
        DAscaled(i,1) = 0;
    else DAscaled(i,1) = 0;
    end

    if E3(i) < 0
        E3scaled(i,1) = E3(i)./abs(0 - wood(g1(i),g2(i),g3(i),g12(i),g13(i),g23(i)));
    elseif E3(i) > 0
        E3scaled(i,1) = E3(i)./abs(min([g1(i), g2(i), g3(i)]) - wood(g1(i),g2(i),g3(i),g12(i),g13(i),g23(i)));
    elseif E3(i) == 0
        E3scaled(i,1) = 0;
    else E3scaled(i,1) = 0;
    end
end

end


end % end of function