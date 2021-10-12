clear; clc; close all

j = 1;
results = James2v3_2(j);
resultsOld = results;
resultsTotal = results;
for j = 2:13
    results = James2v3_2(j);
    resultsTotal = vertcat(resultsTotal, results);
end
