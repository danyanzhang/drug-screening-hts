% go up to previous directory
function data = importData()
addpath(['..', filesep, 'DataDerived'])
data = readtable('james_CL_invert.csv');
rmpath(['..', filesep, 'DataDerived']);
end