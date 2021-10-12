% create a slice of data that is a ternary plot slice simulating Kendall's device for James' data
data = readtable('data-HTS-invert.csv');
conc1 = sort(unique(data.c1));
conc2 = sort(unique(data.c2));
conc3 = sort(unique(data.c3));