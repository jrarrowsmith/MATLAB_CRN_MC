function [a_range, a, a_histogram] = normal_sample_pdf(a_mean, a_std, num_std, n, truncate)
%This function computes the normally distributed pdf
%It assumes that the parameter is a normal distribution centered on a_mean
%with num_std standard deviations (elev_std) and samples it n times.
%JRA October 20, 2016
%Updated October 23, 2018 fot the truncation
%updated and generalized Oct. 14, 2020
%variables that are returned:
% a_range = range of sample values
% a = samples
% a_histogram = histogram output over the a_range
%truncate = 1 yes cut it at the 2 sigma or no if = other


[increment, a_range] = compute_increment_range(a_mean+num_std.*a_std, a_mean-num_std.*a_std, n);
a = a_mean+a_std.*randn(length(a_range),1);
%here we are going to truncate the gaussian at 2*std
if truncate==1
locs_minus = find(a<a_mean+num_std.*a_std);
a=a(locs_minus);
locs = find(a>a_mean-num_std.*a_std);
a=a(locs);
end
%but that cuts out some that we actually need to have the same number as we
%had initially input

nn=length(a); %Length of the concatenated list of rates
numsamples=n; %Sample it the same number of times as the first sampling
m=ceil(rand(numsamples,1).*nn); %Choose randomly and evenly across the length of the concatenated list of rates
sampled_a=[];
for i = 1:numsamples
    sampled_a(i)=a(m(i)); %Sample the composite pdf; THIS IS WHAT WE USE GOING FORWARD
end
a=sampled_a';
a_histogram = hist(a, length(a_range));
end