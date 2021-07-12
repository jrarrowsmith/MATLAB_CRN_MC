clc
clear all
close all

%Script to calculate paleoerosion rates from [10Be] - EEZ
%modified by JRA to try Monte Carlo for error propagation
%_______________________________________________________________________
%Set the MC values
n=1000000; %number of samples
num_std=4; %number of standard deviations to sample for the different variables. 
%We assume normal distributions in general
truncate = 1; %truncate = 1 yes cut the distribution at the num_std sigma or no if = other
min_erosion_rate = 0.053; %just to help with plotting [mm/yr]
max_erosion_rate = 0.062; %just to help with plotting [mm/yr]
erosion_rate_step = 0.0001; %bin size for the final histogram [mm/y]
sigmas=0.99; %0.6827 is one sigma and 0.9545 is two sigma range to cut pdf 1.0 takes all 0.99 cleans it
simple_error = 0.05; %simple error assumption for input values in the absence of other information
sample_name=' CHB11';
calc_name=strcat(date,sample_name);

%_______________________________________________________________________
%Set constant variables for erosion rate calculations
%(sd is standard deviation)
atten_spal = 165;                               %spallation attenuation length [gm/cm^2]
atten_spal_sd = atten_spal.*simple_error;       %spallation attenuation length sd
atten_fast = 1500;            %fast muons attenuation length [gm/cm^2]
atten_fast_sd = atten_fast.*simple_error;       %fast muons attenuation length sd
atten_slow = 5300;            % slow muons atten length [gm/cm^2]
atten_slow_sd = atten_slow.*simple_error;       %slow muons atten length sd
lambda = 4.998e-7;           %10Be decay constant (yr^-1)
lambda_sd = lambda.*simple_error;       %10Be decay constant sd
rho_hill = 2.65;             %density of hillslope material [g/cm^3]
rho_hill_sd = rho_hill.*simple_error;       %density of hillslope material sd
rho_burial = 1.9;            %density of burial overburden [g/cm^3]
rho_burial_sd = rho_burial.*simple_error;       %density of burial overburden

%and sample for normal distribution for each 
[atten_spal_range, atten_spal_samples, atten_spal_histogram] = normal_sample_pdf(atten_spal, atten_spal_sd , num_std, n, truncate);
[atten_fast_range, atten_fast_samples, atten_fast_histogram] = normal_sample_pdf(atten_fast, atten_fast_sd , num_std, n, truncate);
[atten_slow_range, atten_slow_samples, atten_slow_histogram] = normal_sample_pdf(atten_slow, atten_slow_sd , num_std, n, truncate);
[lambda_range, lambda_samples, lambda_histogram] = normal_sample_pdf(lambda, lambda_sd , num_std, n, truncate);
[rho_hill_range, rho_hill_samples, rho_hill_histogram] = normal_sample_pdf(rho_hill, rho_hill_sd , num_std, n, truncate);
[rho_burial_range, rho_burial_samples, rho_burial_histogram] = normal_sample_pdf(rho_burial, rho_burial_sd , num_std, n, truncate);

figure(1)
clf
ncols=2; nrows=3;
hold on
subplot(nrows, ncols, 1)
histogram(atten_spal_samples)
xlabel('atten\_spal [gm/cm^2]')
subplot(nrows, ncols, 2)
histogram(atten_fast_samples)
xlabel('atten\_fast [gm/cm^2]')
subplot(nrows, ncols, 3)
histogram(atten_slow_samples)
xlabel('atten\_fast [gm/cm^2]')
subplot(nrows, ncols, 4)
histogram(lambda_samples)
xlabel('10Be decay constant \lambda (yr^-1)')
subplot(nrows, ncols, 5)
histogram(rho_hill_samples)
xlabel('\rho hill [g/cm^3]')
subplot(nrows, ncols, 6)
histogram(rho_burial_samples)
xlabel('\rho burial [g/cm^3]')


%Set sample-specific variables (Using CHB11 values here)______________________

%Production rate values from time-dependent Lal-Stone ('Lm') - different
%rates for each sample
%Need to ssociate an uncertainty with each of these production rate values
P_hill = 9.20;                    %Total hillslope production rate [atoms/g/yr]
P_hill_sd = P_hill.*simple_error;       %Total hillslope production rate sd
Pspal_burial = 4.578;        %Spallation production at burial site [atoms/g/yr]
Pspal_burial_sd = Pspal_burial.*simple_error;       %Spallation production sd
Pmf_burial =  0.03;               %Fast muons production at burial [atoms/g/yr]
Pmf_burial_sd = Pmf_burial.*simple_error;       %Fast muons production sd
Pms_burial  = 0.06;               %Slow muon production at burial [atoms/g/yr]
Pms_burial_sd = Pms_burial.*simple_error;       %Slow muon production sd

N_mes = 108584;        %Measured [10Be] concentration
N_mes_E = 2770;        %Measured [10Be] uncertainty assume normal dist'n

t = 112082;             %Depositional age of sample [yr]
t_E = 10565;            %1 sigma depositional age uncertainty [yr]

AR = 0.034;             %Sediment accumulation rate above sample [cm/yr]
AR_E = AR.*simple_error;       %1 sigma sed accumulation uncertainty

%and sample for normal distribution for each (production rates)
[P_hill_range, P_hill_samples, P_hill_histogram] = normal_sample_pdf(P_hill, P_hill_sd , num_std, n, truncate);
[Pspal_burial_range, Pspal_burial_samples, Pspal_burial_histogram] = normal_sample_pdf(Pspal_burial, Pspal_burial_sd , num_std, n, truncate);
[Pmf_burial_range, Pmf_burial_samples, Pmf_burial_histogram] = normal_sample_pdf(Pmf_burial, Pmf_burial_sd , num_std, n, truncate);
[Pms_burial_range, Pms_burial_samples, Pms_burial_histogram] = normal_sample_pdf(Pms_burial, Pms_burial_sd , num_std, n, truncate);


figure(2)
clf
ncols=2; nrows=2;
hold on
subplot(nrows, ncols, 1)
histogram(P_hill_samples)
xlabel('P\_hill [atoms/g/yr]')
subplot(nrows, ncols, 2)
histogram(Pspal_burial_samples)
xlabel('Pspal\_burial [atoms/g/yr]')
subplot(nrows, ncols, 3)
histogram(Pmf_burial_samples)
xlabel('Pmf\_burial [atoms/g/yr]')
subplot(nrows, ncols, 4)
histogram(Pms_burial_samples)
xlabel('Pms\_burial [atoms/g/yr]')

%and sample for normal distribution for each (other important parameters)
[N_mes_range, N_mes_samples, N_mes_histogram] = normal_sample_pdf(N_mes, N_mes_E , num_std, n, truncate);
[t_range, t_samples, t_histogram] = normal_sample_pdf(t, t_E , num_std, n, truncate);
[AR_range, AR_samples, AR_histogram] = normal_sample_pdf(AR, AR_E , num_std, n, truncate);

figure(3)
clf
ncols=1; nrows=2;
hold on
subplot(nrows, ncols, 1)
histogram(N_mes_samples)
xlabel('Measured [10Be] concentration')
subplot(nrows, ncols, 2)
histogram(t_samples)
xlabel('Depositional age of sample [yr]')
subplot(nrows, ncols, 2)
histogram(AR_samples)
xlabel('Sediment accumulation rate [cm/yr]')



%Calculate post-depositional 10Be accumulation___________________________

%Set grouped variables (accumulation rate * burial density / attenuation length)
ADA_spal = AR_samples .* rho_burial_samples ./ atten_spal_samples;
ADA_fast = AR_samples .* rho_burial_samples ./ atten_fast_samples;
ADA_slow = AR_samples .* rho_burial_samples ./ atten_slow_samples;
%(decay constant * time [age])
LT = lambda_samples .* t_samples;

%Calculate for spallation
Post_spal = Pspal_burial ./ (lambda - ADA_spal) .* (exp(-t*ADA_spal) - exp(-LT));
%Calculate for fast muons
Post_fast = Pmf_burial ./ (lambda - ADA_fast) .* (exp(-t*ADA_fast) - exp(-LT));
%Calculate for slow muons
Post_slow = Pms_burial ./ (lambda - ADA_slow) .* (exp(-t*ADA_slow) - exp(-LT));

%Total post-depositional accumulation
N_post = Post_spal + Post_fast + Post_slow;


%Calculate nuclide loss to radioactive decay____________________________

N_dec = ((N_mes - N_post) ./ exp(-LT)) - (N_mes - N_post);


%Calculate initial nuclide concentration ______________________________
N_0 = N_mes + N_dec - N_post;

figure(4) %nuclides
clf
ncols=1; nrows=3;
hold on
subplot(nrows, ncols, 1)
histogram(N_post)
xlabel('Post depositional nuclides')
subplot(nrows, ncols, 2)
histogram(N_dec)
xlabel('Decayed nuclides')
subplot(nrows, ncols, 3)
histogram(N_0)
xlabel('Initial nuclides')


%Calculate erosion rate [cm/yr]_____________________________________
E_cm = (atten_spal ./ rho_hill) * ((P_hill ./ N_0) - lambda);

%Convert to [mm/yr]
E_mm = E_cm .* 10;

edges = min_erosion_rate:erosion_rate_step:max_erosion_rate;


figure(5) %erosion rate
clf
%ncols=1; nrows=3;
hold on
%subplot(nrows, ncols, 1)
h=histogram(E_mm, edges);
xlabel('Hillslope erosion rate (mm/yr)')
title(calc_name)

figure(6)
clf
[most_common_rate, max_rate, min_rate] = clean_and_process_pdf(h, sigmas, calc_name, 0.0585)
title(calc_name)
xlabel('Hillslope erosion rate (mm/yr)')


