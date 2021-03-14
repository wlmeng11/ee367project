% MeasureH.m
%
% William Meng
% March 13, 2021
% EE 367 Final Project

clearvars; clc; close all;

% Initialize Field
field_init(-1)

% Transducer Parameters (temporal)
fc = 3e6; % Carrier Frequency [Hz]
c = 1540; % Speed of Sound [m/s]
lambda = c/fc; % Wavelength [m]
fBW = 0.66; % Fractional Bandwidth
fs = 160e6; % Sampling Frequency [Hz]

% Transducer Parameters (spatial)
D = 30e-3; % Aperture Size
kerf = 25e-6; % Kerf [m]
width = lambda; % Width of each element in x-direction [m]
height = 10e-3; % Width of each element in y-direction [m]
tx_elem = ceil(D/width);  % Number of Elements
rx_elem = tx_elem;  % Number of Elements
elevfoc = 20e-3; % Radius of elevation focus
subx = 1; % Number of subdivisions for x elements
suby = 5; % Number of subdivisions for y elements
focus = [0 0 30e-3]; % Focal Point [m]

fprintf('Wavelength λ = %g mm\n', lambda*1e3);
fprintf('Aperture size D = %g mm\n', D*1e3);
fprintf('# elements = %g\n', tx_elem);

% Set Field Parameters
set_field('fs',fs);

% Define Impulse Response and Excitation (to get transmit response)
tc = gauspuls('cutoff', fc, fBW, -6, -40);
t = -tc:1/fs:tc;
impulseResponse = gauspuls(t,fc,fBW);
excitationPulse = 1;

% Sample correction for length of impulse response
tshift = (size(conv(impulseResponse,conv(impulseResponse, ...
    excitationPulse)),2))/fs;


% Tx Aperture
Tx = xdc_focused_array(tx_elem,width,height,kerf,elevfoc,subx,suby,focus);
xdc_impulse(Tx, impulseResponse);
xdc_excitation(Tx, excitationPulse);
xdc_center_focus(Tx, [0 0 0]);
xdc_focus(Tx, 0, focus);

% Rx Aperture
Rx = xdc_focused_array(rx_elem,width,height,kerf,elevfoc,subx,suby,focus);
xdc_impulse(Rx,impulseResponse);
xdc_center_focus(Rx, [0 0 0]);
xdc_focus(Rx,0,focus);

%% Generate pseudorandom delay mask
rng('default'); s = rng; % set seed=0 for RNG

% Delay mask
c_plastic = 2750; % [m/s]
sigma = 1/2 * lambda/c_plastic; % [s]
delay_mask = sigma .* randn(tx_elem, 1); % Gaussian distribution

figure(1);
subplot(1, 2, 1);
bar(delay_mask);
axis square;
xlabel('Element #');
ylabel('Relative Delay (s)');
title('Delay Mask (Simulation)');

% Corresponding delay mask as varying-thickness plastic layer
plastic_mask = c_plastic .* delay_mask;
offset = -1.5 * min(plastic_mask);
plastic_mask = plastic_mask + offset;

subplot(1, 2, 2);
bar(plastic_mask * 1e3);
axis square;
xlabel('x (a.u.)');
ylabel('Mask Thickness (mm)');
title('Plastic Mask (Physical)');

set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 800 400]);
saveas(gcf, 'Mask.png');


% Set the Tx delays

% Set the Rx delays
% TO DO: should it be reversed in time?


%% Plot transmitted field


%% Compute pulse-echo response by autoconvolution of the transmitted field


