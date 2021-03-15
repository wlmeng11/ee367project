% MeasureH.m
%
% William Meng
% March 13, 2021
% EE 367 Final Project

clearvars; clc; close all;

% Initialize Field
field_init(-1)

% Transducer Parameters (temporal)
fc = 5e6; % Carrier Frequency [Hz]
c = 1540; % Speed of Sound [m/s]
lambda = c/fc; % Wavelength [m]
fBW = 0.66; % Fractional Bandwidth
fs = 160e6; % Sampling Frequency [Hz]

% Transducer Parameters (spatial)
D = 50*lambda; % Aperture Size
kerf = 25e-6; % Kerf [m]
width = lambda; % Width of each element in x-direction [m]
height = 10e-3; % Width of each element in y-direction [m]
tx_elem = ceil(D/width);  % Number of Elements
rx_elem = tx_elem;  % Number of Elements
elevfoc = 20e-3; % Radius of elevation focus
subx = 1; % Number of subdivisions for x elements
suby = 5; % Number of subdivisions for y elements
focus = [0 0 10e-3]; % Focal Point [m]

% Voxels
x = -D/2:width:D/2-width;
y = zeros(size(x));
zmin = 10*lambda;
zstep = lambda;
zmax = 50*lambda;
z = repmat(zmin:zstep:zmax, length(x), 1);
N = numel(z); % how many pixels in image

fprintf('Wavelength λ = %g mm\n', lambda*1e3);
fprintf('Aperture size D = %g mm\n', D*1e3);
fprintf('# elements = %g\n', tx_elem);
fprintf('N = %g\n', N);

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
rng('default');
seed = 0; % seed for RNG
s = rng(seed);

% Delay mask as varying-thickness plastic layer
c_plastic = 2750; % [m/s]
lambda_plastic = c_plastic/fc; % [m]
sigma = lambda_plastic/2; % [m]
plastic_mask = sigma .* randn(tx_elem, 1); % Gaussian distribution
offset = -1.5 * min(plastic_mask);
plastic_mask = plastic_mask + offset;

figure(1);
subplot(1, 2, 1);
bar(x*1000, plastic_mask * 1e3);
axis square tight;
xlabel('x (mm)');
ylabel('Mask Thickness (mm)');
title('Plastic Mask (Physical)');

% Delay mask as temporal delays for each transducer element
delay_mask = plastic_mask ./ c_plastic; % Gaussian distribution

subplot(1, 2, 2);
bar(1:tx_elem, delay_mask * 1e9);
axis square tight;
xlabel('Element #');
ylabel('Relative Delay (ns)');
title('Delay Mask (Simulation)');

sgtitle('Delay Mask')
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 800 400]);
saveas(gcf, 'Mask.png');

%% Set the Tx and Rx delays
xdc_focus_times(Tx, 0, transpose(delay_mask));
xdc_focus_times(Rx, 0, transpose(delay_mask));
% Question: should Rx delays be reversed in time?
Tx_timeline = xdc_get(Tx, 'focus');
Rx_timeline = xdc_get(Rx, 'focus');

figure(2);
subplot(1, 2, 1);
plot(Tx_timeline);
title('Tx time line');
xlabel('Element #');
ylabel('Delay (s)');
axis square tight;

subplot(1, 2, 2);
plot(Rx_timeline);
title('Rx time line');
xlabel('Element #');
ylabel('Delay (s)');
axis square tight;

sgtitle('Tx and Rx timelines');
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 800 400]);
saveas(gcf, 'TxRxTimelines.png');

%% Plot transmitted field and pulse-echo field
hp = simulate_and_plot_Tx(Tx, x, y, z);
hhp = simulate_and_plot_pulse_echo(Tx, Rx, x, y, z);
K = size(hhp, 1);
t_array = 1/fs * (0:K-1);

%%
figure(3);
clf;
hold on;
for n = 1:500:N
    plot(t_array, hhp(:, n), 'DisplayName', sprintf('n=%d', n));
end
hold off;
xlabel('Time (s)');
title('Pulse-echo responses at a few example pixels');
axis square tight;
legend();
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 500 500]);
saveas(gcf, 'PulseEchoWaveforms.png');

%% Generate the scene to be imaged
 

%% Function definitions
% Based on examples from RAD 235 class

function hp = simulate_and_plot_Tx(Tx, x, y, z)
    % Hydrophone Points
    pos = [repmat(x',size(z,2),1) repmat(y',size(z,2),1) reshape(z,[],1)];

    % Calculate Tx field    
    [hp,start_hp] = calc_hp(Tx,pos);
    
    % Show image of Actual Tx beam
    hp2 = sum(abs(hilbert(hp)),1);
    % Normalize Tx field to see -6dB Tx beam
    txfield = reshape(hp2,length(x),size(z,2))';
    for xx = 1:size(txfield,1)
      txfield(xx,:) = txfield(xx,:)./max(txfield(xx,:));
    end
    figure(5)
    clf
    imagesc(x*1e3,z(1,:)*1e3,db(txfield))
    axis image
    title('Normalized Transmitted Beam')
    xlabel('Lateral Position (mm)')
    ylabel('Depth (mm)')
    axis equal tight;
    set(gcf, 'Color', 'w');
end

function hhp = simulate_and_plot_pulse_echo(Tx, Rx, x, y, z)
    % This code is mostly from beam_example.m

    % Hydrophone Points
    pos = [repmat(x',size(z,2),1) repmat(y',size(z,2),1) reshape(z,[],1)];

    % Calculate pulse-echo response (B = Tx*Rx)
    [hhp,start_hpp] = calc_hhp(Tx,Rx,pos);

    % Show image of Actual Beam
    hhp2 = sum(abs(hilbert(hhp)),1);
    % Normalize Tx field to see -6dB Tx beam
    beam = reshape(hhp2,length(x),size(z,2))';
    for xx = 1:size(beam,1)
      beam(xx,:) = beam(xx,:)./max(beam(xx,:));
    end
    figure(6)
    clf
    imagesc(x*1e3,z(1,:)*1e3,db(beam))
    axis image
    title('Normalized Pulse-Echo Beam')
    xlabel('Lateral Position (mm)')
    ylabel('Depth (mm)')
    axis equal tight;
    set(gcf, 'Color', 'w');
end