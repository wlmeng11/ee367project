% SingleRotation.m
%
% William Meng
% March 13, 2021
% EE 367 Final Project
%
% Perform 2D compressed sensing ultrasound imaging with a single
% transducer, using a single "rotation" (ie. only one coded aperture).

clearvars; clc; close all;

% Initialize Field
field_init(-1)

% Transducer Parameters (temporal)
fc = 5e6; % Carrier Frequency [Hz]
c = 1540; % Speed of Sound [m/s]
lambda = c/fc; % Wavelength [m]
fBW = 0.66; % Fractional Bandwidth
fs = fc*4; % Sampling Frequency [Hz]

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
zmax = 20*lambda;
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

%% Generate coded aperture
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

% Set the Tx and Rx delays
xdc_focus_times(Tx, 0, transpose(delay_mask));
xdc_focus_times(Rx, 0, transpose(delay_mask));
Tx_timeline = xdc_get(Tx, 'focus');
Rx_timeline = xdc_get(Rx, 'focus');

%% Measure pulse-echo response
% (This code is adapted from the RAD 235 homework)

% Sweep a hydrophone through each pixel in the image region
hydrophone_positions = [repmat(x',size(z,2),1) repmat(y',size(z,2),1) reshape(z,[],1)];
% Measure transmitted field
[hp, start_hp] = calc_hp(Tx, hydrophone_positions);
% Measure pulse-echo field
[hhp, start_hpp] = calc_hhp(Tx, Rx, hydrophone_positions);

% transmitted energy at each pixel
hp2 = sum(abs(hilbert(hp)),1);
txfield = reshape(hp2,length(x),size(z,2))';
for xx = 1:size(txfield,1) % Normalize Tx field to see -6dB Tx beam
  txfield(xx,:) = txfield(xx,:)./max(txfield(xx,:));
end

% pulse-echo energy at each pixel
hhp2 = sum(abs(hilbert(hhp)),1); 
beam = reshape(hhp2,length(x),size(z,2))';
for xx = 1:size(beam,1) % Normalize Tx field to see -6dB Tx beam
  beam(xx,:) = beam(xx,:)./max(beam(xx,:));
end

% Plot spatial distribution of energy Tx and pulse-echo fields
figure(2);
subplot(1, 2, 1);
imagesc(x*1e3,z(1,:)*1e3,db(txfield))
axis image
c = colorbar;
c.Label.String = 'Normalized Energy (dB)';
title('Transmitted Field')
xlabel('Lateral Position (mm)')
ylabel('Depth (mm)')

subplot(1, 2, 2);
imagesc(x*1e3,z(1,:)*1e3,db(beam))
axis image
c = colorbar;
c.Label.String = 'Normalized Energy (dB)';
title('Pulse-Echo Field')
xlabel('Lateral Position (mm)')
ylabel('Depth (mm)')

set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 800 400]);
saveas(gcf, 'EnergyFields.png');

%% Image Formation Matrix
H = hhp;
K = size(H, 1); % how many time samples in pulse-echo data
R = 1; % number of rotations
M = K*R;
assert(N == size(H, 2));
compression = N/M;
t_array = 1/fs * (0:K-1);

figure(3);
clf;
hold on;
for n = 1:N/10:N
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


%% Generate a scene to be imaged
% Coordinates of point targets
xmin = x(1);
xmax = x(end);
targets.x = [xmin/2, 0, xmax/2];
targets.y = zeros(size(targets.x));
zrange = zmax - zmin;
targets.z = [zmin + 0.3*zrange, zmin + 0.5*zrange, zmin + 0.8*zrange];
% Format used by calc_scat:
% points should have 3 columns (x,y,z) and a row for each scatterer
% amplitudes should be a column vector with one entry for each scatterer
points = transpose([targets.x; targets.y; targets.z]);
amplitudes = transpose(ones(size(targets.x)));

% Generate v in matrix and vector forms
Nx = tx_elem;
Nz = N/Nx;
scene = zeros(Nz, Nx); % matrix form
% populate with scatterers
for i = 1:length(targets.x)
    index.x = ceil((targets.x(i) - xmin)/width);
    index.z = ceil((targets.z(i) - zmin)/zstep);
    %fprintf('indices: %g, %g\n', index.x, index.z);
    scene(index.z, index.x) = amplitudes(i);
end
v = scene(:); % convert matrix to column vector

figure(4);
subplot(1, 2, 1);
imagesc(scene);
axis equal tight;
colormap gray;
xlabel('x location');
ylabel('z location');
title('v (matrix form)');

subplot(1, 2, 2);
stem(v);
axis square tight;
title('v (vector form)');

sgtitle('True image')
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 800 400]);
saveas(gcf, 'TrueImage.png');

%% Image formation model
Hv = H * v;

% Add Gaussian noise
rng(s);
electronic_SNR = 100;
noise_sigma = max(Hv)/electronic_SNR;
n = noise_sigma * randn(size(Hv));
u = Hv + n;

figure(6);
subplot(1, 3, 1);
plot(t_array, Hv, 'DisplayName', 'Hv');
axis square tight;
xlabel('t (s)');
title('Measurements without noise');
yl = ylim;
legend();

subplot(1, 3, 2);
plot(t_array, n, 'DisplayName', 'n');
axis square tight;
xlabel('t (s)');
title(sprintf('Additive Gaussian noise\nElectronic SNR = %g', electronic_SNR));
ylim(yl); % set same ylim as first plot
legend();

subplot(1, 3, 3);
plot(t_array, u, 'DisplayName', 'u = Hv + n');
axis square tight;
xlabel('t (s)');
title('Measurements with noise');
ylim(yl); % set same ylim as first plot
legend();

set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 1200 400]);
saveas(gcf, 'Hv+n.png');

%% Reconstruction (single rotation)
% redefine them cuz they get overwritten by ADMM code
v = scene(:);
u = Hv + n;

A = H;
At = transpose(H);
Afun = @(x) A*x;
Atfun = @(x) At*x;
AAtfun  = @(x) reshape(Afun(Atfun( x )), [M 1]);
b = u; % measurement
I = scene; % true image

% Least Norm solution with Predetermined Conjugate Gradient (PCG)
maxItersCG = 50;
x = pcg(AAtfun, b(:), 1e-12, maxItersCG);
x_pcg = Atfun(x);
x_pcg2D = reshape(x_pcg, size(scene)); % reshape into 2D image

MSE_lnorm  = mean( (x_pcg - v).^2 );
PSNR_lnorm = 10*log10(1/MSE_lnorm);
fprintf('\nLeast Norm (PCG):\nMSE = %g\nPSNR = %g dB\n', MSE_lnorm, PSNR_lnorm);

% Least Norm solution with Moore-Penrose Pseudo-inverse
AAtinv = pinv(A*At);
x_pinv = Atfun(AAtinv * b);
x_pinv2D = reshape(x_pinv, size(scene));

MSE_pinv  = mean( (x_pinv - v).^2 );
PSNR_pinv = 10*log10(1/MSE_pinv);
fprintf('\nLeast Norm (Pseudo-inverse):\nMSE = %g\nPSNR = %g dB\n', MSE_pinv, PSNR_pinv);

% Plot true image and reconstructed images
figure(7);
subplot(1, 3, 1);
imagesc(scene);
axis equal tight;
colormap gray;
colorbar;
title('Ground Truth');

subplot(1, 3, 2);
imagesc(x_pcg2D);
axis equal tight;
colormap gray;
colorbar;
title(sprintf('Least Norm (PCG) Solution\nPSNR = %g dB', PSNR_lnorm));

subplot(1, 3, 3);
imagesc(x_pinv2D);
axis equal tight;
colormap gray;
colorbar;
title(sprintf('Least Norm (Pseudo-inverse) Solution\nPSNR = %g dB', PSNR_pinv));

sgtitle(sprintf('Compressed Sensing Reconstruction\n# Rotations = %g\nCompression = %g\nElectronic SNR = %g', R, compression, electronic_SNR));
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 1200 400]);
saveas(gcf, 'Reconstructed.png');


%% ADMM with TV regularization
% This code is adapted from the EE 367 homework, but I redefined Afun and
% Atfun to work with my image formation model.

imageResolution = size(scene);

Afun = @(x) A*x(:);
Atfun = @(x) reshape(At*x, imageResolution);
AtAfun  = @(x) Atfun(Afun(x));
opDtDx  = @(x) opDtx(opDx(x));

numItersADMM    = 25;  % number of ADMM iterations
rho             = 1;
lambda          = 0.01;

x = zeros(imageResolution);
z = zeros([imageResolution(1) imageResolution(2) 2]);
u = zeros([imageResolution(1) imageResolution(2) 2]);

Hfun = @(x) reshape(AtAfun(reshape(x,imageResolution)) + rho.*opDtDx(reshape(x,imageResolution)), [prod(imageResolution) 1]); 

PSNR        = zeros([numItersADMM 1]);
residuals   = zeros([numItersADMM 1]);

figure(8);
for k=1:numItersADMM

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x update 
    
    v = z-u;

    bb = Atfun(b) + rho*opDtx(v);

    maxItersCG = 25;
    x = pcg(Hfun,bb(:),1e-12,maxItersCG,[],[],x(:));    
    x = reshape(x, imageResolution);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % z update - soft shrinkage
    
    kappa = lambda/rho;
    v = (opDx(x)+u); 
    z = max(1 - kappa/abs(v),0) .* v;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u update    
    u = u+opDx(x)-z;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PSNR & residual
    
    r1 = b-Afun(x);
    r2 = opDx(x); 
    residuals(k) = 0.5*sum(r1(:).^2) + lambda.*sum( abs(r2(:)) );  
    
    MSE     = mean( (x(:)-I(:)).^2 );
    PSNR(k) = 10*log10(1/MSE);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot
    
    subplot(2,3,[1 4]);
    imshow(I);    
    title('Target Image');
    
    subplot(2,3,[2 5]);
    imshow(x);    
    title(['Compression=' num2str(compression) ', PSNR=' num2str(PSNR(k),'%3.2f') 'dB, \lambda=' num2str(lambda), ' \rho=' num2str(rho) ', noise \sigma=' num2str(sigma)]);    

    subplot(2,3,3);
    plot(1:numItersADMM, PSNR, 'LineWidth', 2, 'color', [1 0 1]);
    title('PSNR');
    xlabel('Iteration');
    ylabel('PSNR in dB');
    grid on;
    ylim([10 45]);
    
    subplot(2,3,6);
    plot(1:numItersADMM,log(residuals), 'LineWidth', 2);
    title('log(residual)');
    xlabel('Iteration');
    ylabel('Value');
    grid on;
    ylim([6 9]);
    
    drawnow;
end

sgtitle('ADMM TV');
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 1200 600]);
saveas(gcf, 'ADMM_TV.png');
