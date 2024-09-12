% Modified by Yuexin 2024-05
% To replicate Havlicek et Uludag, 2020. Figure 5
% Steady-state example demonstrating variable relative CBF across depths 
% with different amounts of baseline CBV increase of asceding veins (AV)
% towards the gray matter (GM) surface. 

close all; clear all;

set(0,'DefaultAxesFontSize', 14, ...
      'defaultLineLineWidth', 2, ...
      'defaultLineMarkerSize',15,...
      'DefaultAxesTitleFontWeight', 'normal');      

K = 6;                       % We will consider six depths
P{1} = LBR_parameters(K);    % Get parameter structure for LBR model
                             % By default we consider 40 sec stimulus
                             % duration in order to reach steady-state (i.e P{1}.T = 40)
% Define laminar profile relative CBF (considering six depth) with 60% at 
% the top and low depths and 30% in the middle depths 
cbf{1}      = kron([1.6,1.6,1.76,1.76,1.6,1.6],ones(P{1}.T/P{1}.dt,1));
P{1}.s_d = 0.3;  % increase of baseline CBV in the AV towards the surface
P{1}.alpha_v = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15]';
%P{1}.x_v = [1.25, 1.25, 1.25, 1.25, 1.25, 1.25]';
P{1}.E0v = 0.35;
P{1}.E0d = 0.35;
P{1}.E0p = 0.35; 
% Call the LBR model to generate the laminar BOLD profile:
[LBR{1},Y{1}] = LBR_model(P{1},cbf{1});

P{2}     = P{1};
cbf{2}      = cbf{1};
P{2}.alpha_v = [0.15, 0.15, 0.25, 0.25, 0.15, 0.15]';
[LBR{2},Y{2}] = LBR_model(P{2},cbf{2});

% Repeat the same simulations but with a constant baseline CBV across
% depths in the AV:
%P{2}     = P{1};
%P{2}.x_v = [1.25, 1.25, 1.75, 1.75, 1.25, 1.25]';
%cbf{2}      = kron([1.6,1.6,1.4,1.4,1.6,1.6],ones(P{1}.T/P{1}.dt,1));
%[LBR{2},Y{2}] = LBR_model(P{2},cbf{2});

% Repeat the same simulations but with different (larger) increase of baseline CBV
% in the AV towards the surface:
P{3}     = P{1};
%P{3}.x_v = [1.25, 1.25, 1.75, 1.75, 1.25, 1.25]';
P{3}.alpha_v = [0.15, 0.15, 0.35, 0.35, 0.15, 0.15]';
%cbf{3}      = kron([1.6,1.6,1.6,1.6,1.6,1.6],ones(P{1}.T/P{1}.dt,1));
cbf{3}  = cbf{1};
[LBR{3},Y{3}] = LBR_model(P{3},cbf{3});

% Display results:
figure(1)

for i = 1:length(LBR)
    subplot(141),
    plot(P{i}.l,flipud(cbf{i}(end,:)'),'.-'); hold on; xlim([0 100]); ylim([0 2.2]); 
    xlabel('1 - Cortical depth (%)'); ylabel('relative CBF (-)'); axis square; title('Laminar CBF profile')
    subplot(142),
    plot(P{i}.l,flipud(Y{i}.V0vq*100),'.-'); hold on; xlim([0 100]); ylim([0 2.2]); 
    xlabel('1 - Cortical depth (%)'); ylabel('Baseline CBV (%)'); axis square; title('Laminar baseline CBV profile')
    subplot(143),
    plot(P{i}.l,flipud(cbf{i}(end,:)').*flipud(Y{i}.F0v),'.-'); hold on; xlim([0 100]); ylim([0 0.8]);
    xlabel('1 - Cortical depth (%)'); ylabel('Absolute CBF (%)'); axis square; title('Laminar Absolute CBF profile')
    subplot(144),
    plot(P{i}.l,flipud(LBR{i}(end,:)'),'.-'); hold on; xlim([0 100]); ylim([0 6]); title('Laminar BOLD profile')
    xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
end;
%legend('neuronal', 'vascular', 'mixed');
legend('0', '0.1', '0.2');
hold off;