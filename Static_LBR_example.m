% Havlicek et Uludag, 2020. Figure 3 and Figure 4
% Steady-state example demonstrating linear scaling of depth-dependent BOLD
% response with 'stimulus strenght' or increase of relarive CBF;
close all; clear all;

set(0,'DefaultAxesFontSize', 14, ...
      'defaultLineLineWidth', 2, ...
      'defaultLineMarkerSize',15,...
      'DefaultAxesTitleFontWeight', 'normal');      

K = 6;                       % We will consider six depths
P{1} = LBR_parameters(K);    % Get parameter structure for LBR model
                             % By default we consider 40 sec stimulus
                             % duration in order to reach steady-state (i.e P{1}.T = 40)
% Define laminar profile relative CBF (considering six depth) constant at all layers 
cbf_original   = [1 1 1 1 1 1];
% define total amount of baseline CBV in the GM (in mL), default = 2.5 mL
% define fraction of microvasculature (venules) with respect to the AV,
% default = 0.5
P{1}.s_v = 0;
P{1}.s_d  = 0.3;  % increase of baseline CBV in the AV towards the surface
P{1}.E0v = 0.35;
 
% Call the LBR model to generate the laminar BOLD profile:
%[LBR{1},Y{1}] = LBR_model(P{1},cbf{1});


% Repeat the same simulations but with a constant baseline CBV across
% depths in the AV:
P{2}          = P{1};
%P{2}.s_v = 0;
P{2}.s_d  = 1.2;
%P{2}.w_v = 0.333;
%cbf{2}          = kron(cbf_original*1.6,ones(P{1}.T/P{1}.dt,1));
%[LBR{2},Y{2}] = LBR_model(P{2},cbf{2});


% Repeat the same simulations but with a constant baseline CBV across
% depths in the AV:
%P{3}          = P{1};
%P{3}.s_v = 0.3;
%P{3}.s_d  = 0;
%P{3}.w_v = 0.666;
%cbf{3}          = kron(cbf_original*1.6,ones(P{1}.T/P{1}.dt,1));
%[LBR{3},Y{3}] = LBR_model(P{3},cbf{3});

% Display results:
figure(1)

%{
%% Figure 4
cbf_case = [0.2, 0.4, 0.6, 0.8];
for i = 1:4
    for j = 1:K
        vec_add = zeros(1,6);
        vec_add(1,j) = cbf_case(i);
        cbf_add = kron(cbf_original.*vec_add,ones(P{1}.T/P{1}.dt,1));
        cbf{j} = cbf_original + cbf_add;
        [LBR{1},Y{1}] = LBR_model(P{1},cbf{j});
        lbr_case1(j,:) = flipud(LBR{1}(end,:)');
        [LBR{2},Y{2}] = LBR_model(P{2},cbf{j});
        lbr_case2(j,:) = flipud(LBR{2}(end,:)');
    end
    plot(P{1}.l,lbr_case1(:,:)','x-'); hold on; xlim([0 100]); ylim([0 3]); title('Laminar BOLD profile')
    set(gca, 'ColorOrderIndex', 1);
    plot(P{2}.l,lbr_case2(:,:)','.--'); hold on; xlim([0 100]); ylim([0 3]); title('Laminar BOLD profile')
    xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
    scale(1) = LBR{1}(end,1)./LBR{1}(end,end); % ratio of layer 1 to layer 6
    legend({'Case1','','','','','','Case2'})
    hold off;
    
    for j = 2:K
        pttr_case1 = max(lbr_case1, [], 2)./mean(lbr_case1(j,K-j+2:end));
        pttr_case2 = max(lbr_case2, [], 2)./mean(lbr_case2(j,K-j+2:end));
    end
    
    pttr_avg(1, i) = mean(pttr_case1);
    pttr_avg(2, i) = mean(pttr_case2);
end

figure(2)
plot(P{1}.l(2:end-1),pttr_avg(1,:)','x-'); hold on; xlim([0 100]); ylim([0 3.5]); title('Average PTT ratio (-)')
plot(P{1}.l(2:end-1),pttr_avg(2,:)','.--'); hold on; xlim([0 100]); ylim([0 3.5]); title('Average PTT ratio (-)')
xlabel('Relative CBF(%)'); ylabel('Average PTT ratio (-)'); axis square;
legend({'Case1','Case2'})
%}

%% Figure 3
%{
for i = 1:length(LBR)
    subplot(131),
    plot(P{1}.l,flipud(cbf{i}(end,:)'),'.-'); hold on; xlim([0 100]); ylim([1 2]); 
    xlabel('1 - Cortical depth (%)'); ylabel('relative CBF (-)'); axis square; title('Laminar CBF profile')
    subplot(132),
    plot(P{i}.l,flipud(Y{i}.V0dq*100),'.-'); hold on; xlim([0 100]); ylim([0 3]); 
    xlabel('1 - Cortical depth (%)'); ylabel('Baseline CBV (%)'); axis square; title('Laminar baseline CBV profile')
    subplot(133),
    plot(P{i}.l,flipud(LBR{i}(end,:)'),'.-'); hold on; xlim([0 100]); ylim([0 3]); title('Laminar BOLD profile')
    xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
    scale(i) = LBR{i}(end,1)./LBR{i}(end,end); % ratio of layer 1 to layer 6
end;
legend({'Case1','Case2','Case3'})
hold off;
%}

%% Other
%{
figure(2)
bar(scale), title('ratio between upper and lower depth'); % if flat across different stimulus strengths
ylim([0 3])% it suggests linear

figure(3) % relation ship between CBF and BOLD is still nonlinear
          % we also know that the relatinship between CBF and neuronal
          % response is often nonlinear (especially for more sustained responses)
          % combination of these two nonlinearities often leads to
          % observation that relationship between neuronal and BOLD
          % response is more or less linear.
plot([cbf{1}(end,1),cbf{2}(end,1),cbf{3}(end,1),cbf{4}(end,1)]',...
     [LBR{1}(end,1),LBR{2}(end,1),LBR{3}(end,1),LBR{4}(end,1)]); xlim([1 2]); hold on; % lower layer
plot([cbf{1}(end,end),cbf{2}(end,end),cbf{3}(end,end),cbf{4}(end,end)]',...
      [LBR{1}(end,end),LBR{2}(end,end),LBR{3}(end,end),LBR{4}(end,end)]); xlim([1 2]); hold off % upper layer
 xlabel('relative CBF (-)'), ylabel('BOLD (%)'); title('CBF vs BOLD'); legend({'Upper','Lower'})
%}