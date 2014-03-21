%This Matlab script can be used to generate Figure 5 of the paper:
%
%Emil Björnson, Michail Matthaiou, Mérouane Debbah, "A New Look at Dual-Hop
%Relaying: Performance Limits with Hardware Impairments" accepted for
%publication in IEEE Transcations on Communications.
%
%Download article: http://arxiv.org/pdf/1311.2634
%
%This is version 1.2 (Last edited: 2014-03-21).
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

clear all
close all


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));


%Relaying with transceiver hardware impairments

%Select fading parameters with Nakagami-m fading
alphaRange = 1:8; %Range of alpha fading parameters at first and second hop
beta = 1;  %beta fading parameter at first and second hop

R = 2;   %Transmit information rate (per 2 channel uses)
x = 2.^R-1;  %Corresponding outage thresholds

N1 = 1;      %Normalized noise variance at relay
N2 = 1;      %Normalized noise variance at destination

SNR_dB = 20; %SNR range at the strongest hop


mu = [0.2 1 5]; %Different ratios of SNR1/SNR2 in the simulation


kappa1 = 0.1; %Error Vector Magnitude (EVM) at first hop
kappa2 = 0.1; %Error Vector Magnitude (EVM) at second hop
d = kappa1^2+kappa2^2+kappa1^2*kappa2^2; %Recurring function of kappa1 and kappa2


%Placeholders for analytic results: Amplify-and-Forward (AF) case
OP_ideal_af_f = zeros(length(alphaRange),length(mu)); %Outage probabilities, AF fixed gain, ideal hardware
OP_ideal_af_v = zeros(length(alphaRange),length(mu)); %Outage probabilities, AF variable gain, ideal hardware

OP_nonideal_af_f = zeros(length(alphaRange),length(mu)); %Outage probabilities, AF fixed gain, non-ideal hardware
OP_nonideal_af_v = zeros(length(alphaRange),length(mu)); %Outage probabilities, AF variable gain, non-ideal hardware

%Placeholders for analytic results: Decode-and-Forward (DF) case
OP_ideal_df = zeros(length(alphaRange),length(mu));    %Outage probabilities, DF, ideal hardware
OP_nonideal_df = zeros(length(alphaRange),length(mu)); %Outage probabilities, DF, non-ideal hardware


%Run simulations for different values of alpha1 and alpha2
for k = 1:length(alphaRange)
    
    %Extract propagation parameters for first and second hop
    alpha1 = alphaRange(k);
    alpha2 = alphaRange(k);
    beta1 = beta;
    beta2 = beta;
    
    %Extract transmit power for current propagation parameters
    P_dB = SNR_dB - 10*log10(alphaRange(k)*beta); %Corresponding transmit power range at the strongest hop (in dB)
    P = 10.^(P_dB/10); %Corresponding transmit power range at the strongest hop
    
    %Consider different values of the ratio mu=SNR1/SNR2
    for m = 1:length(mu)
        
        %Set the power at each hop, based on the current transmit power at
        %the strongest hop and the ratio mu=SNR1/SNR2
        if mu(m) >= 1
            P1 = P;
            P2 = P / mu(m);
        else
            P1 = P * mu(m);
            P2 = P;
        end
        

        %AF fixed gain relaying with ideal hardware. Outage probability is
        %given by [32]
        G_ideal_f = sqrt(P2./(P1*alpha1*beta1+N1));
        b1 = 0;
        b2 = N1./P1;
        c = N2./(P1.*G_ideal_f.^2);
        d = 0;
        
        OP_ideal_af_f(k,m) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %AF variable gain relaying with ideal hardware. Outage probability is
        %given by [32]
        b1 = N2./P2;
        b2 = N1./P1;
        c = N1*N2./(P1.*P2);
        d = 0;
        
        OP_ideal_af_v(k,m) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %AF fixed gain relaying with non-ideal hardware. Outage probability is
        %given by Theorem 1
        G_nonideal_f = sqrt(P2./(P1*alpha1*beta1*(1+kappa1^2)+N1));
        b1 = 0;
        b2 = (1+kappa2^2)*N1./P1;
        c = N2./(P1.*G_nonideal_f.^2);
        d = kappa1^2+kappa2^2+kappa1^2*kappa2^2;
        
        OP_nonideal_af_f(k,m) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %AF variable gain relaying with non-ideal hardware. Outage probability
        %is given by Theorem 1
        b1 = (1+kappa1^2)*N2./P2;
        b2 = (1+kappa2^2)*N1./P1;
        c = N1*N2./(P1.*P2);
        d = kappa1^2+kappa2^2+kappa1^2*kappa2^2;
        
        OP_nonideal_af_v(k,m) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %DF relaying with ideal hardware. Outage probability is given in [33, Eq. (21)]
        prefactor1 = (N1*x/beta1);
        prefactor2 = (N2*x/beta2);
        part1 = zeros(size(P1));
        part2 = zeros(size(P1));
        
        for j = 0:alpha1-1
            part1 = part1 + exp(-prefactor1./P1) .* (prefactor1./P1).^j / factorial(j);
        end
        
        for j = 0:alpha2-1
            part2 = part2 + exp(-prefactor2./P2) .* (prefactor2./P2).^j / factorial(j);
        end
        
        OP_ideal_df(k,m) = 1 - part1 .* part2;
        
        
        %DF relaying with non-ideal hardware. Outage probability is given by Theorem 2
        prefactor1 = (N1*x/(1-kappa1^2*x)/beta1);
        prefactor2 = (N2*x/(1-kappa2^2*x)/beta2);
        
        part1 = zeros(size(P1));
        part2 = zeros(size(P1));
        
        for j = 0:alpha1-1
            part1 = part1 + exp(-prefactor1./P1) .* (prefactor1./P1).^j / factorial(j);
        end
        
        for j = 0:alpha2-1
            part2 = part2 + exp(-prefactor2./P2) .* (prefactor2./P2).^j / factorial(j);
        end
        
        delta = max([kappa1^2 kappa2^2]);
        
        if delta*x>1
            OP_nonideal_df(k,m) = ones(size(b1));
        else
            OP_nonideal_df(k,m) = 1 - part1 .* part2;
        end
        
    end
    
end


%Plot the results
figure; hold on; box on;

for m = 1:length(mu)
    
    plot(alphaRange,OP_nonideal_af_f(:,m),'bs-.','LineWidth',1);
    plot(alphaRange,OP_ideal_af_f(:,m),'ko-','LineWidth',1);
    
end

set(gca,'yscale','log')

legend('Impairments','Ideal Hardware','Location','SouthWest')

xlabel('Shape Parameter of Gamma Distribution');
ylabel('Outage Probability (OP)');

axis([1 8 1e-6 1])
