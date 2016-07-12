%This Matlab script can be used to generate all figures in the article:
%
%Emil Björnson, Eduard Jorswieck, Mérouane Debbah, Björn Ottersten,
%"Multi-Objective Signal Processing Optimization: The Way to Balance
%Conflicting Metrics in 5G Systems," IEEE Signal Processing Magazine, vol.
%31, no. 6, pp. 14-23, November 2014.
%
%Download article: http://arxiv.org/pdf/1406.2871
%
%This is version 1.0 (Last edited: 2016-07-12)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Initialization
close all;
clear;


%% Define parameters that define resource bundle

%Maximal number of BS antennas
N_max = 500;

%Maximal emitted power per BS antenna (W)
P_max = 20;


%% Define channel conditions

B = 10e6; %Bandwidth (10 MHz)

sigma2 = 10^(-13); %Noise power (W)

B_coherence = 200e3; %Coherence bandwidth (Hz)
tau_coherence = 5e-3; %Coherence time (s)

Upsilon = B_coherence*tau_coherence; %Length of channel coherence in symbols

ISD = 0.25; %Inter-site distance (km)

A = ISD^2; %Area per cell (km^2)

Lambda1 = 1.7166e+09; %Average inverse channel loss
Lambda2 = 0.5419; %Average strength of inter-cell inteference


%% Hardware parameters

%Static hardware power (W)
C_0 = 10;

%Hardware power consumed per transmit antenna (W)
C_N = 0.5;

%Hardware power per user (W)
C_K = 0.2;

%Computational efficiency (flops/W)
L = 12.8e9;

%Efficiency of the power amplifiers at the BS
eta = 0.31;



%% Generate sampling points according to Approach 1

%Prepare to save sample points
UserRate_samples = [];
AreaRate_samples = [];
EE_samples = [];


%Generate ranges of number of antennas and users to consider (note that
%non-integer values are accepted here to get smoother curves)
stepSize = 0.25;
Nrange = 1:stepSize:N_max;
Krange = 1:stepSize:N_max/2;


%Part 1: Fix the user rate and generate different area rates and EE values

%Number of fixed values of the user rate
nbrOfStepsUserRate = 100;

%Generate range of user rate values
UserRate_fixed = linspace(2,97,nbrOfStepsUserRate)*1e6;

%Prepare to save the largest EE for every user rate
EE_fixedUserRate = zeros(size(UserRate_fixed));

%Go through the range of user rates
for n = 1:length(UserRate_fixed)
    
    %Output simulation progress
    disp(['Part 1: ' num2str(n) ' user rates out of ' num2str(length(UserRate_fixed))]);
    
    %Prepare to save new sample points
    EE_values = zeros(length(Nrange),length(Krange));
    AreaRate_values = zeros(length(Nrange),length(Krange));
    
    %Go through range of BS antennas
    for Nind = 1:length(Nrange)
        
        N = Nrange(Nind);
        
        %Go through range of users
        for Kind = 1:length(Krange)
            
            K = Krange(Kind);
            
            %Check if K has a feasible value
            if K<=min([N/2,Upsilon])
                
                %Compute cost of precoding
                C_precoding = 3*K^2*N*B/Upsilon;
                
                %Compute transmit required to attain the current user rate
                P = (2^(UserRate_fixed(n)/(B*( 1 - K/Upsilon))) - 1)*(sigma2*Lambda1) / ((N/K-1) - Lambda2*(2^(UserRate_fixed(n)/(B*( 1 - K/Upsilon))) - 1));
                
                %Check that the power value is feasible and if so, store
                %the corresponding values
                if P<=P_max*N && P>=0
                    EE_values(Nind,Kind) = K*UserRate_fixed(n)/( P/eta + N*C_N + K*C_K + C_precoding/L + C_0);
                    AreaRate_values(Nind,Kind) = K*UserRate_fixed(n)/A;
                end
                
            end
            
        end
        
    end
    
    %Normalize the sample points, keep the feasible ones and sort them
    Z = [AreaRate_values(:)/1e9 EE_values(:)/1e6 ];
    Z = Z(Z(:,1)>0,:);
    Zsort = sortrows(Z,-1);
    
    %Go through the sample points and remove those that are dominated by
    %other sample points (i.e,. cannot be at the Pareto boundary)
    ind = 1;
    
    while ind<length(Zsort)
        
        keep = Zsort(ind,2)<Zsort(:,2);
        keep(1:ind) = true;
        
        Zsort = Zsort(keep,:);
        
        ind = ind + 1;
        
    end
    
    
    %Save the remaining sample points
    UserRate_samples = [UserRate_samples; ones(size(Zsort,1),1)*UserRate_fixed(n)/1e6];
    AreaRate_samples = [AreaRate_samples; Zsort(:,1)];
    EE_samples = [EE_samples; Zsort(:,2)];
   

    %Extract the maximum EE, to be used in Figure 6
    EE_fixedUserRate(n) = max(Zsort(:,2))*1e6;
    
end



%Part 2: Fix the area rate and generate different user rates and EE values

%Number of fixed values of the user rate
nbrOfStepsAreaRate = 50;

%Generate range of sum rate values (obtained from area rate by multiplying
%with the area). 
SumRate_fixed = linspace(2,49,nbrOfStepsAreaRate)*A*1e9;

%Prepare to save the largest EE for every are rate
EE_fixedAreaRate = zeros(size(SumRate_fixed));

%Go through the range of area rates
for n = 1:length(SumRate_fixed)
    
    %Output simulation progress
    disp(['Part 2: ' num2str(n) ' area rates out of ' num2str(length(SumRate_fixed))]);
    
    %Prepare to save new sample points
    EE_values = zeros(length(Nrange),length(Krange));
    UserRate_values = zeros(length(Nrange),length(Krange));
    
    %Go through range of BS antennas
    for Nind = 1:length(Nrange)
        
        N = Nrange(Nind);
        
        %Go through range of users
        for Kind = 1:length(Krange)
            
            K = Krange(Kind);
            
            %Check if K has a feasible value
            if K<=min([N/2,Upsilon])
                
                %Compute cost of precoding
                C_precoding = 3*K^2*N*B/Upsilon;
                
                %Compute transmit required to attain the current area rate
                P = (2^(SumRate_fixed(n)/(K*B*( 1 - K/Upsilon))) - 1)*(sigma2*Lambda1) / ((N/K-1) - Lambda2*(2^(SumRate_fixed(n)/(K*B*( 1 - K/Upsilon))) - 1));
                
                %Check that the power value is feasible and if so, store
                %the corresponding values
                if P<=P_max*N && P>=0
                    EE_values(Nind,Kind) = SumRate_fixed(n)/( P/eta + N*C_N + K*C_K + C_precoding/L + C_0);
                    UserRate_values(Nind,Kind) = SumRate_fixed(n)/K;
                end
                
            end
            
        end
        
    end
    
    %Normalize the sample points, keep the feasible ones and sort them
    Y = [UserRate_values(:)/1e6 EE_values(:)/1e6];
    Y = Y(Y(:,1)>0,:);
    Ysort = sortrows(Y,-1);
    
    %Go through the sample points and remove those that are dominated by
    %other sample points (i.e,. cannot be at the Pareto boundary)
    ind = 1;
    
    while ind<length(Ysort)
        
        keep = Ysort(ind,2)<Ysort(:,2);
        keep(1:ind) = true;
        
        Ysort = Ysort(keep,:);
        
        ind = ind + 1;
    end
    
    
    %Save the remaining sample points
    UserRate_samples = [UserRate_samples; Ysort(:,1)];
    AreaRate_samples = [AreaRate_samples; ones(size(Ysort,1),1)*SumRate_fixed(n)/1e9/A];
    EE_samples = [EE_samples; Ysort(:,2)];
    
    %Extract the maximum EE, to be used in Figure 7
    EE_fixedAreaRate(n) = max(Ysort(:,2))*1e6;
    
end



%% Generate Figure 6

%Compute utopia point
utopiaPoint = [max(UserRate_fixed) max(EE_fixedUserRate)];

%Compute max weighted sum, max weighted product, and max weighted Chebyshev
[~,maxSumInd] = max(UserRate_fixed/utopiaPoint(1)+EE_fixedUserRate/utopiaPoint(2));
[~,maxProdInd] = max(UserRate_fixed.*EE_fixedUserRate);
[~,maxChebyInd] = max(min([UserRate_fixed/utopiaPoint(1); EE_fixedUserRate/utopiaPoint(2)],[],1));

%Plot the figure
figure(6); hold on; box on;
plot(utopiaPoint(1)/1e6,utopiaPoint(2)/1e6,'k.','MarkerSize',30);
plot(UserRate_fixed(maxSumInd)/1e6,EE_fixedUserRate(maxSumInd)/1e6,'ko','MarkerSize',8);
plot(UserRate_fixed(maxProdInd)/1e6,EE_fixedUserRate(maxProdInd)/1e6,'ks','MarkerSize',8);
plot(UserRate_fixed(maxChebyInd)/1e6,EE_fixedUserRate(maxChebyInd)/1e6,'k*','MarkerSize',8);
plot([0 UserRate_fixed UserRate_fixed(end)]/1e6,[0 EE_fixedUserRate 0]/1e6,'k','LineWidth',1);

xlabel('Average User Rate [Mbit/s/user]');
ylabel('Energy Efficiency [Mbit/Joule]');
legend('Utopia point','Max Weighted sum','Max Weighted product','Max Weighted Chebyshev','Location','Best');



%% Generate Figure 7

%Compute utopia point
utopiaPoint = [max(SumRate_fixed) max(EE_fixedAreaRate)];

%Compute max weighted sum, max weighted product, and max weighted Chebyshev
[~,maxSumInd] = max(SumRate_fixed/utopiaPoint(1)+EE_fixedAreaRate/utopiaPoint(2));
[~,maxProdInd] = max(SumRate_fixed.*EE_fixedAreaRate);
[~,maxChebyInd] = max(min([SumRate_fixed/utopiaPoint(1); EE_fixedAreaRate/utopiaPoint(2)],[],1));

%Plot the figure
figure(7); hold on; box on;
plot(utopiaPoint(1)/1e9/A,utopiaPoint(2)/1e6,'k.','MarkerSize',30);
plot(SumRate_fixed(maxSumInd)/1e9/A,EE_fixedAreaRate(maxSumInd)/1e6,'ko','MarkerSize',8);
plot(SumRate_fixed(maxProdInd)/1e9/A,EE_fixedAreaRate(maxProdInd)/1e6,'ks','MarkerSize',12);
plot(SumRate_fixed(maxChebyInd)/1e9/A,EE_fixedAreaRate(maxChebyInd)/1e6,'k*','MarkerSize',8);
plot([0 SumRate_fixed SumRate_fixed(end)]/1e9/A,[0 EE_fixedAreaRate 0]/1e6,'k','LineWidth',1);

xlabel('Average Area Rate [Gbits/km^2]');
ylabel('Energy Efficiency [Mbit/Joule]');
legend('Utopia point','Max Weighted sum','Max Weighted product','Max Weighted Chebyshev','Location','Best');




%% Generate Figure 8

%Create the mesh grid of (user rate, area rate) points to plot
[xq,yq] = meshgrid(min(UserRate_samples):2:max(UserRate_samples), min(AreaRate_samples):2:max(AreaRate_samples));

%Interpolate the point cloud generated by Approach 1
F = TriScatteredInterp(UserRate_samples,AreaRate_samples,EE_samples);

%Compute the EE at the grid points
vq = F(xq,yq);

%Plot the figure
figure(8)
surf(yq,xq,vq);
colormap(hot);

xlabel('Average Area Rate [Gbits/km^2]');
ylabel('Average User Rate [Mbits/user]');
zlabel('Energy Efficiency [Mbits/Joule]');
view(-110,15);
