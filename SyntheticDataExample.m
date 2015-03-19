%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% SyntheticDataExample.m
%
% Make a synthetic temperature field with sinusoidal
% oscillations, and resample with simulated MP or CTD at different speeds.
%
%
% 19 Mar 2015 - A. Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% * modify for your path *
addpath /Users/Andy/Cruises_Research/SimProfiler/

% make a synthetic temperature field
t=0 : 2/60/24 : 1 ;

% ~~ Use a linear temperature profile
z=0:10:2000;
tmean=linspace(2,8,length(z));
tmean=flipud(tmean');

[T,Z]=meshgrid(t,z);

% ~~ add a sinusoidal variability
Tsyn=12 % period of oscillations (hr)
om=2*pi*24/Tsyn; % frequency of oscillations
A=1;
temp=nan*ones(size(T));
temp=repmat(tmean,1,length(t))+ A.*sin(om.*T);

% Now simulate sampling

%~ Sampling parameters
clear SP
SP.z_range=[z(1) z(end)];
SP.dz_samp=10;
SP.w_samp=0.5;
SP.t_start=0.0;     % start time
SP.tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
SP.time_range=0.9;

%~ the 'true' data to sample
data_real=temp;
t_real=t;
z_real=z;

REsamp=ResampleFieldGeneral(data_real,t_real,z_real,SP)

figure(1);clf

ax1=subplot(211)
ezpc(T,Z,temp)
hold on
%contour(T,Z,temp,[10:30],'k')
contour(T,Z,temp,tmean(1:10:end),'k')
cb=colorbar
cb.Label.String='Temperature'
xlabel('Time [days]','fontsize',16)
ylabel('Depth [m] ','fontsize',16)
caxis([2 8])

whc=10
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')

ax2=subplot(212)
ezpc(REsamp.tgrid(whc,:),REsamp.z,REsamp.data_resamp(:,:,whc))
hold on
contour(REsamp.tgrid(whc,:),REsamp.z,REsamp.data_resamp(:,:,whc),tmean(1:10:end),'k')
cb=colorbar
cb.Label.String='Temperature'
xlabel('Time [days]','fontsize',16)
ylabel('Depth [m] ','fontsize',16)
caxis([2 8])

linkaxes([ax1 ax2])

%%