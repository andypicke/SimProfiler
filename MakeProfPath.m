function [tmat,zmat]=MakeProfPath(t0,z_range,dz_samp,time_range,w_samp,plotit)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%function [tmat,zmat]=MakeProfPath(t0,z_range,dz_samp,time_range,w_samp,plotit)
%
% Make simulated depth-time path for a profiling instrument (i.e. a Moored
% Profiler or CTD) traveling up and down at a fixed speed.
%
% INPUTS
%   t0          : Start time (yday)
%   z_range     : Depth range to sample over (m)
%   dz_samp     : Vertical sampling interval (m)
%   time_range  : Total time (days) of sampling
%   w_samp      : Vertical speed of profiler (m/s)
%   plotit      : Option to plot depth-time path
%
% OUTPUTS
%   tmat, zmat  : Matrices of time and depth along sampling path. Each
%   column is one up or down profile. Time increases down each column,
%   depth can decrease or increase depending on whether profile is up/down.
%
%   To plot sampling path, plot(tmat,zmat)
%
% Orginally modified from part of ResampMPnew12Nov.m
%
% 3 Feb 2015 - A. Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%~~~~~~ Everything below calculated from above chosen parameters

dt_samp=dz_samp/w_samp;

% time to do 1 profile
T_samp=range(z_range)/w_samp; % sec
T_samp_mins=T_samp/60; % min
T_samp_day=T_samp/86400; %days

% how many profiles we want to do (each up or down ==1 profile)
%    Nprof=300;
Nprof=round(time_range/T_samp_day);

aZvec=z_range(1):dz_samp:z_range(end);
%~~ # of points in one profile
Nz=length(aZvec);

dt_samp_day=dt_samp/86400;
%~~ make time vector of sampling for all profiles
t_samp=t0 : dt_samp_day : t0 + dt_samp_day*((Nz*Nprof)-1 );

%t_samp_mat=reshape
zmat=nan*ones(Nz,Nprof);
tmat=nan*ones(Nz,Nprof);
data_resamp=nan*ones(Nz,Nprof);
t_grid=nan*ones(1,Nprof);

for whp=1:Nprof
    clear tvev idt T Z datai
    tvec=t_samp( 1+ (whp-1)*Nz : Nz*whp);
    
    clear zvec
    if iseven(whp)
        zvec=fliplr(aZvec);
    else
        zvec=aZvec;
    end
    
    t_grid(whp)=nanmean(tvec);
    zmat(:,whp)=zvec;
    tmat(:,whp)=tvec;
    
end
%
if plotit==1
    figure(1);clf
    plot(tmat,zmat,'k')
    xlabel('Yearday')
    ylabel('Depth')
    title(['w=' num2str(w_samp) 'm/s'])
    axis ij
    grid on
    xlim([t0 t0+time_range])
    ylim(z_range)
end

return
%%