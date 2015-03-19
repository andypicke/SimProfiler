function REsamp=ResampleFieldGeneral(data_real,t_real,z_real,SP)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% function REsamp=ResampleFieldGeneral(data_real,t_real,z_real,SP)
%
% ResampleFieldGeneral.m
%
% Resample a depth-time field of data, simulating sampling by a profiling
% instrument such as a moored profiler or CTD going up and down at a
% specified (constant) speed.
%
% Does an ensemble of phases (ie does CombineREsamp except of having
% different function.)
%
%
% INPUT
% data_real : `Real' data to be resampled (Depth X time) [M X N]
% z_real    : Depths of data_real [M X 1] (m)
% t_real    : Times of data_real [1 X N] (days)
% SP : (optional) Structure with sampling parameters. If not specified, use defaults
%
% OUTPUT
% Structure `REsamp' with fields:
% Nshift        : # of different realizations to cover one profile period
% data_resamp   :[MM X NN X Nshift]
% t_grid        :[Nshift X NN] Average time of each profile
% zmat          :[MM X NN X Nshift]
% tmat          :[MM X NN X Nshift]
%
% Dependencies:
% MakeProfPath.m
%
% 17 Mar. 2015 - A. Pickering - andypicke@gmail.com
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% sampling parameters

%clear z_range dz_samp dztot Tprof Nshift


if ~exist('SP','var')
    SP=struct();
end

%SP.w_samp=0.5;
%SP.dz_samp=10;
%SP.z_range=z_range;
%SP.t_start=t_start;
%SP.time_range=time_range;
%SP.plotit=0;

if ~isfield(SP,'w_samp')
    SP.w_samp=0.25;
end

if ~isfield(SP,'dz_samp')
    SP.dz_samp=nanmean(diff(z_real));
end

if ~isfield(SP,'z_range')
    SP.z_range=[ nanmin(z_real) nanmax(z_real) ] ;
end

if ~isfield(SP,'t_start')
    SP.t_start=nanmin(t_real);
end

if ~isfield(SP,'time_range')
    SP.time_range=floor(range(t_real))-3;
end

if ~isfield(SP,'tshift')
    SP.tshift=(2/60/24)% 2 mins
end

if ~isfield(SP,'plotit')
    SP.plotit=0
end

dztot=range(SP.z_range);
Tprof=dztot/SP.w_samp/86400; % time for a profile (1 up or down) in days
Nshift=round(2*Tprof/SP.tshift)% number of shifts to go through all phases

%%

% Results will be saved to structure 'REsamp'
REsamp=struct();

plotit=0
% make one sampling path to get sizes for arrays
[tmat,zmat]=MakeProfPath(SP.t_start,SP.z_range,SP.dz_samp,SP.time_range,SP.w_samp,SP.plotit);
REsamp.w_samp=SP.w_samp;
REsamp.Nshift=Nshift;
EmptyMat=nan*ones(size(tmat,1),size(tmat,2),Nshift);
REsamp.data_resamp=EmptyMat;
REsamp.data_real=data_real;
REsamp.tsamp=EmptyMat;
REsamp.zsamp=EmptyMat;
REsamp.tgrid=nan*ones(Nshift,size(tmat,2));
REsamp.zreal=z_real;
REsamp.treal=t_real;

clear tmat zmat

%% Now sample for each phase of ensemble

hb=waitbar(0,'resampling data');
makefig=0
% cycle through each case, creating an ensemble of phases
for whcase=1:Nshift
    waitbar(whcase/Nshift,hb)
    
    clear t0 time_range tmat zmat Nprof Nz
    t0=SP.t_start + (whcase-1)*SP.tshift;
    %    [tmat,zmat]=MakeProfPath(t0,z_range,dz_samp,time_range,w_samp,plotit);
    [tmat,zmat]=MakeProfPath(t0,SP.z_range,SP.dz_samp,SP.time_range,SP.w_samp,SP.plotit);
    Nprof=size(zmat,2);
    Nz=size(zmat,1);
    
    % make empty arrays
    clear data_resamp t_grid
    data_resamp=nan*ones(Nz,Nprof);
    t_grid=nan*ones(1,Nprof);
    
    % resample along each profile segment
    for whp=1:Nprof
        clear tvec idt datai
        tvec=tmat(:,whp);
        idt=isin(t_real,[tvec(1) tvec(end)]);
        
        if numel(idt)>2
        datai= interp2(t_real(idt),z_real,data_real(:,idt),tvec,zmat(:,whp));
        
        if nanmean(diff(zmat(:,whp)))<0
            data_resamp(:,whp)=flipud(datai(:));
        else
            data_resamp(:,whp)=datai(:);
        end
        t_grid(whp)=nanmean(tvec);
        
        else
            % no good data, just keep nans
        end
        
    end
    %
    idtall=isin(t_real,[t_grid(1) t_grid(end)]);
    
    if makefig==1
        figure(1);clf
        
        ax1=subplot(211);
        ezpc(t_real(idtall),z_real,data_real(:,idtall))
        hold on
        plot(tmat,zmat,'w-')
        
        ax2=subplot(212);
        ezpc(t_grid,zmat(:,1),data_resamp)
        xlabel('Yearday')
        ylabel('Depth')
        linkaxes([ax1 ax2])
    end
    
    REsamp.data_resamp(:,:,whcase)=data_resamp;
    REsamp.tsamp(:,:,whcase)=tmat;
    REsamp.zsamp(:,:,whcase)=zmat;
    REsamp.tgrid(whcase,:)=t_grid;
    REsamp.z=sort(zmat(:,1));
    
end

delete(hb)
%
REsamp.MakeInfo=['Made ' datestr(now) ' w/ ResampleFieldGeneral.m in vers ' version]
REsamp.SP=SP;

return
%%

% figure(2);clf
% ezpc(REsamp.tgrid(10,:),REsamp.z,REsamp.data_resamp(:,:,10))

%%

% REsamp=AddOverturnsToREsamp(REsamp)

%%
