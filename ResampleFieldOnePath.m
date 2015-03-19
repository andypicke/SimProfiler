function [data_resamp, t_grid]=ResampleFieldOnePath(t_real,z_real,data_real,tmat,zmat,makefig)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% function [data_resamp, t_grid]=ResampleFieldOnePath(t_real,z_real,data_real,tmat,zmat,makefig)
%
% Sample field 'data_real' along a depth-time path specified by tmat,zmat
% (made in MakeProfPath.m).
%
% INPUT
%
% data_real : Data to resample (depth X time gridded), [M X N]
% t_real    : Time of data_real, [1 X N]
% z_real    : Depths of data_real, [M X 1]
% tmat      : Times of sampling path, [MM X NN]
% zmat      : Depths of sampling paths, [MM X NN]
% makefig   : Option to make a plot
%
% OUTPUT
%
% data_resamp : Resampled data , [MM X NN]
% t_grid  : 'Gridded' time of resampled profiles (average time of each profile) [1 X NN]
%
%
% 19 Mar. 2015 - A. Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

[Nz,Nprof]=size(zmat);

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

%%