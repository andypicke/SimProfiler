function REsamp=AddOverturnsToREsamp(REsamp,Params)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% function REsamp=AddOverturnsToREsamp(REsamp,Params)
%
% AddOverturnsToREsamp.m
%
% Compute overturns for resampled data made in ResampleFieldGeneral.m
%
% INPUT
% REsamp structure  (as made in ResampleFieldGeneral.m)
% Params :(optional) Parameters for overturn calculations
%
% OUTPUT
% REsamp structure with overturn fields added
%
% 17 Mar 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

if ~exist('Params','var')
    %~ default overturns parameters
    minOT=50;
    Params.lat=20.5;
    Params.usetemp=0;
    Params.minotsize=minOT;
    Params.sigma=1e-5;
    Params.runlmin=0;
    Params.plotit=0;
else
end

hb=waitbar(0,'Adding overturns to REsamp')
for whc=1:REsamp.Nshift
    whc;
    waitbar(whc/REsamp.Nshift,hb)
    % Make structure 'xOT' to do overturns
    clear xOT
    xOT.T=REsamp.data_resamp(:,:,whc);
    % assume constant salinity
    xOT.S=34.604-.045*(xOT.T-2.5);
    xOT.z=REsamp.z;
    %
    addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/
    addpath /Users/Andy/Cruises_Research/Overturns_AP/
    
    %~ Now compute overturns for these periods!
    xOT=AddOverturnsToGrid(xOT,Params);
    
    
    REsamp.eps(:,:,whc)=xOT.eps;
    REsamp.Lot(:,:,whc)=xOT.Lot;
    REsamp.Lttot(:,:,whc)=xOT.Lttot;
    REsamp.d(:,:,whc)=xOT.d;
    REsamp.Lot(:,:,whc)=xOT.Lot;
    REsamp.Otnsq_out(:,:,whc)=xOT.Otnsq_out;
    REsamp.n2(:,:,whc)=xOT.n2;
    
    REsamp.Lot_each(:,:,whc)=xOT.Lot_each;
    REsamp.Lt_each(:,:,whc)=xOT.Lt_each;
    REsamp.Otnsq_each(:,:,whc)=xOT.Otnsq_each;
    REsamp.eps_each(:,:,whc)=xOT.eps_each;
    
end % whc

REsamp.OTinfo=['Overturns added ' datestr(now) ' w/ ' mfilename '.m'];
REsamp.OTparams=Params;
delete(hb)

return
%%