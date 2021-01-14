function RatesInStates = GetRatesInStates_LR(force)
%
% Lisa Roux Mai 2018
% June 2, 2020: Lisa update - get data from sleep scoring in new data
% format
% July 27: update to automatically retrieve: baseName,namestatefile
%   adds firing rates during whole session
%   deals with cases when there is no epoch in REM, NREM or WAKE
% 
% retrieves firing rate of each units in the different brain states
% (obtained with TheStateEditor)
% Note: uses buzcode format and [baseName,'.spikes.cellinfo.mat']
%
% INPUT
% force: 1 if you want to save FiringRates.mat again
% 
% OUTPUTS
% RatesInStates: result matrix with one row per unit and in column 
% 1- Firing rate in NREM
% 2- Firing rate in REM
% 3- Firing rate in Wake
    % 4- nan for now
    % 5- nan for now
% 6- Firing rate over the entire session 
%
% saves:
% FiringRates.mat in session folder with: RatesInStates totREM totSWS totWAKE
% 
% Dependency: AccumulateTimeInt (included below this function), bz_BasenameFromBasepath
%
% Requires to have done the sleep scoring before with TheStateEditor

%%
 basepath = pwd;
 basename = bz_BasenameFromBasepath(basepath);

%  namestatefile = [basename,'_merged-states.mat'];% Pascal
  namestatefile = [basename,'-states.mat'];
 baseName = basename;
%ex baseName = 'TC03_Intan_S03_08012017';
% ex namestatefile = 'TC03_Intan_S03_08012017_merged-states.mat'


%% 1- Compute duration spent in each brain state

% % namestatefile = 'TC03_Intan_S03_08012017_merged-states.mat'
% cd StateEditor
% statesint = GetStatesInt_LR(namestatefile);
% cd ..

if ~exist('StatesIntervals.mat') % june 2 2020 update Lisa
    statesint = GetStatesInt_LR(namestatefile);
else
    load('StatesIntervals.mat');
end
    % statesint:    structure which stores the intervals for each states
                    % statesint.sws = NREMint;
                    % statesint.wake = WAKint;
                    % statesint.rem = REMint;
                    % statesint.drowzy = Dint;
                    % statesint.inter = Iint;


totREM = AccumulateTimeInt(statesint.rem); % in sec
totSWS = AccumulateTimeInt(statesint.sws);
totWAKE = AccumulateTimeInt(statesint.wake);
% todo : other states

%% 2- Retrieve rates 

if ~exist('FiringRates.mat') || force==1

load([basename,'.sessionInfo.mat'])
load([baseName,'.spikes.cellinfo.mat'])
% load('TC03_Intan_S03_08012017.spikes.cellinfo.mat')
% spikes

Ncells = size((spikes.UID),2);
Nstates = 5 + 1;
Res = nan(Ncells,Nstates);


% NREM
int = statesint.sws;
durState = totSWS;
StateID = 1;

if durState ~=0
    for UID = 1:Ncells
        tsp = spikes.times{UID};
        tspState = Restrict(tsp,int);
        totspkState = size(tspState,1);
        rateState = totspkState/durState;
        Res(UID,StateID)= rateState;
    end
else
    Res(1:Ncells,StateID) = nan;    
end


% REM
int = statesint.rem;
durState = totREM;
StateID = 2;

if durState ~=0
    for UID = 1:Ncells
        tsp = spikes.times{UID};
        tspState = Restrict(tsp,int);
        totspkState = size(tspState,1);
        rateState = totspkState/durState;
        Res(UID,StateID)= rateState;
    end
else
    Res(1:Ncells,StateID) = nan;    
end


% Wake
int = statesint.wake;
durState = totWAKE;
StateID = 3;

if durState ~=0
    for UID = 1:Ncells
        tsp = spikes.times{UID};
        tspState = Restrict(tsp,int);
        totspkState = size(tspState,1);
        rateState = totspkState/durState;
        Res(UID,StateID)= rateState;
    end
else
    Res(1:Ncells,StateID) = nan;    
end

    
% all states
durState = sessionInfo.duration;
StateID = 6;

for UID = 1:Ncells
    tsp = spikes.times{UID};
%     tspState = Restrict(tsp,int);
    totspkState = size(tsp,1);
    rateState = totspkState/durState;
    Res(UID,StateID)= rateState;
end


RatesInStates = Res;

save FiringRates RatesInStates totREM totSWS totWAKE

else
    load('FiringRates.mat')
    disp('Loading existing data from FiringRates')

end
end

function [TotTime] = AccumulateTimeInt(Intervals)

% Intervals should be a list of [start stop] intervals


% Intervals = MovIntPre;

IntNB = length(Intervals(:,1));

TotTime = 0;

for ii = 1:IntNB
    TimeInt = Intervals(ii,2)-Intervals(ii,1);
    TotTime = TotTime + TimeInt;
end

end