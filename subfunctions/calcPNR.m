function PNR = calcPNR(MUPulses,IPT,fsamp)

if ~isempty(MUPulses)
    
    % get time instants ranging 2 samples before and after each discharge instant
    tmpMUPulses = [];
    for t = -2:2
        tmpMUPulses = [tmpMUPulses MUPulses+t];
    end
    
    nInd = setdiff(MUPulses(1):MUPulses(end),tmpMUPulses); %only select instants between the first and last discharges
    IPT = IPT/mean(IPT(MUPulses)); %normalize by the average amplitude of the peaks
    
    tmpn = IPT(nInd); %get IPT in the section between first and last discharges
    tmpn = tmpn(~isnan(tmpn)); %remove NaN
    tmpn = tmpn(tmpn>=0); % get only values > 0
    PNR = round(10*10*log10(  mean(IPT(MUPulses).^2) / mean(tmpn.^2)))/10;
else
    PNR = 0;
end
