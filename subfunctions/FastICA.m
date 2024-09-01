function [Spikes,PT,ACT,B] = FastICA(eYW,fsamp,B,G,DG,ACT,MAXCOUNT,ii,ind,PLOTiBSS)


TolX = 1e-5;                    % tolerance fixed point !!!!!!!!!!!!!!!!!

PT_Func = @(x) x.*abs(x);

   %%
    % [~,ind] = max(abs(ACT));
    % ii=ii+1;

    ACT(ind-10:ind+10)=0;

    wT = FixedPointEMG(eYW,ind,fsamp,TolX,MAXCOUNT,B,G,DG,PLOTiBSS);
    B(:,ii) = wT;                                   
  
    s = wT'*eYW;                                    % estimation of the source s
    s = s/max(s);

    PT = PT_Func(s);
    % PT = PT/norm(PT);

    Spikes = peakseek(PT,0.001*fsamp,0.1); %MinDistance = 0.001*fsamp, MinHeight 0.1
    % if PLOTiBSS, hold off, plot(PT,'k'), hold on, scatter(Spikes,PT(Spikes),10,'r','filled'), xlabel('Sample'), end

    %Delete spikes in ACT
    LW = 10e-4;
    auxidx=Spikes;
    auxidx(auxidx<=round(LW*fsamp))=[];
    auxidx(auxidx>=length(s)-round(LW*fsamp))=[];
    for is=1:length(auxidx)
        ACT(auxidx(is)-round(LW*fsamp):auxidx(is)+round(LW*fsamp)) = 0;
    end

end
%%

% [~,ii]=min(ERs);
% 
% wT=B(:,ii);                                   
% 
% s = wT'*eYW;                                    % estimation of the source s
% PT = s.^3;
% PT=PT/norm(PT);
% 
% index = peakseek(PT,0.001*fsamp,0.005); %MinDistance = 0.001*fsamp, MinHeight 0.1
% 
% if ERs(ii)>1
%     wT=mean(eYW(:,index),2);
%     wT = wT/norm(wT);
%     s = wT'*eYW;                                    % estimation of the source s
%     PT = s.^3;
%     PT=PT/norm(PT);
% 
%     index = peakseek(PT,0.001*fsamp,0.005); %MinDistance = 0.001*fsamp, MinHeight 0.1
% 
% end
% 
%  [RoA, FP, FN, ER]=RoA_Error_calculation(RefPulses,index,fsamp,length(eYW)/fsamp,0);
% 
%  fprintf('final: ER: %1.2f, RoA: %1.2f\n', ER,RoA);
% 
% end


%% UTILITIES
function eYW = PreProcessObservations(Data,extFact)
% Extend Observations
eYT = extension(Data,extFact); 
eYT = eYT - repmat(mean(eYT,2),1,size(eYT,2));

% whitening using SVD
[U,S,V] = svd(eYT*eYT'/length(eYT));    % SVD correlation matrix
FACT = prctile(diag(S),25);         % remove eigenvectors
SI = 1./sqrt(diag(S) + FACT);           % inverse eigenvalues
WM = U * diag(SI) * V';                 % whitening matrix

eYW = WM * eYT;                         % whitened extended signals

end

function eY = extension(Y,extfact)
% funzione per estendere le misure
% Y = signals
eY = zeros(size(Y,1)*extfact,size(Y,2)+extfact);
for index = 1:extfact
    eY((index-1)*size(Y,1)+1:index*size(Y,1),[1:size(Y,2)]+(index-1)) = Y;
end

eY=eY(:,1:end-extfact);
end



function wT = FixedPointEMG(eYW,ind,fsamp,TolX,MAXCOUNT,B,G,DG,PLOT_PAUSE)
   

n = size(eYW,1);
m = size(eYW,2);


w = eYW(:,ind);         % Initial projection vector
w = w/norm(w,2);        % normalization

w0 = randn(n, 1);
w0 = w0/norm(w0, 2);

% counter forfixed point
counter = 0;

% fixed point with tolerance TolX or max number of iterations
while abs(abs(w0'*w)-1) > TolX && counter<MAXCOUNT

    w0 = w;                                     % for stopping at TolX

    temp = w'*eYW;                              % estimation source
    w = eYW*G(temp)'/m - sum(DG(temp))*w/m;     %iteration

    %w = w/norm(w, 2);

    w = w - (B*(B'*w));                         % Deflation (decorrelation)

    w = w/norm(w);                           % normalization

    if PLOT_PAUSE ==1
        figure(1),hold off,plot([0:length(temp)-1]/fsamp,(temp.^2)/norm(temp.^2),'k');
        hold on,scatter(ind/fsamp,temp(ind).^2/norm(temp.^2),10,'red','filled');
        ylim([0 0.25]), title(sprintf('Fixed Point Iteration: %d',counter))
        drawnow
                    % pause(0.5)
    end

    counter = counter+1;                        % increment counter
end

wT = real(w);                                   % save latest projecting vector

end

function [locs pks]=peakseek(x,minpeakdist,minpeakh)
% Alternative to the findpeaks function.  This thing runs much much faster.
% It really leaves findpeaks in the dust.  It also can handle ties between
% peaks.  Findpeaks just erases both in a tie.  Shame on findpeaks.
%
% x is a vector input (generally a timecourse)
% minpeakdist is the minimum desired distance between peaks (optional, defaults to 1)
% minpeakh is the minimum height of a peak (optional)
%
% (c) 2010
% Peter O'Connor
% peter<dot>ed<dot>oconnor .AT. gmail<dot>com
if size(x,2)==1, x=x'; end
% Find all maxima and ties
locs=find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>=x(3:end))+1;
if nargin<2, minpeakdist=1; end % If no minpeakdist specified, default to 1.
if nargin>2 % If there's a minpeakheight
    locs(x(locs)<=minpeakh)=[];
end
if minpeakdist>1
    while 1
        del=diff(locs)<minpeakdist;
        if ~any(del), break; end
        pks=x(locs);
        [garb mins]=min([pks(del) ; pks([false del])]); %#ok<ASGLU>
        deln=find(del);
        deln=[deln(mins==1) deln(mins==2)+1];
        locs(deln)=[];
    end
end
if nargout>1,
    pks=x(locs);
end
end


function [ RoA, FPositive, FNegative, ErrorRate ] = RoA_Error_calculation( Spikes1, Spikes2, fsamp,TotalLength, PLOT )
%Mambrito and De Luca 1984
%Negro 2016
%RoA = cj/(cj+Aj+Bj)*100%
%cj: number of discharges of the jth motor unit spike train identified by
%    both decompositions (with tolerance +- 0.5 ms)
%Aj: number of discharges identified only by one of the two decompositions
%Bj: number of discharges identified only by the other decomposition



TolxX = round(0.5e-3 *fsamp); %it has to be high to account on shifts in the action potential waveform

FPositive = ones(size(Spikes2));
FNegative = ones(size(Spikes1));

for k=-TolxX:TolxX
    xcor(k+TolxX+1)=length(intersect(Spikes1,Spikes2+k));
    FPositive(ismember(Spikes2+k,intersect(Spikes1,Spikes2+k)))=0;
    FNegative(ismember(Spikes1,intersect(Spikes1,Spikes2+k)))=0;
end
% figure(), plot(-TolxX:TolxX,xcor)
c=sum(xcor);
% A=length(Spikes1);
% B=length(Spikes2);
A=max([0 length(Spikes1)-c]);
B=max([0 length(Spikes2)-c]);

RoA = c/(c+A+B)*100;

%False-negative rate is defined as the number of false negatives
%divided by the number of spikes in the ground truth recording
% FNrate = sum(FNegative)/length(Spikes1);
%False-positive rate is defined as the number of false positives
%divided by the number of spikes in the spike train
ErrorRate = (sum(FPositive)+sum(FNegative))/(length(Spikes1)+length(Spikes2))*100;
% ErrorRate = mean([FNrate, FPrate])*100; %eLife!
% FPrate = sum(FPositive)/length(Spikes2);
% ErrorRate = (1-sum(~FPositive)/length(Spikes2))*100;

if PLOT
% n=5;
% Spikes1=SNeuronsAll{n,1}; %Our
% Spikes2=SNeuronsAll{n,2}; %Ref
firing1=zeros(1,TotalLength*fsamp);
firing2=zeros(1,TotalLength*fsamp);
firing1(Spikes1)=1;
firing2(Spikes2)=1;
t=(1:length(firing1))/fsamp;

figure('units','inch','position',[0 0 10 2.5]), hold on
plot(t,firing2+1,'g'),plot(t,firing1);
plot(Spikes2(FPositive==1)/fsamp,ones(1,sum(FPositive)),'ro')
plot(double(Spikes1(FNegative==1))/fsamp,ones(1,sum(FNegative)),'ko')
% xlim([0 length(firing2)])


title(sprintf('RoA = %2.0f (Spikes = %d, FP = %d, FN = %d)',RoA, length(Spikes2), sum(FPositive), sum(FNegative)))
% ha(1)=subplot(2,1,1);hold on, plot(firing2), plot(firing1+1,'g');
legend({'Spk1','Spk2'}), ylim([0.9998 1.0002]),%title('Spike Train')
set(gca,'FontSize',12,'FontWeight','bold')
set(gca,'YTick',[]);
xlabel('Time (s)')
% ha(2)=subplot(2,1,2); plot(DATAall(BestChIdx(Same(n,2)),:)), title('Signal'), box off
% linkaxes(ha,'x')
end

end


