function wT = FixedPoint(eYW,b0,fsamp,B,G,DG,PLOT_PAUSE)
% eYW = whiten observations
% b0 = initialization of separation vector
% fsamp = sample frequency
% B = separation matrix
% G = first derivate of contrast function
% DG = second derivate of contrast function

TolX=1e-4; %convergence tolerance
MAXCOUNT=40; %max number of iterations

n = size(eYW,1);
m = size(eYW,2);


% w = eYW(:,ind);         % Initial projection vector
w=b0;
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

    w = w - (B*(B'*w));                         % Deflation (decorrelation)

    w = w/norm(w);                           % normalization
    if PLOT_PAUSE ==1
        figure(1),hold off,plot([0:length(temp)-1]/fsamp,(temp.^2)/norm(temp.^2));
        % hold on,plot(ind/fsamp,temp(ind).^2/norm(temp.^2),'Or');
        ylim([0 0.25]), title(sprintf('iteration %d',counter))
        drawnow
        % pause(0.5)
    end

    counter = counter+1;                        % increment counter
end

wT = real(w);                                   % save latest projecting vector

end