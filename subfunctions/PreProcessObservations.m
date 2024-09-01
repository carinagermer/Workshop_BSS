function eYW = PreProcessObservations(Data,extFact,Whiten)
% Extend Observations
eYT = ExtendObservations(Data,extFact); 
eYT = eYT - repmat(mean(eYT,2),1,size(eYT,2));

if Whiten
    % % whitening using SVD
    % [U,S,V] = svd(eYT'*eYT/length(eYT));    % SVD correlation matrix
    % FACT = prctile(diag(S),25);             % remove eigenvectors
    % SI = 1./sqrt(diag(S) + FACT);           % inverse eigenvalues
    % WM = U * diag(SI) * V';                 % whitening matrix

    % whitening using covmat
    % covariance matrix
    covmat = eYT*eYT' / (length(eYT)-1);
    % eigendecomposition
    [evecs,evals] = eig(covmat);
    [evals,sidx] = sort(diag(evals),'descend');
    evecs = evecs(:,sidx);
    % Compute the whitening matrix
    FACT = prctile(evals,25);             % remove eigenvectors
    WM = evecs * diag(1 ./ sqrt(evals+FACT)) * evecs';

    eYW = WM * eYT;                         % whitened extended signals
else 
    eYW = eYT;
end

end

