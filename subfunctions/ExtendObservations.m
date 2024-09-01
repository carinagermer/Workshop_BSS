function eY = ExtendObservations(Y,extfact)
% Y = signals
eY = zeros(size(Y,1)*extfact,size(Y,2)+extfact);
for index = 1:extfact
    eY((index-1)*size(Y,1)+1:index*size(Y,1),(1:size(Y,2))+(index-1)) = Y;
end

eY=eY(:,1:end-extfact);
end