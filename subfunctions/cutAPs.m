function [AP, APs, VR] = cutAPs(Pulses,Y,len,PLOT,FILT)
            % extracts cosecutive APs out of signal Y
            % INPUTS:
            %   - Pulses      triggering positions (in samples) for rectangular window used in extraction of APs (firring patterns);
            %   - len           radius of rectangular window (window length = 2*len+1)
            %   - Y             single signal channel (raw vector containing a single channel of a recorded signals)
            %
            % OUTPUTS:
            %   - MUAPs         row-wise matrix of extracted MUAPs (aligned signal intervals of length 2*len+1)

            if isempty(Y)
                AP=nan(1,2*len+1);
                APs=[];
                VR=[];
            else
                c=length(Pulses);

                edgeLen = round(len/2);
                tmp = gausswin(2*edgeLen)';
                win = [tmp(1:edgeLen) ones(1,2*len-2*edgeLen+1) tmp(edgeLen+1:end) ];
                APs=zeros(c-1,1+2*len);


                for k=1:c
                    if FILT
                        APs(k,1:1+2*len)= win .*[ zeros(1,max(Pulses(k)-len,1)-(Pulses(k)-len)) ...
                            (Y(max(Pulses(k)-len,1):min(Pulses(k)+len,length(Y)))) ...
                            zeros(1,Pulses(k)+len-min(Pulses(k)+len,length(Y)))];
                    else
                        APs(k,1:1+2*len)= [ zeros(1,max(Pulses(k)-len,1)-(Pulses(k)-len)) ...
                            (Y(max(Pulses(k)-len,1):min(Pulses(k)+len,length(Y)))) ...
                            zeros(1,Pulses(k)+len-min(Pulses(k)+len,length(Y)))];
                    end
                end


                AP=mean(APs,1);
                GrandMean=mean(AP);

                T=(size(APs,2)); %number of time points
                N=size(APs,1);  %number of AP in the train

                Numerator = sum(sum((APs-repmat(AP,[N 1])).^2))/(T*(N-1)); % mean squared differences between APs
                Denominator =  sum(sum((APs-GrandMean).^2))/(T*N-1);

                VR=Numerator/Denominator;


                if PLOT
                    j=bone(15);
                    plot(APs','color',j(10,:))
                    hold on
                    plot(mean(APs),'color',j(2,:),'LineWidth',2)
                end

            end

end