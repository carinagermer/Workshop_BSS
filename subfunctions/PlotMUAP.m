function [] = PlotMUAP(Spikes,SIG,fsamp)

    MU_form = cellfun(@(x) cutAPs(Spikes,x,30,0,0),SIG,'Uni',0);

    % figure(), hold on

    maxp2p = max(max(cellfun(@peak2peak,MU_form)));
    % plot_MU_matrix2( MU_form, maxp2p, 'k' )

    tap=(0:length(MU_form{1})-1)/fsamp;

    st_x=1.1;
    st_y=1;

    x=0:tap(end)*st_x:(size(MU_form,2)-1)*tap(end)*st_x;
    y=0:maxp2p*st_y:(size(MU_form,1)-1)*maxp2p*st_y;

    xCell=num2cell(repmat(x,[size(MU_form,1) 1]));
    yCell=num2cell(repmat(flip(y)',[1 size(MU_form,2)]));

    % figure()
    hold on
    cellfun(@(x,y,M) plot(x+tap,M+y,'-','color','k','LineWidth',1.5),xCell,yCell,MU_form)

    xlim = ([0 (size(MU_form,2))*tap(end)*st_x]);
    ylim = ([-maxp2p (size(MU_form,1))*maxp2p*st_y]);

end

