function plot_groundtracks(drhist,crhist,althist,xf_dr,xf_cr,num2plot,id)

    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    for i = 1:length($drhist)
        px = $drhist{i};
        py = $crhist{i};
        if i < ($num2plot +1)
            plot(px,py,'Color',rgb1 + drgb*(i-1)/($num2plot),'linewidth',3)
        end
        %plot(px(1),py(1),'r.','markersize',20)
    end
    plot($xf_dr,$xf_cr,'g.','markersize',20)
    xlabel('downrange (km)')
    ylabel('crossrange (km)')
    hold off
    %saveas(gcf,'range.png')
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    matlab2tikz(strcat($id,'_crdr.tex'))
    %close all
    "

    # this one is for plotting
    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    for i = 1:length($althist)
        px = $drhist{i};
        alt = $althist{i};
        if i < ($num2plot +1)
            colo = drgb*(i-1)/$num2plot;
            plot(px,alt,'Color',rgb1 + colo,'linewidth',3)
        end
        %plot(px(1),alt(1),'r.','markersize',20)
    end
    plot([0,800],ones( 2,1)*10,'r' )
    plot($xf_dr,10,'g.','markersize',20)
    xlim([0 650])
    xlabel('downrange (km)')
    ylabel('altitude (km)')
    hold off
    %saveas(gcf,'alt.png')
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    %matlab2tikz('bankaoa_alt.tex')
    %close all
    "
    return nothing
end
