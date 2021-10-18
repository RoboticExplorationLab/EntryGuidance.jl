function plot_groundtracks(drhist,crhist,althist,xf_dr,xf_cr,num2plot,id)

    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    cmap = [];
    for i = 1:round($num2plot)
        px = $drhist{i};
        py = $crhist{i};
        plot(px,py,'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
        cmap = [cmap;rgb1 + drgb*(i)/($num2plot)]
        %plot(px(1),py(1),'r.','markersize',20)
    end
    plot($xf_dr,$xf_cr,'g.','markersize',20)
    colormap(cmap);
    pp = colorbar;
    pp.Ticks = 0:(1/$num2plot):1;
    pp.Location = 'northoutside';
    pp.TickLabels = 0:round($num2plot);
    xlabel('downrange (km)')
    ylabel('crossrange (km)')
    xlim([250,700])
    ylim([0,20])
    hold off
    fleg = legend('figure()');
    set(fleg,'visible','off')
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    %matlab2tikz(strcat('cpeg_examples/bank_angle/tikz/',$id,'_crdr.tex'))
    %close all
    "

    # this one is for plotting
    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    cmap = [];
    for i = 1:round($num2plot)
        px = $drhist{i};
        alt = $althist{i};
        plot(px,alt,'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
        cmap = [cmap;rgb1 + drgb*(i)/($num2plot)]
    end
    plot($xf_dr,$xf_cr,'g.','markersize',20)
    colormap(cmap);
    pp = colorbar;
    pp.Ticks = 0:(1/$num2plot):1;
    pp.Location = 'northoutside';
    pp.TickLabels = 0:round($num2plot);
    plot([0,800],ones( 2,1)*10,'r' )
    plot($xf_dr,10,'g.','markersize',20)
    xlim([400,700])
    ylim([8,30])
    xlabel('downrange (km)')
    ylabel('altitude (km)')
    hold off
    %saveas(gcf,'alt.png')
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    fleg = legend('figure()');
    set(fleg,'visible','off')
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    matlab2tikz(strcat('cpeg_examples/bank_angle/tikz/',$id,'_altdr.tex'))
    %matlab2tikz('bankaoa_alt.tex')
    %close all
    "
    return nothing
end
