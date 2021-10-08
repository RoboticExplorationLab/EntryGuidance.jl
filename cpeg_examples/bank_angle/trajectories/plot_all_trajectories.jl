
let

# l1_traj = jldopen("cpeg_examples/bank_angle/trajectories/L1_bank.jld2")
l1_traj = load_object("cpeg_examples/bank_angle/trajectories/L1_bank.jld2")
l2_traj = load_object("cpeg_examples/bank_angle/trajectories/L2_bank.jld2")
fl_traj = load_object("cpeg_examples/bank_angle/trajectories/full_lift.jld2")


drhist = [l1_traj.dr, l2_traj.dr, fl_traj.dr]
crhist = [l1_traj.cr, l2_traj.cr, fl_traj.cr]
althist = [l1_traj.alt, l2_traj.alt, fl_traj.alt]
xf_dr = l1_traj.xf_dr
xf_cr = l1_traj.xf_cr

@show norm(l1_traj.dr - l2_traj.dr)
mat"
figure
hold on
rgb1 = [29 38 113]/255;
rgb2 = 1.3*[195 55 100]/255;
drgb = rgb2-rgb1;
for i = 1:length($drhist)
    px = $drhist{i};
    py = $crhist{i};
    plot(px,py,'linewidth',3)
end
plot($xf_dr,$xf_cr,'g.','markersize',20)
xlabel('downrange (km)')
ylabel('crossrange (km)')
xlim([250,650])
ylim([0,11.5])
hold off
legend('L1 bank angle','L2 bank angle','AoA + bank angle','target','location','northwest')
addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
matlab2tikz(strcat('cpeg_examples/bank_angle/tikz/all_crdr.tex'))
close all
"

# this one is for plotting
mat"
figure
hold on
for i = 1:length($althist)
    px = $drhist{i};
    alt = $althist{i};
    plot(px,alt,'linewidth',3)
end
plot([0,800],ones( 2,1)*10,'r' )
plot($xf_dr,10,'g.','markersize',20)
xlim([0 650])
xlabel('downrange (km)')
ylabel('altitude (km)')
hold off
xlim([400,650])
ylim([8,30])
%saveas(gcf,'alt.png')
addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
legend('L1 bank angle','L2 bank angle','AoA + bank angle','target altitude','location','northeast')
addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
matlab2tikz(strcat('cpeg_examples/bank_angle/tikz/all_altdr.tex'))
%matlab2tikz('bankaoa_alt.tex')
%close all
"

end
