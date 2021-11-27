using LinearAlgebra


let


    function ρf(a)
        r0 =3396.2
        ρ0 = 5.25*1e7
        H = 7.295
        ρs = ρ0*exp(-(a)/H)
    end
    function ρf_p(a)
        r0 =3396.2
        ρ0 = 5.25*1e7
        H = 7.295
        ρs = -ρ0*exp(-(a)/H)/H
    end
    function ρf_2p(a)
        r0 =3396.2
        ρ0 = 5.25*1e7
        H = 7.295
        ρs = ρ0*exp(-(a)/H)/H^2
    end

    alts = 1:125

    ρss = ρf.(alts)
    ρs_ps = ρf_p.(alts)
    ρs_2ps = ρf_2p.(alts)
    # @show alts
    # @show ρs
    mat"
    figure
    hold on
    plot($ρss)
    set(gca,'YScale','log')
    hold off
    "

    mat"
    figure
    hold on
    plot($ρs_ps)
    set(gca,'YScale','log')
    hold off
    "

    mat"
    figure
    hold on
    plot($ρs_2ps)
    set(gca,'YScale','log')
    hold off
    "
end
