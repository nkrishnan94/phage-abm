import numpy as np

def cm2inch(size):
    return size/2.54

def linear_reg(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    X = np.vander(x, 2)
    coeffs, res, rank, sing_vals = np.linalg.lstsq(X, y)
    mx = x.sum()/len(x)
    sx = float(((x - mx)**2).sum())
    if len(x) > 2:
        r2 = 1. - res/(y.size*np.var(y))
    else:
        r2 = 0
    return coeffs, r2

def fit_heterozygosity(x_het, het, surv):
    X_INIT = min(min(range(len(het)), key=lambda i:abs(het[i] - 0.1)), 1000) #find initial fitting point
    X_FIN = 0
    while (surv[X_FIN] > 0.05 and X_FIN < len(het) - 1):
        X_FIN += 1 #find final fitting point

    if (X_FIN - X_INIT <= 400):
        X_INIT = max(0, X_FIN - 400)

    X_THIRD = int((X_FIN - X_INIT)/3)
    xh = x_het[X_INIT:X_FIN]
    het_fit = het[X_INIT:X_FIN]
    coeffs, r_sq = linear_reg(xh, np.log(het_fit))#Do fit
    return xh, het_fit, coeffs, r_sq


def make_het_comparison(het_data, ax, label, axis_fontsize):
    fstr_comparison = np.load('data/fstr_neff_det_comparison.npy')
    het_data = np.flipud(het_data)

    ax.set_yscale('log')
    ax.ticklabel_format(style='sci', scilimit=(-2,2), axis='x')
    ax.set_xlabel('time, t', fontweight='bold', fontsize=axis_fontsize)
    ax.set_ylabel('heterozygosity, H', fontweight='bold', fontsize=axis_fontsize)
    ax.text(-3000, 2.1E0, label, fontweight='bold', fontsize=12)
    ax.text(13000, 3E-1, '$H \sim e^{-\Lambda t}$', fontsize=10, color='k')
    ax.set_xticks([0, 10000, 20000])
    ax.set_yticks([1E-1, 1])
    ax.set_ylim([1E-1, 1.8E0])

    number_of_points = 20
    reduced_het_indices = (1000/number_of_points)*np.arange(number_of_points)

    for elem in het_data:
        deme = elem[0][0]
        fstr = elem[0][1]
        x_het = elem[1].T[0]
        het = elem[1].T[1]
        surv = elem[1].T[2]

        x_plot = [x_het[i] for i in reduced_het_indices]
        het_plot = [het[i] for i in reduced_het_indices]
                
        xh, het_fit, coeffs, r_sq = fit_heterozygosity(x_het, het, surv)
        fit = np.poly1d(coeffs)
        x_fit = [x_het[50], 1.1*x_het[-1]] #choose range for plotting fit
        est = np.exp(fit(x_fit))

        #Plot results
        if fstr > -0.5:
            clr = 'b'
            lbl = 'pushed'
        else:
            clr = 'r'
            lbl = 'pulled'
        ax.scatter(x_plot, het_plot, s=20, edgecolor=clr, facecolor='none', lw=1, label=lbl)
        ax.plot(x_fit, est, c=clr, lw=1)
    ax.set_xlim([0, 22000])

    legend_properties={'weight':'normal', 'size':6}
    ax.legend(loc='best', prop=legend_properties, scatterpoints=1)