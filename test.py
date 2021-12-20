from src.Modeler import *
import matplotlib.pyplot as plt

def uv_load(uv_data, uv_model, filename_ms, telescp, scaling):

    u,v,freqs = uv_data
    uvdist = (u**2 + v**2)**0.5*1e-3
    
    MS = MsReader(filename_ms, telescp)
    _, uv_real, uv_imag, uv_wghts = MS.uvdata_loader()

    if len(scaling)>1:
        if telescp == 'com12m':
            uv_real *= scaling['scale1']
        elif telescp == 'com07m':
            uv_real *= scaling['scale0']
    else:
        uv_real *= scaling['scale0']

    r_uv_real = uv_real - uv_model.real*1e3
    r_uv_imag = uv_imag - uv_model.imag*1e3
    return r_uv_real, r_uv_imag, uv_wghts, uvdist


def bin_it(UVreal, UVimag, uvwghts, uvdist, nbins = 7):
    bins = np.logspace(np.log10(np.nanmin(uvdist)), np.log10(np.nanmax(uvdist)), nbins+1)

    UVrealbinned = np.empty(nbins)
    UVrealerrors = np.empty(nbins)

    UVimagbinned = np.empty(nbins)
    UVimagerrors = np.empty(nbins)

    for i in range(nbins):
        mask = (uvdist>bins[i]) & (uvdist <= bins[i+1])
        UVrealbinned[i], UVrealerrors[i] = np.average(UVreal[mask], weights = uvwghts[mask], returned = True) #weighted average
        UVimagbinned[i], UVimagerrors[i] = np.average(UVimag[mask], weights = uvwghts[mask], returned = True)

        UVrealerrors[i] = UVrealerrors[i]**(-0.5)
        UVimagerrors[i] = UVimagerrors[i]**(-0.5)
        
    return UVrealbinned, UVrealerrors, UVimagbinned, UVimagerrors, bins


def main(filename_samples, filename_ms, telescp): 
    
    fig, ax = plt.subplots(1,2)
    fig.set_figwidth(10)
    fig.set_figheight(4)
    
    uvdist = np.empty(0)
    uvreal = np.empty(0)
    uvimag = np.empty(0)
    uvwght = np.empty(0)
    uvtype = np.empty(0)
    
    for i in range(len(telescp)):

        obj = Modeler(filename_samples, filename_ms[i], telescp[i], save = False)
        uvmodel, scaling, uvdata = obj.run()

        r_uv_real, r_uv_imag, uv_wghts, uv_dist = uv_load(uvdata, uvmodel, filename_ms[i], telescp[i], scaling)

        uvtype = np.append(uvtype, np.ones_like(r_uv_imag)*i)
        uvdist = np.append(uvdist, uv_dist)
        uvreal = np.append(uvreal, r_uv_real)
        uvimag = np.append(uvimag, r_uv_imag)
        uvwght = np.append(uvwght, uv_wghts)

    UVrealbinned, UVrealerrors, UVimagbinned, UVimagerrors, bins = bin_it(uvimag, uvreal, uvwght, uvdist, nbins = 7)        

    ax[0].scatter(uvdist[uvtype == 0], uvreal[uvtype == 0], s = .1, label = telescp[0])
    ax[0].scatter(uvdist[uvtype == 1], uvreal[uvtype == 1], s = .1, label = telescp[1])
    ax[0].errorbar((bins[1:]+bins[:-1])/2, UVrealbinned, xerr = (bins[1:] - bins[0:-1])/2, yerr = UVrealerrors, c = 'C2', ls = '', marker = 'o')
    ax[0].axhline(0, ls = 'dashed', c = 'gray')
    ax[0].axis(ymin = -3, ymax = 5)

    ax[0].set_ylabel('Re(V)$_{data}$ - Re(v)$_{model}$ [mJy]')
    ax[0].set_xlabel('k$\lambda$')
    ax[0].set_xscale('log')
    ax[0].legend()

    ax[1].scatter(uvdist[uvtype == 0], uvimag[uvtype == 0], s = .1, label = telescp[0])
    ax[1].scatter(uvdist[uvtype == 1], uvimag[uvtype == 1], s = .1, label = telescp[1])
    ax[1].errorbar((bins[1:]+bins[:-1])/2, UVimagbinned, xerr = (bins[1:] - bins[:-1])/2, yerr = UVimagerrors, c = 'C2', ls = '', marker = 'o')
    ax[1].axhline(0, ls = 'dashed', c = 'gray')
    ax[1].axis(ymin = -5e0, ymax = 5e-0)
    ax[1].set_ylabel('Im(V)$_{data}$ - Im(V)$_{model}$ [mJy]')
    ax[1].set_xlabel('k$\lambda$') 
    ax[1].set_xscale('log')
    ax[1].legend()

    plt.tight_layout()
    plt.savefig('/scigarfs/home/jvmarrewijk/eszee/plots/vis_residue_pointandSZ_ACA_ALMA.pdf')
    plt.show()
    
if __name__ == '__main__':
    names              = ['RXC_J2014.8_ALMA_pointonly', 
                          'RXC_J2014.8_ALMA_bothpoint_and_SZ', 
                          'RXC_J2014.8_ACA_and_ALMA_null_test_pointsource',
                          'RXC_J2014.8_ACA_and_ALMA_extraspw_up']
    
    name               = names[2]
    filename_samples   = '/scigarfs/home/jvmarrewijk/eszee/outputdir/'+name+'_pickle'
    telescp            = ['com12m', 'com07m']
    filename_ms        = [
        '/scigarfs/home/jvmarrewijk/eszee/RXC_J2014.8_Science_Data/RXC_J2014.8_ALMA_timebin30s_SZ_concat.ms',
        '/scigarfs/home/jvmarrewijk/eszee/RXC_J2014.8_Science_Data/RXC_J2014.8_ACA_timebin30s_SZ_concat.ms'        
    ]
    
    main(filename_samples, filename_ms, telescp)
