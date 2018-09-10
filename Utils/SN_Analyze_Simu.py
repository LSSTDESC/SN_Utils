import numpy as np
import h5py
import argparse
from astropy.table import Table
import pylab as plt
from SN_Rate import SN_Rate

def Plot_Parameters(fieldname,fieldid,tab,season):

    idx = tab['season']==season
    sel=tab[idx]
    thesize=15
    toplot=['X1','Color','DayMax','z']
    fig, ax = plt.subplots(ncols=2, nrows=2,figsize=(10,9))
    fig.suptitle(fieldname+'  - Fieldid '+str(fieldid)+' - season '+str(season+1), fontsize=thesize)
    
    for i,var in enumerate(toplot):
        ix=int(i/2)
        iy=i%2
        axis=ax[ix][iy]
        if var != 'z':
            axis.hist(sel[var])
        axis.set_xlabel(var,fontsize=20)
        axis.set_ylabel('Number of entries',fontsize=thesize)
        axis.tick_params(axis='x', labelsize=thesize)
        axis.tick_params(axis='y', labelsize=thesize)
        if var == 'z':
            n, bins, patches =axis.hist(sel[var])
            bin_center = (bins[:-1] + bins[1:]) / 2
            dz=bins[1]-bins[0]
            sn_rate=SN_Rate(rate='Perrett',H0=72,Om0=0.3,survey_area=9.6)
            zmin=np.min(sel['z'])
            zmax=np.max(sel['z'])
            duration=np.max(sel['DayMax'])-np.min(sel['DayMax'])
            zz,rate, err_rate,nsn, err_nsn=sn_rate(zmin=zmin-dz/2.,zmax=zmax,dz=dz,duration=duration)
            axis.plot(zz,np.cumsum(nsn))
            axis.plot(bin_center,np.cumsum(n))
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Simu file analysis')
    parser.add_argument('--fieldname',
                        default='DD',
                        help='fieldname - DD or WFD')
    parser.add_argument('--fieldid',
                        default=100,
                        help='fieldid of the field to be processed')
    parser.add_argument('--dir_input',
                        default='/home/philippe/LSST/SN_pipelines/SN_Catalog_Simulations/Output_Simu/',
                        help='input files directory')
    parser.add_argument('--season',
                        default=-1,
                        help='Season to process (start at 0/-1 = all)')
    
    args = parser.parse_args()


    #load the hdf5 file
    fieldname=args.fieldname
    fieldid=args.fieldid
    file_name=args.dir_input+'/Simu_'+fieldname+'_'+str(fieldid)+'_'+str(args.season)+'.hdf5'
    f = h5py.File(file_name,'r')
    print(f.keys())
    for i,key in enumerate(f.keys()):
        summary=Table.read(file_name, path=key)
        print(len(summary),summary.dtype)
        for season in np.unique(summary['season']):
            Plot_Parameters(fieldname,fieldid,summary,season)
