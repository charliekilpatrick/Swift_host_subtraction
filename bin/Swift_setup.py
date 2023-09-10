from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from astropy.table import Column
from astropy.time import Time
from astropy import units as u
import numpy as np
import os
import sys
import argparse

def get_uvot_image_table(coord, radius=10 * u.arcmin, verbose=True,
    init_date=57500.0):

    heasarc = Heasarc()
    mission = 'swiftuvlog'
    table = heasarc.query_region(coord, mission=mission, radius='5 arcmin')

    table['START_TIME']=table['START_TIME'].data.astype('float')
    table.sort('START_TIME')

    table['FILTER']=np.array([str(r).strip().upper() for r in table['FILTER']])

    mask = np.array([str(r).strip().upper()=='IMAGE' for r in table['OPERATION_MODE']])
    table = table[mask]

    filtids={'V':'uvv','B':'ubb','U':'uuu','UVW1':'uw1','UVM2':'um2',
        'UVW2':'uw2'}

    rel_date = table['START_TIME']-init_date
    table.add_column(Column(rel_date, name='REL_DATE'))
    table.add_column(Column(['sw'+r['OBSID']+filtids[r['FILTER']]+'_sk.img.gz'
        for r in table], name='FILENAME'))

    table=table['TARGET_ID','OBSID','RA','DEC','START_TIME','EXPOSURE','FILTER',
        'REL_DATE','FILENAME']

    return(table)

def add_sci_tmpl(table, max_date=200.0):

    if 'REL_DATE' not in table.keys():
        print('WARNING: table needs a relative date column (REL_DATE) to '+\
            'assign science and template tags.')
        return(table)

    mask = (table['REL_DATE']>0.0) & (table['REL_DATE']<max_date)
    tags = np.array(['template']*len(table))
    tags[mask]='science'

    table.add_column(Column(tags, name='TAG'))

    return(table)

def download_image(row, verbose=True, dryrun=False):

    baseurl='https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/'
    datestr = Time(row['START_TIME'],format='mjd').datetime.strftime('%Y_%m')
    basefilename = row['FILENAME']
    outfilename=os.path.join('raw',basefilename)

    fullurl = os.path.join(baseurl, datestr, row['OBSID'], 'uvot', 'image',
        basefilename)

    cmds = ['wget','-q','-nH','--no-check-certificate','--no-clobber',
        '--cut-dirs=5','-r','-l0','-c','-np','-R','\'index*\'',
        '-erobots=off','--retr-symlinks','--show-progress',
        '-O',outfilename,fullurl]

    if verbose: print(' '.join(cmds))
    if not dryrun: os.system(' '.join(cmds))

def mk_swift_reduction_files(table, coord, radius=5*u.arcsec, 
    inner=100*u.arcsec, outer=200*u.arcsec):

    tmp_file = open('tmpl.dat','w')
    sci_file = open('sci.dat','w')

    for row in table:
        if row['TAG']=='science':
            sci_file.write(os.path.join('raw',row['FILENAME'])+'\n')
        elif row['TAG']=='template':
            tmp_file.write(os.path.join('raw',row['FILENAME'])+'\n')

    ra_hms, dec_dms = coord.to_string(style='hmsdms',sep=':').split()

    with open('sn.reg', 'w') as f:
        rad = radius.to(u.arcsec).value
        f.write(f'fk5;circle({ra_hms},{dec_dms},{rad}")')

    with open('bkg.reg', 'w') as f:
        inner_rad = inner.to(u.arcsec).value
        outer_rad = outer.to(u.arcsec).value
        f.write(f'fk5;annulus({ra_hms},{dec_dms},{inner_rad}",{outer_rad}")')



if __name__=='__main__':
    params = argparse.ArgumentParser()
    params.add_argument('ra', default=None, help='Right ascension for reduction')
    params.add_argument('dec', default=None, help='Declination for reduction')
    params.add_argument('--init-date', type=str, default='2020-01-01', help='Initial date for transient event to define light curve reduction.')
    params.add_argument('--max-date', type=float, default=200.0, help='Maximum date from initial to use for science reductions.')
    params.add_argument('--verbose', default=False, action='store_true', help='Verbose output for download and analysis.')
    params.add_argument('--dry-run', default=False, action='store_true', help='Do a dry run of download to display available data.')
    args = params.parse_args()

    coord = SkyCoord(args.ra, args.dec, unit='deg')
    init_date = Time(args.init_date).mjd
    table = get_uvot_image_table(coord, init_date=init_date)
    table = add_sci_tmpl(table, max_date=args.max_date)

    for row in table:
        download_image(row, verbose=args.verbose, dryrun=args.dry_run)

    mk_swift_reduction_files(table, coord)

    print('Swift_photom_host.py sci.dat tmpl.dat -s sn.reg -b bkg.reg -a -d 3')



