#!/usr/bin/env python3
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.table import unique, Column
from astropy.io import fits
from astropy import wcs
import os
import glob
import warnings
import sys
warnings.filterwarnings('ignore')

def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

def parse_coord(ra, dec):
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in str(ra) and ':' not in str(dec))):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in str(ra) and ':' in str(dec)):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)

def get_swift_data(ra, dec, radius=30.0 * u.arcmin, discovery_date='2023-10-01',
    max_delta_date=365.25 * u.day):

    heasarc = Heasarc()

    coord = parse_coord(ra, dec)
    mission='swiftuvlog'
    radius = radius.to(u.arcmin).value

    table = heasarc.query_region(coord, mission=mission, radius=f'{radius} arcmin')
    table['OBSID'] = [str(i) for i in table['OBSID']]
    table.sort('START_TIME')

    # Get only unique values of OBSID for downloading
    table = unique(table, keys='OBSID', keep='first')

    # Create new row to decide if data are science or template
    discovery_date = Time(discovery_date)
    dt = TimeDelta(max_delta_date)
    obs_type = []
    for row in table:
        t = Time(row['START_TIME'], format='mjd')

        if t < discovery_date or t > discovery_date + dt:
            obs_type.append('template')
        else:
            obs_type.append('science')

    table.add_column(Column(obs_type,name='OBS_TYPE'))

    return(table)

def download_swift_data(obstable, outdir='.'):

    obstable.sort('OBSID')
    for row in obstable:

        date = Time(row['START_TIME'], format='mjd')
        month = date.datetime.strftime('%m')
        year = date.datetime.strftime('%Y')
        obsid = row['OBSID']

        cmd = 'wget -q -nH --no-clobber --no-check-certificate --cut-dirs=5 -r '
        cmd += '-l0 -c -np -R \'index*\' -erobots=off --retr-symlinks '
        url = f'https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{year}_{month}/{obsid}/uvot/'

        cmd += url

        print(cmd)

        os.system(cmd)


def create_run_files(ra, dec, obstable, outdir='.', phot_radius=5.0 * u.arcsec,
    bkg_radius=10.0*u.arcsec, verbose=False):

    # Create source and background files
    coord = parse_coord(ra, dec)

    phot_radius = phot_radius.to(u.arcsec).value

    sn_file = os.path.join(outdir, 'sn.reg')
    with open(sn_file, 'w') as f:
        ra_hms, dec_dms = coord.to_string(style='hmsdms', precision=2, sep=':').split()
        f.write(f'fk5;circle({ra_hms},{dec_dms},{phot_radius}\")\n')

    bkg_file = os.path.join(outdir, 'bkg.reg')
    with open(bkg_file, 'w') as f:
        ra_hms, dec_dms = coord.to_string(style='hmsdms', precision=2, sep=':').split()
        inner_radius = 2*phot_radius
        outer_radius = 4*phot_radius
        f.write(f'fk5;annulus({ra_hms},{dec_dms},{inner_radius}\",{outer_radius}\")\n')

    coord = parse_coord(ra, dec)
    globstr = os.path.join(outdir,'*','uvot','image','*_sk.img.gz')

    science_file = os.path.join(outdir, 'science')
    template_file = os.path.join(outdir, 'template')

    science = open(science_file,'w')
    template = open(template_file, 'w')

    science_data = []
    template_data = []

    for file in glob.glob(globstr):

        hdu = fits.open(file)

        exptime = hdu[0].header['TSTOP']-hdu[0].header['TSTART']

        w = wcs.WCS(hdu[1].header)

        naxis1 = hdu[1].header['NAXIS1']
        naxis2 = hdu[1].header['NAXIS2']

        x,y = w.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)

        if x<0 or x>naxis1 or y<0 or y>naxis2:
            continue

        obsid = hdu[0].header['OBS_ID']
        mask = obstable['OBSID']==obsid

        obs_type = obstable[mask][0]['OBS_TYPE']

        if obs_type=='science':
            science.write(file+'\n')
            science_data.append({'filter':hdu[0].header['FILTER'].strip(),
                'exptime':exptime,'mjd':Time(hdu[0].header['DATE-OBS']).mjd})

        elif obs_type=='template':
            template.write(file+'\n')
            template_data.append({'filter':hdu[0].header['FILTER'].strip(),
                'exptime':exptime,'mjd':Time(hdu[0].header['DATE-OBS']).mjd})

    if verbose:
        science_data = sorted(science_data, key=lambda x: x['mjd'])

        templates = {}
        for t in template_data:
            filt = t['filter']
            if filt in templates.keys():
                templates[filt]+=t['exptime']
            else:
                templates[filt]=t['exptime']

        for sci in science_data:
            filt = sci['filter']
            if filt not in templates.keys():
                template_exptime = 0.0
            else:
                template_exptime = templates[filt]

            template_fraction = template_exptime/sci['exptime']

            if template_fraction==0.0:
                template_msg = 'NO TEMPLATE'
            elif template_fraction < 1.0:
                template_msg = 'SHALLOW TEMPLATE'
            elif template_fraction > 1.0:
                template_msg = 'GOOD TEMPLATE'

            print('%5.4f'%sci['mjd'],sci['filter'].ljust(5),
                str('%.3f'%sci['exptime']).rjust(10),
                '%.3f'%float(template_exptime/sci['exptime']),
                template_msg)

    science.close()
    template.close()

    return(sn_file, bkg_file, science_file, template_file)


if __name__=='__main__':

    ra = sys.argv[1]
    dec = sys.argv[2]
    discovery_date = sys.argv[3]
    if len(sys.argv)>4:
        max_date = sys.argv[4] * u.day
    else:
        max_date = 365.25 * u.day

    obstable = get_swift_data(ra, dec, discovery_date=discovery_date,
        max_delta_date=max_date)
    download_swift_data(obstable)

    sn,bkg,sci,tmpl = create_run_files(ra, dec, obstable, verbose=True)

    cmd = f'Swift_photom_host.py {sci} {tmpl} -s {sn} -b {bkg} -a -d 3'
    print(cmd)







