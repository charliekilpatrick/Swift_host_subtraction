import os,json
import astropy.io.fits as pf
import SwiftPhotom.errors
import SwiftPhotom.commands as sc
import numpy as np
import matplotlib.pyplot as plt

ZP={'V':[17.88,0.01],'B':[18.98,0.02],'U':[19.36,0.02],'UVW1':[18.95,0.03],'UVM2':[18.54,0.03],'UVW2':[19.11,0.03]}
Vega={'V':-0.01,'B':-0.13,'U':1.02,'UVW1':1.51,'UVM2':1.69,'UVW2':1.73}

mjdref = 51910.0  #  jd reference used by swift

def sort_filters(_filt_list):
    '''
    Function to recognize filter-related keywords
    '''
    
    full_filter_list=['V','B','U','UVW1','UVM2','UVW2']
    if _filt_list=='ALL':
        return full_filter_list
    elif _filt_list=='OPT':
        return ['V','B','U']
    elif _filt_list=='UV':
        return ['UVW1','UVM2','UVW2']
    else:
        out_filter_list=[]
        for _f in _filt_list.split(','):
            if _f.upper() not in full_filter_list:
                print('WARNING - Filter %s not recognized. Skipped.\n' % _f)
                continue
            out_filter_list.append(_f.upper())
        
        if len(out_filter_list)==0:
            raise SwiftPhotom.errors.FilterError
        
        return out_filter_list

def load_obsid(_obsid_string):
    '''
    Looking for sky frames with the given ObsID.
    The search is done in the working directory
    and in all subdirectories.
    It will look for file with the conventional
    naming
    
    sw[obsID][obsIdx]u[filter]_sk.img.gz
    
    or without the .gz
    '''

    #in the file path the ObsID have 8 digits, with leading zeros
    obsid=_obsid_string.zfill(8)

    out_file_list=[]

    for root, dirs, files in os.walk("."):
        for file in files:
            if file.startswith('sw'+obsid) and (file.endswith('_sk.img.gz') or file.endswith('_sk.img')):
                sky_frame = os.path.join(root,file)
                #skipping files in the products folder
                if 'products' not in os.path.normpath(sky_frame).split(os.sep):
                    out_file_list.append(sky_frame)
    
    return out_file_list


def interpret_infile(_infile):
    '''
    Interpret what type of input the user is providing
    
    _infile will be a list of strings with either one
    or two elements
    
    The output will be file_list, a list containing
    2 lists:
    file_list[0] will be the object file list
    file_list[1] will be the template file list
    
    Creating a single list, I can fill them both
    with a for loop, even if only the object list is
    provided.
    '''
    
    file_list=[[],[]]
    
    for i in range(len(_infile)):
        if os.path.isfile(_infile[i]):
            file_exist=1
        else:
            file_exist=0
        
        if file_exist:
            #single file is provided
            if _infile[i].endswith('_sk.img.gz') or _infile[i].endswith('_sk.img'):
                file_list[i].append(_infile[i])
            else:
                with open(_infile[i]) as inp:
                    for line in inp:
                        ff=line.strip('\n')
                        if os.path.isfile(ff):
                            file_list[i].append(ff)
                        else:
                            list_from_obsid = load_obsid(ff)
                            if len(list_from_obsid) == 0:
                                print(ff+' not found. Skipped.')
                            else:
                                file_list[i] = file_list[i] + list_from_obsid
                
                if len(file_list[i])==0:
                    raise SwiftPhotom.errors.ListError(_infile[i])

        #If no file exists, maybe an ObsID was provided.
        else:
            file_list[i] = load_obsid(_infile[i])
            
            #If we reached this point there is a problem with
            if len(file_list[i])==0:
                raise SwiftPhotom.errors.FileNotFound

    return file_list

def get_aperture_size(_reg):
    '''
    Retrieve the aperture size from the region file.
    '''
    with open(_reg) as inp:
        for line in inp:
            if line[0]=='#':continue
            size=line.strip('\n').split(',')[-1][:-2]
    
    size = size.rstrip('0')
    
    #Very brutal way to check if the size is an int
    if size[-1]=='.':
        size = size[:-1]
    
    return size

def check_aspect_correction(_infile):
    '''
    Simply check if the 'ASPCORR' label in the header
    of each extension is equal to 'DIRECT'. If it isn't,
    it prints out a WARNING, notifying that astrometry
    could be off.
    '''
    
    hdu=pf.open(_infile)
    for i in range(len(hdu)):
        if hdu[i].name=='PRIMARY': continue
        #print(hdu[i].header['FRAMTIME'])
        if hdu[i].header['ASPCORR'] !='DIRECT':
            ext = hdu[i].header['EXTNAME']
            print('WARNING - Extension '+str(i)+ ' '+ext+' of '+_infile+' has not been aspect corrected. You should check the astrometry.')
    hdu.close()
    del hdu

def sort_file_list(_flist):
    '''
    Organize the file list into a dictionary
    sorted by filter
    '''
    
    out_file_list={}
    for file in _flist:
        filter=pf.getheader(file)['FILTER']
        if filter not in out_file_list:
            out_file_list[filter]=[]
        out_file_list[filter].append(file)

    return out_file_list

def combine(_list,_outfile):
    '''
    Combines images with fappend
    '''
    
    for i,img in enumerate(_list):
        if i==0:
            sc.fcopy(img,_outfile)
        else:
            sc.fappend(img,_outfile)

def create_product(_flist,_filter,template=0,no_combine=0):
    '''
    If a file has multiple extensions, these are
    first combined into a single image. If the
    FRAMTIMEs are different, than no merging is
    not possible. The merging is carried out by
    uvotimsum, and the combined images are saved
    in the 'mid-products' folder, created for each
    filer inside the 'reduction' directory.
    
    If the user prefer to extract every extension
    separately, they can use the --no_combine flag.
    
    All the files are then appended to a single
    'product file', which end up containing all
    extensions of all images.
    
    During this step, the aspect correction for each
    file is also checked.
    '''
    
    out_dir=os.path.join('reduction',_filter,'mid-products')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    fig_dir=os.path.join('reduction',_filter,'figures')
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)

    prod_list=[]
    for file in _flist:
        check_aspect_correction(file)
        
        hdu=pf.open(file)
        
        if no_combine:
            for i in range(len(hdu)):
                if hdu[i].name=='PRIMARY': continue
                prod_list.append(file+'['+str(i)+']')
            continue
        
        obsID=hdu[0].header['OBS_ID']
        out_file=os.path.join(out_dir,obsID+'_'+_filter+'.fits')
        if os.path.isfile(out_file):
            os.remove(out_file)
        
        framtime=[]

        for i in range(len(hdu)):
            if hdu[i].name=='PRIMARY': continue
            framtime.append(hdu[i].header['FRAMTIME'])
        if len(set(framtime))==1:
            sc.uvotimsum(file,out_file,_exclude='none')
            prod_list.append(out_file)
        else:
            print('WARNING - extensions of '+file+' have different FRAMTIMEs. Left unmerged.')
            for i in range(len(hdu)):
                if hdu[i].name=='PRIMARY': continue
                prod_list.append(file+'['+str(i)+']')
        hdu.close()
        del hdu

    if not template:
        prod_out_file= os.path.join('reduction',_filter,pf.getheader(file)['OBJECT']+'_'+_filter+'.img')
    else:
        prod_out_file= os.path.join('reduction',_filter,'templ_'+_filter+'.img')
    if os.path.isfile(prod_out_file):
        os.remove(prod_out_file)

    combine(prod_list,prod_out_file)

    return prod_out_file

def run_uvotmaghist(_prod_file,_sn_reg,_bg_reg,_filter):
    fig_dir=os.path.join('reduction',_filter,'figures')
    photo_out=_prod_file[:-4]+'_phot.fits'
    gif_out=os.path.join(fig_dir,_prod_file.split('/')[-1][:-4]+'_phot.gif')
    if os.path.isfile(photo_out):
        os.remove(photo_out)
    if os.path.isfile(gif_out):
        os.remove(gif_out)
    sc.uvotmaghist(_prod_file,_sn_reg,_bg_reg, photo_out,gif_out)

    return photo_out

def extract_photometry(_phot_file, _ab, _det_limit, _ap_size, _templ_file=None):

    ####
    #The default aperture of uvotmaghist is 5 arcsec.
    #Corrections to retrieve this fluxes are provide.
    #We can compare the 2 apertures.
    ####
    
    if _templ_file!=None:
        template=1

    if _ab==1:
        Vega_corr={'V':0,'B':0,'U':0,'UVW1':0,'UVM2':0,'UVW2':0}
        mag_sys = 'AB'
    else:
        Vega_corr=Vega
        mag_sys = 'Vega'

    user_ap = _ap_size+'_arcsec'

    col={user_ap:'b','5_arcsec':'r'}
    mag={user_ap:[],'5_arcsec':[]}

    BCR_temp={}
    BCRe_temp={}

    for i,file in enumerate([_templ_file,_phot_file]):
        if file==None:
            #In case there is no template, nothing will be subtracted.
            BCR_temp={user_ap:0.,'5_arcsec':0.}
            BCRe_temp={user_ap:0.,'5_arcsec':0.}
            template=0
            continue
    
        hdu=pf.open(file)
        dd=hdu[1].data
        filter=dd['FILTER'][0]
        hdu.close()

        #The long-term detector sensitivity correction factor.
        SC=dd['SENSCORR_FACTOR']

        #This is the count rate for the user defined aperture
        S3BCR=dd['COI_SRC_RATE']*SC
        
        #adding 3% error on count rate of the source in quadrature to poission error
        if template and i==1:
            S3BCRe=np.sqrt((dd['COI_SRC_RATE_ERR'])**2+(S3BCR*0.03)**2)
        else:
            S3BCRe=dd['COI_SRC_RATE_ERR']
        
        #These is the count rate for 5arcsec
        S5CR=dd['RAW_STD_RATE']*dd['COI_STD_FACTOR']
        S5CRe=dd['RAW_STD_RATE_ERR']*dd['COI_STD_FACTOR']
        S5BCR=((dd['RAW_STD_RATE'] *dd['COI_STD_FACTOR'] * SC) - (dd['COI_BKG_RATE'] * SC * dd['STD_AREA']))
        if template and i==1:
            S5BCRe=np.sqrt((S5CRe)**2+(S5CR*0.03)**2+(dd['COI_BKG_RATE_ERR']*dd['STD_AREA'])**2)
        else:
            S5BCRe=np.sqrt((S5CRe)**2+(dd['COI_BKG_RATE_ERR']*dd['STD_AREA'])**2)
        
        
        
        fig_dir=os.path.join('reduction',filter,'figures')
        fig=plt.figure()
        ax=fig.add_subplot(111)
        
        if template:
            fig_mag=plt.figure()
            ax_mag=fig_mag.add_subplot(111)
        
        fig_sub=plt.figure()
        ax_sub=fig_sub.add_subplot(111)

        for BCR,BCRe,label in [[S3BCR,S3BCRe,user_ap] , [S5BCR,S5BCRe,'5_arcsec']]:
        



            if i==0:
                epochs=range(len(BCR))
                ax.errorbar(epochs,BCR,yerr=BCRe,marker='o', color=col[label],label=label)
            
                #weighed avarage the template fluxes
                BCR_temp[label]=np.sum(BCR/BCRe**2)/np.sum(1./BCRe**2)
                BCRe_temp[label]=np.sqrt(1./np.sum(1./BCRe**2))
                
                xx=[0,len(BCR)-1]
                ax.plot(xx,[BCR_temp[label]]*2,color=col[label],lw=2)
                ax.fill_between(xx, [BCR_temp[label]+BCRe_temp[label]]*2,  [BCR_temp[label]-BCRe_temp[label]]*2, color=col[label],lw=2, alpha=0.2)
            
                ax.set_xlabel('Epoch')
            
                print('Galaxy count rates in '+label.split('_')[0]+'" aperture: %.3f +- %.3f' %(BCR_temp[label],BCRe_temp[label]))

            else:
                all_point=[]
                mjd=mjdref+(dd['TSTART']+dd['TSTOP'])/2./86400.
                ax.errorbar(mjd,BCR,yerr=BCRe,marker='o', color=col[label],label=label)
                
                ax.set_xlabel('MJD')
                
                
                #subtract galaxy, propogate error
                BCGR=(BCR-BCR_temp[label])
                BCGRe=np.sqrt((BCRe)**2+(BCRe_temp[label])**2)

                if label==user_ap:
                    # apply aperture correction
                    BCGAR=BCGR*dd['AP_FACTOR']
                    BCGARe=BCGRe*dd['AP_FACTOR_ERR']
                    BCAR=BCR*dd['AP_FACTOR']
                    BCARe=BCRe*dd['AP_FACTOR_ERR']
                else:
                    BCGAR=BCGR
                    BCGARe=BCGRe
                    BCAR=BCR
                    BCARe=BCRe
                
                if template:
                    #Not subtracted magnitudes
                    or_mag=-2.5*np.log10(BCAR)+ZP[filter][0]-Vega_corr[filter]
                    or_mage=np.sqrt(((2.5/np.log(10.))*((BCRe/BCAR)))**2+ZP[filter][1]**2)
                
                    mag_host=-2.5*np.log10(BCR_temp[label])+ZP[filter][0]-Vega_corr[filter]
                    mag_hoste=np.sqrt(((2.5/np.log(10.))*((BCRe_temp[label]/BCR_temp[label])))**2+ZP[filter][1]**2)
                
                    ax_mag.errorbar(mjd,or_mag,yerr=or_mage,marker='o', color=col[label],label=label)

                    ax_mag.plot([min(mjd),max(mjd)],[mag_host]*2,color=col[label],lw=2)
                    ax_mag.fill_between([min(mjd),max(mjd)], [mag_host+mag_hoste]*2,  [mag_host-mag_hoste]*2, color=col[label],lw=2, alpha=0.2)

                
            



                # determine significance/3 sigma upper limit
                BCGARs=BCGAR/BCGARe
                BCGAMl=(-2.5*np.log10(3.*BCGARe))+ZP[filter][0]-Vega_corr[filter]

                #convert rate,err to magnitudes"
                for j,CR in enumerate(BCGARs):
                    mag[label].append({
                    'filter':filter,
                    'aperture_correction':float(dd['AP_FACTOR'][j]),
                    'coincidence_loss_correction':float(dd['COI_STD_FACTOR'][j]),
                    'mag_sys':mag_sys,
                    'mag_limit':float(BCGAMl[j]),
                    'template_subtracted':bool(template),
                    'mjd':float(mjd[j])
                    })
                    
                    #detection
                    if BCGARs[j]>=_det_limit:
                        BCGAM=-2.5*np.log10(BCGAR[j])+ZP[filter][0]-Vega_corr[filter]
                        BCGAMe=np.sqrt(((2.5/np.log(10.))*((BCGARe[j]/BCGAR[j])))**2+ZP[filter][1]**2)
                    
                        mag[label][-1]['upper_limit']=False
                        print('%.2f\t%.3f\t%.3f' % (mjd[j],BCGAM,BCGAMe))
                    #non detection
                    else:
                        BCGAM=BCGAMl[j]
                        BCGAMe=0.2
                        
                        mag[label][-1]['upper_limit']=True
                        print('%.2f\t> %.3f (%.2f)' % (mjd[j],BCGAM,np.fabs(BCGARs[j])))
                    
                    mag[label][-1]['mag'] = float(BCGAM)
                    mag[label][-1]['mag_err'] = float(BCGAMe)
                    all_point.append([mjd[j],BCGAM])
 
                detections = [[ep['mjd'], ep['mag'], ep['mag_err']] for ep in mag[label] if not ep['upper_limit']]
                non_detections = [[ep['mjd'], ep['mag'], ep['mag_err']] for ep in mag[label] if ep['upper_limit']]
                
                if len(detections)>0:
                    xx,yy,ee=zip(*sorted(detections))
                    ax_sub.errorbar(xx,yy,yerr=ee,marker='o', color=col[label],label=label,ls='')

                
                if len(non_detections)>0:
                    xx,yy,ee=zip(*sorted(non_detections))
                    ax_sub.errorbar(xx,yy, yerr=[-1.*np.array(ee),[0]*len(ee)], uplims=True, marker='o', color=col[label],label=label,ls='')
                
                #since I'm plotting detections and non detections separately
                #this is done only to have a line connecting all epochs
                xx,yy=zip(*sorted(all_point))
                ax_sub.plot(xx,yy, color=col[label])

                if label==user_ap:
                    print('\n'+_ap_size+'" aperture done!\n')
                else:
                    print('\n5" aperture done!\n')
        
        if i==1:
        
            ax_sub.set_title(file.split('/')[-1])
            ax_sub.set_xlabel('MJD')
            ax_sub.set_ylabel('Mag')
            ax_sub.invert_yaxis()
            ax_sub.legend()
            
            
            out_fig=os.path.join(fig_dir,file.split('/')[-1][:-5]+'_mag_final.png')
            if os.path.isfile(out_fig):
                os.remove(out_fig)

            fig_sub.savefig(out_fig)
            plt.close(fig_sub)
            del fig_sub
            del ax_sub
            
            
            if template:
                ax_mag.set_title(file.split('/')[-1])
                ax_mag.set_xlabel('MJD')
                ax_mag.set_ylabel('Mag')
                ax_mag.invert_yaxis()
                ax_mag.legend()


                out_fig=os.path.join(fig_dir,file.split('/')[-1][:-5]+'_mag.png')
                if os.path.isfile(out_fig):
                    os.remove(out_fig)

                fig_mag.savefig(out_fig)
                plt.close(fig_mag)
                del fig_mag
                del ax_mag

        ax.set_title(file.split('/')[-1])
        ax.set_ylabel('Coincident-corrected count rates')
        ax.legend()
        
        

        out_fig=os.path.join(fig_dir,file.split('/')[-1][:-5]+'_counts.png')
        if os.path.isfile(out_fig):
            os.remove(out_fig)

        fig.savefig(out_fig)
        plt.close(fig)
        del fig
        del ax
        


    return mag

def output_mags(_mag,_ap_size):
    user_ap = _ap_size+'_arcsec'

    with open(os.path.join('reduction',_ap_size+'_arcsec_photometry.json'),'w') as out:
        out.write(json.dumps(_mag[user_ap], indent = 4))
    
    with open(os.path.join('reduction','5_arcsec_photometry.json'),'w') as out:
        out.write(json.dumps(_mag['5_arcsec'], indent = 4))

    print('5-arcsec output photometry:\n')
    print('#'*80+'\n\n')
    _mag['5_arcsec'] = sorted(_mag['5_arcsec'], key=lambda x: x['mjd'])
    for photom in _mag['5_arcsec']:
        mjd = photom['mjd']
        filt = photom['filter']
        if photom['upper_limit']:
            mag = photom['mag_limit']
            magerr = 0.0
        else:
            mag = photom['mag']
            magerr = photom['mag_err']

        mjd = '%.5f'%mjd
        mag = '%.4f'%mag
        magerr = '%.4f'%magerr

        print(mjd,filt,mag,magerr)

