from awips.dataaccess import DataAccessLayer
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from math import exp, log
import numpy as np
from metpy.calc import wind_components, lcl, lfc, el, dry_lapse, parcel_profile, dewpoint, wet_bulb_temperature, mean_pressure_weighted
from metpy.calc import most_unstable_parcel, parcel_profile_with_lcl, bulk_shear, storm_relative_helicity, lifted_index
from metpy.calc import wind_speed, wind_direction, thermo, vapor_pressure, bunkers_storm_motion, pressure_to_height_std
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
import metpy.calc as mpcalc
from base64 import b64encode, b64decode
import codecs
import sys
from datetime import datetime
import datetime as dt
from metpy.units import units
import metpy.calc as mpcalc
from matplotlib.gridspec import GridSpec
from matplotlib import gridspec
import sharppy.sharptab.profile as profile
from sharppy.sharptab import utils, winds, params, interp, thermo, watch_type, fire
import sharppy.sharptab as tab
from sharppy.sharptab.constants import *
from sharppy.sharptab.constants import MISSING, TOL

# make unique directory to store output
def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

#grabbing data from NOMADS
startTime=datetime.now()

year = startTime.year

if startTime.month <10:
    month = '0'+str(startTime.month)
else:
    month = str(startTime.month)

if startTime.day <10:
    day = '0'+str(startTime.day)
else:
    day = str(startTime.day)

if startTime.hour <10:
    hour = '0'+str(startTime.hour)
else:
    hour = str(startTime.hour)
    

def get_init_hr(hour):
    if int(hour) <6:
        init_hour = '00'
    elif int(hour) <12:
        init_hour = '06'
    elif int(hour) <17:
        init_hour = '12'
    elif int(hour) <22:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)

mdate = str(year)+str(month)+str(day)
# Create new directory to store output
output_dir = str(year)+str(month)+str(day)  #this string names the output directory
mkdir_p(output_dir)
mkdir_p(output_dir+'/Soundings/GFS/CHA') #create subdirectory to store GFS output like this


DataAccessLayer.changeEDEXHost("edex-cloud.unidata.ucar.edu")
request = DataAccessLayer.newDataRequest("modelsounding")
forecastModel = "GFS"
request.addIdentifier("reportType", forecastModel)
request.setParameters("pressure","temperature","specHum","uComp","vComp","omega","cldCvr")

locations = DataAccessLayer.getAvailableLocationNames(request)
locations.sort()
list(locations)



for i in range(0,156):
    request.setLocationNames("KCHA")
    cycles = DataAccessLayer.getAvailableTimes(request, True)
    times = DataAccessLayer.getAvailableTimes(request)

    try:
        fcstRun = DataAccessLayer.getForecastRun(cycles[-1], times)
        list(fcstRun)
        response = DataAccessLayer.getGeometryData(request,[fcstRun[i]])
    except:
        print('No times available')
        exit

    tmp,prs,sh = np.array([]),np.array([]),np.array([])
    uc,vc,om,cld = np.array([]),np.array([]),np.array([]),np.array([])


    for ob in response:
        tmp = np.append(tmp,ob.getNumber("temperature"))
        prs = np.append(prs,ob.getNumber("pressure"))
        sh = np.append(sh,ob.getNumber("specHum"))
        uc = np.append(uc,ob.getNumber("uComp"))
        vc = np.append(vc,ob.getNumber("vComp"))
        om = np.append(om,ob.getNumber("omega"))
        cld = np.append(cld,ob.getNumber("cldCvr"))


    print("parms    = " + str(ob.getParameters()))
    print("site     = " + str(ob.getLocationName()))
    print("geom     = " + str(ob.getGeometry()))
    print("datetime = " + str(ob.getDataTime()))
    print("reftime  = " + str(ob.getDataTime().getRefTime()))
    print("fcstHour = " + str(ob.getDataTime().getFcstTime()))
    print("period   = " + str(ob.getDataTime().getValidPeriod()))
    fcstHour = str(ob.getDataTime().getFcstTime())
    t = (tmp-273.15) * units.degC
    p = prs/100 * units.mbar
    z = mpcalc.pressure_to_height_std(p)
    print(t)
    print(p)
    print(z)
    sfc_hgt = z[0]
    print(sfc_hgt)

    u,v = uc*1.94384,vc*1.94384 # m/s to knots
    spd = wind_speed(u*units.knots, v*units.knots)
    dir = wind_direction(u*units.knots, v*units.knots) * units.deg

    rmix = (sh/(1-sh)) *1000 * units('g/kg')
    e = vapor_pressure(p, rmix)
    td = dewpoint(e)
    
    #td = metpy.calc.dewpoint_from_relative_humidity(t, rh)

    td2 = dewpoint(vapor_pressure(p, sh))
    
    # Calculate Wetbulb Temperature - need to figure out what to do with nan values
    #wetbulb = wet_bulb_temperature(p, t, td)
    # Sets up the prof object - not sure if this really will work
    '''
    p_sounding = np.sort(np.append(lev, sfc[0,ilat[0][0], ilon[1][0]]))
    ind = np.where(p_sounding >= sfcp[0,ilat[0][0], ilon[1][0]])[0][0]
    hgt_sounding = np.insert(z[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, sfc_hgt[0,ilat[0][0], ilon[1][0]].magnitude) * hgt.units
    T_sounding = (np.insert(t[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, t[0,0,ilat[0][0], ilon[1][0]].magnitude) * t.units).to(tdc.units)
    Td_sounding = np.insert(tdc_up.magnitude, ind, td[0,0,ilat[0][0], ilon[1][0]].magnitude) * tdc_up.units
    u_sounding = np.insert(uc[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, u10[0,0,ilat[0][0], ilon[1][0]].magnitude) * usfc.units
    v_sounding = np.insert(vc[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, v10[0,0,ilat[0][0], ilon[1][0]].magnitude) * usfc.units

    p_skewt = p_sounding[p_sounding <= sfcp[0,ilat[0][0], ilon[1][0]]]
    hgt_skewt = hgt_sounding[p_sounding <= sfcp[0,ilat[0][0], ilon[1][0]]]
    T_skewt = T_sounding[p_sounding <= sfcp[0,ilat[0][0], ilon[1][0]]]
    Td_skewt = Td_sounding[p_sounding <= sfcp[0,ilat[0][0], ilon[1][0]]]
    u_skewt = u_sounding[p_sounding <= sfcp[0,ilat[0][0], ilon[1][0]]].to('kt')
    v_skewt = v_sounding[p_sounding <= sfcp[0,ilat[0][0], ilon[1][0]]].to('kt')
    
    prof = profile.create_profile(profile='default', pres=p, hght=z, tmpc=t, dwpc=td2, wspd=spd, wdir=dir, missing=-9999, strictQC=False)
    prof = profile.create_profile(profile='default', pres=p_skewt[::-1], hght=hgt_skewt[::-1], tmpc=T_skewt[::-1], dwpc=Td_skewt[::-1], wspd=wind_spd[::-1], wdir=wind_dir[::-1], missing=-9999, strictQC=False)                                    

    mupcl = params.parcelx(prof,flag=3)
    '''
    # Ends this block of code
    def __mu(self, prof, **kwargs):
        '''
            Create the most unstable parcel within the lowest XXX hPa, where
            XXX is supplied. Default XXX is 400 hPa.
            
            '''
        self.desc = 'Most Unstable Parcel in Lowest %.2f hPa' % self.presval
        pbot = prof.pres[prof.sfc]
        ptop = pbot - self.presval
        self.pres = most_unstable_level(prof, pbot=pbot, ptop=ptop)
        self.tmpc = interp.temp(prof, self.pres)
        self.dwpc = interp.dwpt(prof, self.pres)
    
    
    def __ml(self, prof, **kwargs):
        '''
            Create a mixed-layer parcel with mixing within the lowest XXX hPa,
            where XXX is supplied. Default is 100 hPa.
            If
            
            '''
        self.desc = '%.2f hPa Mixed Layer Parcel' % self.presval
        pbot = kwargs.get('pbot', prof.pres[prof.sfc])
        ptop = pbot - self.presval
        self.pres = pbot
        mtheta = mean_theta(prof, pbot, ptop, exact=True)
        self.tmpc = thermo.theta(1000., mtheta, self.pres)
        mmr = mean_mixratio(prof, pbot, ptop, exact=True)
        self.dwpc = thermo.temp_at_mixrat(mmr, self.pres)
    
    
    def __user(self, prof, **kwargs):
        '''
            Create a user defined parcel.
            
            '''
        self.desc = '%.2f hPa Parcel' % self.presval
        self.pres = self.presval
        self.tmpc = kwargs.get('tmpc', interp.temp(prof, self.pres))
        self.dwpc = kwargs.get('dwpc', interp.dwpt(prof, self.pres))
    
    
    def __effective(self, prof, **kwargs):
        '''
            Create the mean-effective layer parcel.
            
            '''
        ecape = kwargs.get('ecape', 100)
        ecinh = kwargs.get('ecinh', -250)
        pbot, ptop = effective_inflow_layer(prof, ecape, ecinh)
        if utils.QC(pbot) and pbot > 0:
            self.desc = '%.2f hPa Mean Effective Layer Centered at %.2f' % ( pbot-ptop, (pbot+ptop)/2.)
            mtha = mean_theta(prof, pbot, ptop)
            mmr = mean_mixratio(prof, pbot, ptop)
            self.pres = (pbot+ptop)/2.
            self.tmpc = thermo.theta(1000., mtha, self.pres)
            self.dwpc = thermo.temp_at_mixrat(mmr, self.pres)
        else:
            self.desc = 'Defaulting to Surface Layer'
            self.pres = prof.pres[prof.sfc]
            self.tmpc = prof.tmpc[prof.sfc]
            self.dwpc = prof.dwpc[prof.sfc]
        if utils.QC(pbot): self.pbot = pbot
        else: self.pbot = ma.masked
        if utils.QC(ptop): self.ptop = ptop
        else: self.pbot = ma.masked


    def helicity(prof, lower, upper, stu=0, stv=0, dp=-1, exact=True):
        '''
        Calculates the relative helicity (m2/s2) of a layer from lower to upper.
        If storm-motion vector is supplied, storm-relative helicity, both
        positve and negative, is returned.
        Parameters
        ----------
        prof : profile object
            Profile Object
        lower : number
            Bottom level of layer (m, AGL)
        upper : number
            Top level of layer (m, AGL)
        stu : number (optional; default = 0)
            U-component of storm-motion (kts)
        stv : number (optional; default = 0)
            V-component of storm-motion (kts)
        dp : negative integer (optional; default -1)
            The pressure increment for the interpolated sounding (mb)
        exact : bool (optional; default = True)
            Switch to choose between using the exact data (slower) or using
            interpolated sounding at 'dp' pressure levels (faster)
        Returns
        -------
        phel+nhel : number
            Combined Helicity (m2/s2)
        phel : number
            Positive Helicity (m2/s2)
        nhel : number
            Negative Helicity (m2/s2)
        '''
        if prof.wdir.count() == 0 or not utils.QC(lower) or not utils.QC(upper) or not utils.QC(stu) or not utils.QC(stv):
            return ma.masked, ma.masked, ma.masked

        if lower != upper:
            lower = interp.to_msl(prof, lower)
            upper = interp.to_msl(prof, upper)
            plower = interp.pres(prof, lower)
            pupper = interp.pres(prof, upper)
            if np.isnan(plower) or np.isnan(pupper) or \
                type(plower) == type(ma.masked) or type(pupper) == type(ma.masked):
                return np.ma.masked, np.ma.masked, np.ma.masked
            if exact:
                ind1 = np.where(plower >= prof.pres)[0].min()
                ind2 = np.where(pupper <= prof.pres)[0].max()
                u1, v1 = interp.components(prof, plower)
                u2, v2 = interp.components(prof, pupper)
                u = np.concatenate([[u1], prof.u[ind1:ind2+1].compressed(), [u2]])
                v = np.concatenate([[v1], prof.v[ind1:ind2+1].compressed(), [v2]])
            else:
                ps = np.arange(plower, pupper+dp, dp)
                u, v = interp.components(prof, ps)
            sru = utils.KTS2MS(u - stu)
            srv = utils.KTS2MS(v - stv)
            layers = (sru[1:] * srv[:-1]) - (sru[:-1] * srv[1:])
            phel = layers[layers > 0].sum()
            nhel = layers[layers < 0].sum()
        else:
            phel = nhel = 0

        return phel+nhel, phel, nhel
        
    

    def sr_wind(prof, pbot=850, ptop=250, stu=0, stv=0, dp=-1):
        '''
        Calculates a pressure-weighted mean storm-relative wind through a layer.
        The default layer is 850 to 200 hPa. This is a thin wrapper around
        mean_wind().
        Parameters
        ----------
        prof: profile object
            Profile object
        pbot : number (optional; default 850 hPa)
            Pressure of the bottom level (hPa)
        ptop : number (optional; default 250 hPa)
            Pressure of the top level (hPa)
        stu : number (optional; default 0)
            U-component of storm-motion vector (kts)
        stv : number (optional; default 0)
            V-component of storm-motion vector  (kts)
        dp : negative integer (optional; default -1)
            The pressure increment for the interpolated sounding (mb)
        Returns
        -------
        mnu : number
            U-component (kts)
        mnv : number
            V-component (kts)
        '''
        return mean_wind(prof, pbot=pbot, ptop=ptop, dp=dp, stu=stu, stv=stv)
        
    def get_indices(self):
        '''
        Function to set any additional indices that are included in the 
        thermo window.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        '''
        self.tei = params.tei(self)
        self.esp = params.esp(self)
        self.mmp = params.mmp(self)
        self.wndg = params.wndg(self)
        self.sig_severe = params.sig_severe(self)
        self.dcape, self.dpcl_ttrace, self.dpcl_ptrace = params.dcape(self)
        self.drush = thermo.ctof(self.dpcl_ttrace[-1])
        self.mburst = params.mburst(self)
    
    # Calculate LCL height and plot as black dot
    lcl_pressure, lcl_temperature = lcl(p[0], t[0], td[0])
    
    # Calculate LFC height and plot as yellow dash
    lfc_pressure, lfc_temperature = lfc(p, t, td)
    
    el_pressure, el_temperature = el(p, t, td)
    
    lcl_hgt = np.round(mpcalc.pressure_to_height_std(lcl_pressure), decimals=1).to(units.meter)/1000
    lfc_hgt = np.round(mpcalc.pressure_to_height_std(lfc_pressure), decimals=1).to(units.meter)/1000
    el_hgt = np.round(mpcalc.pressure_to_height_std(el_pressure), decimals=1).to(units.meter)/1000
    
    sb_cape, sb_cin = mpcalc.surface_based_cape_cin(p, t, td)
    ml_cape, ml_cin = mpcalc.mixed_layer_cape_cin(p, t, td)
    mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(p, t, td)
    
    #muparc = mpcalc.most_unstable_parcel(p, t, td)
    #print(muparc)
    #muparc_pressure = np.round(mpcalc.height_to_pressure_std(muparc), decimals=3).to(units.meter)
    #print(muparc_pressure)
    
    mixed_0_3 = mpcalc.mixed_parcel(p, t, td, depth=3000 * units.meter)
    print(mixed_0_3)    

    
    sbcape = np.round(sb_cape, 1)
    sbcin = np.round(sb_cin, 1)
    mlcape = np.round(ml_cape, 1)
    mlcin = np.round(ml_cin, 1)
    mucape = np.round(mu_cape, 1)
    pw = mpcalc.precipitable_water(p, td)
    pw = pw.to(units.inch)
    pw = round(pw, 2)
    
    u_shear01, v_shear01 = mpcalc.bulk_shear(p, u * units('m/s'), v * units('m/s'), depth = 1000 * units.meter)
    shear01 = np.round((np.sqrt(u_shear01**2 + v_shear01**2)), 1)
    shear01 = shear01.to(units.knots)
    shear01 = np.round(shear01)
    u_shear005, v_shear005 = mpcalc.bulk_shear(p, u * units('m/s'), v * units('m/s'), depth = 500 * units.meter)
    shear005 = np.round((np.sqrt(u_shear005**2 + v_shear005**2)),1)
    shear005 = shear005.to(units.knots)
    shear005 = np.round(shear005)
    u_shear015, v_shear015 = mpcalc.bulk_shear(p, u * units('m/s'), v * units('m/s'), depth = 1500 * units.meter)
    shear015 = np.round((np.sqrt(u_shear015**2 + v_shear015**2)),1)
    shear015 = shear015.to(units.knots)
    shear015 = np.round(shear015)
    u_shear02, v_shear02 = mpcalc.bulk_shear(p, u * units('m/s'), v * units('m/s'), depth = 2000 * units.meter)
    shear02 = np.round((np.sqrt(u_shear02**2 + v_shear02**2)), 1)
    shear02 = shear02.to(units.knots)
    shear02 = np.round(shear02)
    u_shear03, v_shear03 = mpcalc.bulk_shear(p, u * units('m/s'), v * units('m/s'), depth = 3000 * units.meter)
    shear03 = np.round((np.sqrt(u_shear03**2 + v_shear03**2)), 1)
    shear03 = shear03.to(units.knots)
    shear03 = np.round(shear03)
    u_shear06, v_shear06 = mpcalc.bulk_shear(p,  u * units('m/s'), v * units('m/s'), depth = 6000 * units.meter)
    shear06 = np.round((np.sqrt(u_shear06**2 + v_shear06**2)), 1)
    shear06 = shear06.to(units.knots)
    shear06 = np.round(shear06)
    rmover, lmover, mean = mpcalc.bunkers_storm_motion(p, u * units('m/s'), v * units('m/s'), z)
    #srh_01_pos, srh_01_neg, srh_01_tot = mpcalc.storm_relative_helicity( u * units('m/s'), v * units('m/s'), z * units('m'), depth = 1000 * units.meter, bottom = z[0], storm_u = lmover[0], storm_v = lmover[1])
    #srh_01 = np.round(srh_01_neg, 1)
    #srh_03_pos, srh_03_neg, srh_03_tot = mpcalc.storm_relative_helicity( u * units('m/s'), v * units('m/s'), z * units('m'), depth = 3000 * units.meter, bottom = z[0], storm_u = lmover[0], storm_v = lmover[1])
    #srh_03 = np.round(srh_03_neg, 1)
    '''
    # Compute the low-level (SFC-1500m) mean wind
    p_1p5km = interp.pres(prof, interp.to_msl(prof, 1500.))
    mnu2, mnv2 = mean_wind_npw(prof, prof.pres[prof.sfc], p_1p5km)

    # Compute the upshear vector
    upu = mnu1 - mnu2
    upv = mnv1 - mnv2

    # Compute the downshear vector
    dnu = mnu1 + upu
    dnv = mnv1 + upv

    #return upu, upv, dnu, dnv
    
    
    ## K Index
    self.k_idx = tab.utils.INT2STR( prof.k_idx )
    ## precipitable water
    self.pwat = prof.pwat
    ## 0-3km agl lapse rate
    self.lapserate_3km = tab.utils.FLOAT2STR( prof.lapserate_3km, 1 )
    ## 3-6km agl lapse rate
    self.lapserate_3_6km = tab.utils.FLOAT2STR( prof.lapserate_3_6km, 1 )
    ## 850-500mb lapse rate
    self.lapserate_850_500 = tab.utils.FLOAT2STR( prof.lapserate_850_500, 1 )
    ## 700-500mb lapse rate
    self.lapserate_700_500 = tab.utils.FLOAT2STR( prof.lapserate_700_500, 1 )
    ## convective temperature
    self.convT = prof.convT
    ## sounding forecast surface temperature
    self.maxT = prof.maxT
    #fzl = str(int(self.sfcparcel.hght0c))
    ## 100mb mean mixing ratio
    self.mean_mixr = tab.utils.FLOAT2STR( prof.mean_mixr, 1 )
    ## 150mb mean rh
    self.low_rh = tab.utils.INT2STR( prof.low_rh )
    self.mid_rh = tab.utils.INT2STR( prof.mid_rh )
    ## calculate the totals totals index
    self.totals_totals = tab.utils.INT2STR( prof.totals_totals )
    self.dcape = tab.utils.INT2STR( prof.dcape )
    self.drush = prof.drush
    self.sigsevere = tab.utils.INT2STR( prof.sig_severe )
    self.mmp = tab.utils.FLOAT2STR( prof.mmp, 2 )
    self.esp = tab.utils.FLOAT2STR( prof.esp, 1 )
    self.wndg = tab.utils.FLOAT2STR( prof.wndg, 1 )
    self.tei = tab.utils.INT2STR( prof.tei )
    '''
    
    right_mover,left_mover,wind_mean = mpcalc.bunkers_storm_motion(p, u * units('m/s'), v * units('m/s'), z)
    wind_mean = np.round(wind_mean)
    wind_mean = wind_mean.to(units.kt)
    wind_mean = np.round(wind_mean)

    pos_SRH,neg_SRH,total_SRH = mpcalc.storm_relative_helicity(z, u * units('m/s'), v * units('m/s'), depth = 3000 * units('m'), bottom = z[0], storm_u = lmover[0], storm_v = lmover[1])
    pos1_SRH,neg1_SRH,total1_SRH = mpcalc.storm_relative_helicity(z, u * units('m/s'), v * units('m/s'), depth = 1000 * units('m'), bottom = z[0], storm_u = lmover[0], storm_v = lmover[1])
    pos2_SRH,neg2_SRH,total2_SRH = mpcalc.storm_relative_helicity(z, u * units('m/s'), v * units('m/s'), depth = 2000 * units('m'), bottom = z[0], storm_u = lmover[0], storm_v = lmover[1])
    pos05_SRH,neg05_SRH,total05_SRH = mpcalc.storm_relative_helicity(z, u * units('m/s'), v * units('m/s'), depth = 500 * units('m'), bottom = z[0], storm_u = lmover[0], storm_v = lmover[1])
    # Not sure if this is really starting from 1km for the  bottom - need to investigate this more
    pos13_SRH,neg13_SRH,total13_SRH = mpcalc.storm_relative_helicity(z, u * units('m/s'), v * units('m/s'), depth = 3000 * units('m'), bottom = z[10], storm_u = lmover[0], storm_v = lmover[1])
    #print(pos_SRH)
    #print(neg_SRH)
    #print(total_SRH)
    # new arrays
    # Need to round these numbers to make the string look pretty
    tot_SRH = np.round(total_SRH)
    tot1_SRH = np.round(total1_SRH)
    tot05_SRH = np.round(total05_SRH)
    tot13_SRH = np.round(total13_SRH)
    tot2_SRH = np.round(total2_SRH)
    
    # Layer Lapse Rates
    lr_05 = concatenate((mean_pressure_weighted(p, t, height=z, depth=500 * units('meter'))))
    lr_13 = concatenate((mean_pressure_weighted(p, t, height=z, depth=2000 * units('meter'), bottom=z[0] + 1000 * units('meter'))))
    lr_36 = concatenate((mean_pressure_weighted(p, t, height=z, depth=3000 * units('meter'), bottom=z[0] + 3000 * units('meter'))))
    lr_05 = np.round(lr_05)
    lr_13 = np.round(lr_13)
    lr_36 = np.round(lr_36)
    
    wind_mean6 = concatenate((mean_pressure_weighted(p, u * units('m/s'), v * units('m/s'), height=z, depth=6000 * units('meter'))))
    # mean wind from sfc-500m
    wind_500m = concatenate(mean_pressure_weighted(p, u * units('m/s'), v * units('m/s'), height=z,depth=500 * units('meter')))
    # mean wind from 5.5-6km
    wind_5500m = concatenate(mean_pressure_weighted(p, u * units('m/s'), v * units('m/s'), height=z, depth=500 * units('meter'), bottom=z[0] + 5500 * units('meter')))
    # mean wind from sfc-1km
    wind_mean1 = concatenate(mean_pressure_weighted(p, u * units('m/s'), v * units('m/s'), height=z,depth=1000 * units('meter')))
    # mean wind from 6-9km
    wind_6_9m = concatenate(mean_pressure_weighted(p, u * units('m/s'), v * units('m/s'), height=z, depth=3000 * units('meter'), bottom=z[0] + 6000 * units('meter')))
    # mean wind from 8.5-9km
    wind_8500m = concatenate(mean_pressure_weighted(p, u * units('m/s'), v * units('m/s'), height=z, depth=500 * units('meter'), bottom=z[0] + 8500 * units('meter')))
    
    shear = wind_mean1 - wind_500m
    shear_cross = concatenate([shear[1], -shear[0]])
    shear_mag = np.hypot(*(arg.magnitude for arg in shear)) * shear.units
    rdev = shear_cross * (7.5 / shear_mag)
    
    # Add the deviations to the layer average wind to get the RM motion
    right_mover1 = wind_mean1 + rdev * units('m/s')
    right_mover1 = right_mover1.to(units.knots)
    right_mover1 = np.round(right_mover1)

    # Subtract the deviations to get the LM motion
    left_mover1 = wind_mean1 - rdev * units('m/s')
    left_mover1 = left_mover1.to(units.knots)
    left_mover1 = np.round(left_mover1)
    
    shear69 = wind_8500m - wind_mean6
    shear_cross69 = concatenate([shear69[1], -shear69[0]])
    shear_mag69 = np.hypot(*(arg.magnitude for arg in shear)) * shear.units
    rdev69 = shear_cross69 * (7.5 / shear_mag)
    
    right_mover69 = wind_mean6 + rdev69 * units('m/s')
    right_mover69 = right_mover69.to(units.knots)
    right_mover69 = np.round(right_mover69)
    
    left_mover69 = wind_mean6 - rdev69 * units('m/s')
    left_mover69 = left_mover69.to(units.knots)
    left_mover69 = np.round(left_mover69)

    # Need to first convert to knots then round to nearest whole number
    wind_mean6 = wind_mean6.to(units.knots)
    wind_mean6 = np.round(wind_mean6)
    wind_500m = wind_500m.to(units.knots)
    wind_500m = np.round(wind_500m)
    wind_5500m = wind_5500m.to(units.knots)
    wind_5500m = np.round(wind_5500m)
    wind_mean1 = wind_mean1.to(units.knots)
    wind_mean1 = np.round(wind_mean1)
    wind_6_9m = wind_6_9m.to(units.knots)
    wind_6_9m = np.round(wind_6_9m)
    
    
    bunk_right_dir = np.round(mpcalc.wind_direction(right_mover[0], right_mover[1]))
    bunk_left_dir = np.round(mpcalc.wind_direction(left_mover[0], left_mover[1]))
    bunk_right_spd = np.round(np.sqrt(right_mover[0]**2 + right_mover[1]**2))
    bunk_right_spd = bunk_right_spd.to(units.knots)
    bunk_right_spd = np.round(bunk_right_spd)
    bunk_left_spd = np.round(np.sqrt(left_mover[0]**2 + left_mover[1]**2))
    bunk_left_spd = bunk_left_spd.to(units.knots)
    bunk_left_spd = np.round(bunk_left_spd)
    bunk_right_dir = np.round(bunk_right_dir)
    bunk_left_dir = np.round(bunk_left_dir)
    #bunk_right_spd = bunk_right_spd * units.knots
    #bunk_left_spd = bunk_left_spd * units.knots
    
    # Calculate composite parameters
    ehi_01 = np.round(np.divide(neg1_SRH * sbcape, 160000 * ((units.m**2 * units.joule)/(units.s**2 * units.kilogram))), 1)
    ehi_03 = np.round(np.divide(neg_SRH * sbcape, 160000 * ((units.m**2 * units.joule)/(units.s**2 * units.kilogram))), 1)
    scp = np.round(np.divide(sbcape, 1000 * units('J/kg')) * np.divide(shear06, 20 * units('m/s')) * np.divide(neg_SRH, 100 * (units.m**2/units.s**2)), 1)
    sig_tor =  np.round((mpcalc.significant_tornado(sbcape, lcl_hgt, neg1_SRH, shear06)), 1)

    # Calculate critical angle, then round to nearest whole number
    critical_angle = mpcalc.critical_angle(p, u * units('m/s'), v * units('m/s'), z, u_storm = lmover[0], v_storm = lmover[1])
    ca = np.round(critical_angle)
    
    ntmp = tmp

    # where p=pressure(pa), T=temp(C), T0=reference temp(273.16)
    rh = 0.263*prs*sh / (np.exp(17.67*ntmp/(ntmp+273.15-29.65)))
    vaps =  6.112 * np.exp((17.67 * ntmp) / (ntmp + 243.5))
    vapr = rh * vaps / 100
    dwpc = np.array(243.5 * (np.log(6.112) - np.log(vapr)) / (np.log(vapr) - np.log(6.112) - 17.67)) * units.degC

    plt.rcParams['figure.figsize'] = (12, 14)
    fig = plt.figure(figsize=(24, 14))

    gs = fig.add_gridspec(ncols=2,nrows=1)

	# identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
    #ax3 = fig.add_subplot(gs[1, 0])
    #ax4 = fig.add_subplot(gs[1, 1])

    # Grid for plots
    #skew = SkewT(fig, rotation=45, subplot=gs[0, 0])
    
    # Grid for plots
    skew = SkewT(fig, rotation=45, subplot=gs[0, 0])



    # Plot the data
    skew.plot(p, t, 'r', linewidth=2)
    skew.plot(p, td, 'b', linewidth=2)
    skew.plot(p, td2, 'y')
    skew.plot(p, dwpc, 'g', linewidth=2)
    skew.plot_moist_adiabats(color='grey',alpha=0.2)
    skew.plot_mixing_lines(color='grey',alpha=0.2)
    skew.plot_dry_adiabats(color='grey',alpha=0.2)
    #skew.plot(p, wetbulb, 'b', linewidth=1)
    skew.plot(lcl_pressure, lcl_temperature, marker="_", color='orange', markersize=30, markeredgewidth=3, label='LCL')
    skew.plot(lfc_pressure, lfc_temperature, marker="_", color='brown', markersize=30, markeredgewidth=2, label='LFC')
    skew.plot(el_pressure, el_temperature, marker="_", color='darkblue', markersize=30, markeredgewidth=2, label='EL')
    skew.ax.text(0.90, lfc_pressure, '-   LFC', verticalalignment='center', color='orange', alpha=0.9)
    skew.ax.text(0.90, lcl_pressure, '-   LCL', verticalalignment='center', color='brown', alpha=0.9)
    skew.ax.text(0.90, el_pressure, '-   EL', verticalalignment='center', color='darkblue', alpha=0.9)
    # Calculate full parcel profile and add to plot as black line
    prof = parcel_profile(p, t[0], td[0]).to('degC')
    prof_mu = parcel_profile(p, t[0], td[0]).to('degC')
    prof_lcl = parcel_profile_with_lcl(p, t, td)
    prof_li = lifted_index(p, t, prof)
    #li = mpcalc.lifted_index(p, t[0], prof)
    print(prof_li)
    
    lr_700_500 = np.round(-1 * np.divide(t[24]-t[14], (z[24]-z[14])),2)
    lr_850_500 = np.round(-1 * np.divide(t[24]-t[12], (z[24]-z[12])),2)
    lr_sfc_3 = np.round(-1 * np.divide(t[18]-t[0], (z[18]-z[0])),2)
    
    # Kinematic Calculations
    bulkshear = bulk_shear(p, u*units.knots, v*units.knots)
    
    skew.plot(p, prof, 'k', linewidth=2)
    #skew.plot(p, prof_mu, 'red', linewidth=2)
    #skew.plot(p, prof_lcl, 'red', linewidth=2)
    
    # Shade areas of CAPE and CIN
    #skew.shade_cin(p, t, prof, td)
    skew.shade_cape(p, t, prof, alpha=0.3)
    skew.shade_cape(p, t, prof_mu, alpha=0.5)

    skew.plot_barbs(p, u, v)
    skew.ax.set_ylim(1000, 100)
    plt.xlabel("Temperature [C]")
    skew.ax.set_xlim(-40, 60)
    plt.ylabel("Height [m above MSL]")
    
    plt.figtext( 0.15, 0.36, 'Levels (km):', fontsize=12)
    plt.figtext( 0.15, 0.35, 'LCL:', fontsize=12)
    plt.figtext( 0.17, 0.35, f'{lcl_hgt}', fontsize=12)
    plt.figtext( 0.15, 0.34, 'LFC:', fontsize=12)
    plt.figtext( 0.17, 0.34, f'{lfc_hgt:~P}', fontsize=12)
    plt.figtext( 0.15, 0.33, 'EL:', fontsize=12)
    plt.figtext( 0.17, 0.33, f'{el_hgt:~P}', fontsize=12)
    plt.figtext( 0.28, 0.20, 'MLLR:', fontsize=14)
    plt.figtext( 0.31, 0.20, f'{lr_700_500:~P}', fontsize=14)
    plt.figtext( 0.28, 0.18, 'LLLR:', fontsize=14)
    plt.figtext( 0.31, 0.18, f'{lr_sfc_3:~P}', fontsize=14)
    plt.figtext( 0.28, 0.16, '0-0.5:', fontsize=14)
    plt.figtext( 0.33, 0.16, f'{lr_05[0]:~P}', fontsize=14)
    plt.figtext( 0.28, 0.14, '1-3:', fontsize=14)
    plt.figtext( 0.33, 0.14, f'{lr_13[0]:~P}', fontsize=14)
    plt.figtext( 0.28, 0.12, '3-6:', fontsize=14)
    plt.figtext( 0.33, 0.12, f'{lr_36[0]:~P}', fontsize=14)
    plt.figtext( 0.37, 0.20, 'SBCAPE:', fontsize=14)
    plt.figtext( 0.41, 0.20, f'{sbcape:~P}', fontsize=14)
    plt.figtext( 0.37, 0.18, 'SBCIN:', fontsize=14)
    
    plt.figtext( 0.41, 0.18, f'{sbcin:~P}', fontsize=14)
    plt.figtext( 0.37, 0.16, 'MLCAPE:', fontsize=14)
    plt.figtext( 0.41, 0.16, f'{mlcape:~P}', fontsize=14)
    plt.figtext( 0.37, 0.14, 'MLCIN:', fontsize=14)
    plt.figtext( 0.41, 0.14, f'{mlcin:~P}', fontsize=14)
    plt.figtext( 0.37, 0.12, 'MUCAPE:', fontsize=14)
    plt.figtext( 0.41, 0.12, f'{mucape:~P}', fontsize=14)
    plt.figtext( 0.45, 0.20, '0-1km EHI:', fontsize=14)
    plt.figtext( 0.50, 0.20, f'{ehi_01:~P}', fontsize=14)
    plt.figtext( 0.45, 0.18, '0-3km EHI:', fontsize=14)
    plt.figtext( 0.50, 0.18, f'{ehi_03:~P}', fontsize=14)
    plt.figtext( 0.45, 0.16, 'SCP:', fontsize=14)
    plt.figtext( 0.49, 0.16, f'{scp:~P}', fontsize=14)
    plt.figtext( 0.45, 0.14, 'STP:', fontsize=14)
    plt.figtext( 0.49, 0.14, f'{sig_tor[0]:~P}', fontsize=14)
    plt.figtext( 0.45, 0.12, 'EL T:', fontsize=14)
    plt.figtext( 0.49, 0.12, f'{el_temperature:.2~P}', fontsize=14)
    #plt.figtext( 0.18, 0.25, f'{el_pressure:.2~P}')
    #plt.figtext( 0.80, 0.42, 'BulkShear:')
    #plt.figtext( 0.80, 0.42, '{0} knots'.format(bulkshear))
    #plt.figtext( 0.25, 0.32, 'LI:')
    #plt.figtext( 0.27, 0.32, f'{li:~P}')
    #plt.figtext( 0.25, 0.31, 'H7-H5 LR:')
    #plt.figtext( 0.32, 0.31, f'{lr_700_500:~P}')
    #plt.figtext( 0.25, 0.30, 'H8-H5 LR:')
    #plt.figtext( 0.32, 0.30, f'{lr_850_500:~P}')
    #plt.figtext( 0.25, 0.29, 'Sfc-3km LR:')
    #plt.figtext( 0.32, 0.29, f'{lr_sfc_3:~P}')
    
    plt.figtext( 0.79, 0.78, 'Bulk Shear', fontsize=14)
    plt.figtext( 0.79, 0.76, '0-0.5km:' ,fontsize=14)
    plt.figtext( 0.84, 0.76, f'{shear005:~P}', fontsize=14, color='purple')
    plt.figtext( 0.79, 0.74, '0-1 km:', fontsize=14)
    plt.figtext( 0.84, 0.74, f'{shear01:~P}', fontsize=14, color='purple')
    plt.figtext( 0.79, 0.72, '0-1.5km:',fontsize=14)
    plt.figtext( 0.84, 0.72, f'{shear015:~P}', fontsize=14, color='purple')
    plt.figtext( 0.79, 0.70, '0-2 km:', fontsize=14)
    plt.figtext( 0.84, 0.70, f'{shear02:~P}', fontsize=14)
    plt.figtext( 0.79, 0.68, '0-3 km:', fontsize=14)
    plt.figtext( 0.84, 0.68, f'{shear03:~P}', fontsize=14)
    plt.figtext( 0.79, 0.66, '0-6 km:', fontsize=14)
    plt.figtext( 0.84, 0.66, f'{shear06:~P}', fontsize=14)
    plt.figtext( 0.79, 0.64, 'SR Wind:', fontsize=14)
    #plt.figtext( 0.84, 0.59, '{}'.format(sr_wind))
    #plt.figtext( 0.65, 0.30, 'SRH 0-1 km:')
    #plt.figtext( 0.8, 0.30, f'{srh_01:~P}')
    plt.figtext ( 0.55, 0.78, 'Critical SRH', fontsize=14)
    plt.figtext ( 0.55, 0.76, 'SRH 0-0.5:', fontsize=14)
    plt.figtext ( 0.60, 0.76, f'{tot05_SRH:~P}', fontsize=14)
    plt.figtext ( 0.55, 0.74, 'SRH 0-1:', fontsize=14)
    plt.figtext ( 0.59, 0.74, f'{tot1_SRH:~P}', fontsize=14)
    plt.figtext (0.55, 0.72, 'SRH 1-3:', fontsize=14)
    plt.figtext (0.59, 0.72, f'{tot13_SRH:~P}', fontsize=14)
    plt.figtext ( 0.55, 0.70, 'SRH 0-2:', fontsize=14)
    plt.figtext ( 0.59, 0.70, f'{tot2_SRH:~P}', fontsize=14)
    plt.figtext ( 0.55, 0.68, 'SRH 0-3:', fontsize=14)
    plt.figtext ( 0.59, 0.68, f'{tot_SRH:~P}', fontsize=14)
    plt.figtext ( 0.55, 0.38, 'Bunkers Rgt:', fontsize=14)
    plt.figtext ( 0.61, 0.38, f'{bunk_right_spd:~P} at {bunk_right_dir:~P}', fontsize=14)
    plt.figtext ( 0.55, 0.36, 'Bunkers Lft:', fontsize=14)
    plt.figtext ( 0.61, 0.36, f'{bunk_left_spd:~P} at {bunk_left_dir:~P}', fontsize=14)
    plt.figtext ( 0.55, 0.34, '0-6km Mean:', fontsize=14)
    plt.figtext ( 0.61, 0.34, f'{wind_mean[0]:~P}', fontsize=14)
    plt.figtext ( 0.79, 0.15, 'Critical Angle:', fontsize=14)
    plt.figtext ( 0.85, 0.15, f'{ca:~P}', fontsize=14)
    plt.figtext ( 0.55, 0.15, 'SR Mean Wind', fontsize=14)
    plt.figtext ( 0.59, 0.13, 'LM', fontsize=12)
    plt.figtext ( 0.63, 0.13, 'MW', fontsize=12)
    plt.figtext ( 0.67, 0.13, 'RM', fontsize=12)
    plt.figtext ( 0.55, 0.12, '0-1km:', fontsize=12)
    plt.figtext ( 0.58, 0.12, f'{left_mover1[0]:~P}', fontsize=12, color='blue')
    plt.figtext ( 0.62, 0.12, f'{wind_mean1[0]:~P}', fontsize=12, color='grey')
    plt.figtext ( 0.66, 0.12, f'{right_mover1[0]:~P}', fontsize=12, color='red')
    plt.figtext ( 0.55, 0.11, '6-9km:', fontsize=12)
    plt.figtext ( 0.58, 0.11, f'{left_mover69[0]:~P}', fontsize=12, color='blue')
    plt.figtext ( 0.66, 0.11, f'{right_mover69[0]:~P}', fontsize=12, color='red')
    plt.figtext ( 0.62, 0.11, f'{wind_6_9m[0]:~P}', fontsize=12, color='grey')
    plt.figtext ( 0.81, 0.23, '0-0.5km Vec', fontsize=10, color='darkblue')
    plt.figtext ( 0.81, 0.22, '0-3km Vec', fontsize=10, color='magenta')
    plt.figtext ( 0.81, 0.21, '0-6km Vec', fontsize=10, color='gold')
    plt.figtext ( 0.67, 0.78, 'SR Mean Wind', fontsize=14)
    plt.figtext ( 0.67, 0.76, '0-1km:', fontsize=14)
    plt.figtext ( 0.70, 0.76, f'{wind_mean1[0]:~P}', fontsize=14)
    plt.figtext ( 0.67, 0.74, '6-9km:', fontsize=14)
    plt.figtext ( 0.70, 0.74, f'{wind_6_9m[0]:~P}', fontsize=14)
    

    
    plt.title( forecastModel + " " \
              + ob.getLocationName() \
#              + "("+ str(ob.getGeometry()) + ")" \
              + ", " + str(ob.getDataTime()) \
              + ", " + "Hour" + " " + str(int(fcstHour)/60/60)
    )

    # An example of a slanted line at constant T -- in this case the 0 isotherm
    l = skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
    m20 = skew.ax.axvline(-20, color='darkblue', linestyle='--', linewidth=2)

    # Draw hodograph
    ax_hod = fig.add_subplot(gs[0, 1])
    h = Hodograph(ax_hod, component_range=spd.max()/units.knots)
    h.add_grid(increment=20, alpha=0.2)
    h.plot_colormapped(u, v, spd)
    origin = np.array([[0, 0, 0],[0, 0, 0]])
    #plt.quiver(*origin, wind_mean, color='grey', scale=21)
    #plt.quiver(*origin,  bunk_left_dir, bunk_left_spd, color='grey', scale=21)
    #plt.quiver(*origin, u_storm, v_storm, color='grey', scale=21)
    
    # Draw vector hodograph
    ax_hod = inset_axes(ax_hod, '25%', '25%', loc=4)
    h = Hodograph(ax_hod, component_range=spd.max()/units.knots)
    h.add_grid(increment=20, alpha=0.4)
    h.plot(u, v, color='grey', alpha=0.0)
    #h.plot(right_mover, 'ko', color='red')
    origin = np.array([[0, 0, 0],[0, 0, 0]]) # origin point
    plt.quiver(*origin, u_shear005, v_shear005, color='darkblue', scale=21)
    plt.quiver(*origin, u_shear03, v_shear03, color='magenta', scale=21)
    plt.quiver(*origin, u_shear06, v_shear06, color='gold', scale=21)
    
     ######## Save the plot
   
    plt.savefig((output_dir+'/Soundings/GFS/CHA/'+str(int(fcstHour)/60/60)+'.png'),bbox_inches='tight',pad_inches=0.1)
    fcst_hr = str(0)
    plt.clf()