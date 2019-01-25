#!/usr/bin/python

#raw_analysis.py performs initial analysis of raw
#data collected in flight behavior rig during sun experiments

###USAGE 
#python raw_analysis.py raw_date_range.txt

#it saves a .pck file for each experiment with 
##summary data for that experiment.

##A text file (e.g. 'raw_date_range.txt' specifies the dates, or range of dates,
#over which to analyze data.
#A single date specified as YYYY.MM.DD
#range specified as YYYY.MM.DD-YYYY.MM.DD


#######################################
import numpy as np
import os
import pylab
import pdb
import pickle
from py_utilities import circ
import time
import datetime
import pytz
from py_utilities import fp_library as fpl
import matplotlib.pyplot as pyplot
#import scipy.stats.binned_statistic as binst
import scipy.stats.morestats as st
import scipy.signal as signal
from py_utilities import tw_plot_library3 as plt
from py_utilities import fly_plot_basics as fpb
from py_utilities import tw_filehandling as fh
from py_utilities import tw_calc_library as calc
#import plot_sum16 as plot_sum16
import matplotlib.gridspec as gridspec
#import von_mises_fit as vmfit
pylab.ion()
fname_list=list()
gain_list=list()
time_list=list()
combgaindata={}
#gainStepsPerRadSec=23375
plotbnds=[]
combgaindata[0]=[]
ANAL_ALL_FLAG=0
SAMPLERATE=40.0
PLOT_SCATTER_POL=1
POL_SENSOR_MINIMUM=75
PAIR_CRITERION_MINUTES=22
FIT_FEEDBACK_FLAG=0
SET_OFFSET_SHIFT_FLAG=0


NUMBER_OF_MOTOR_ADJUSTMENTS=2
CALC_VEC_STRENGTH_BY_CONVOLUTION=1
REALIGN_EACH_HALF=False
OVERWRITE_OFFSET_FLAG=True
OVERWRITE_OFFSET_VL=5.35-np.pi

ARCHIVE_BOUNDARY=pylab.date2num(datetime.datetime(2014,1,1))

input_path='/users/tim/Dropbox/flystuff/uo_data/'
CHECK_CALIBRATION_FLAG=0
CALC_CALIBRATION_FLAG=0
CALC_SPECTRUM=0
#to do, apply gaussian filter, reduce step_len to 0.15 s #experiment with different filter len
#is there a good way of calculating a weighted mean vector strength?
FILTER_LEN_IN_SEC=30.0
STEP_LEN_IN_SEC=0.1
VEC_STRENGTH_THRESHOLD=0.8


ZERO_MOTOR_POSITION=1790.
TRUNC_DATA_FLAG=0
TRUNC_DATA_TIME_IN_SEC=240
MAX_SUN_POSITION=4060.
###############################################
class Sun_Anal():
    #calls prep_data_fxn for each file
    def __init__ (self,input_args):
        #self.anal_slip_flag=0
        #self.RESET_LENGTH_FLAG=0

        self.day_array,self.exp_type,self.anal_type=fh.read_days(input_args[1])

        self.PLOT_FLAG=1

        if len(input_args)<3:
            self.fly_limiter=[]
        else:
            if input_args[2]=='no_plot':
                self.PLOT_FLAG=0
                self.fly_limiter=[]

            else:
                self.fly_limiter=input_args[2]



       


    def plot_traj_by_vec_strength(self,ax,mn_drxn,vec_strength,**kwargs):

        ax.plot(mn_drxn,vec_strength)
        ax.get_yaxis().set_ticks([])
        ax.title.set_visible(False)
        ax.get_xaxis().set_ticklabels([])
        ax.spines['polar'].set_color('none')

        cr_mn=mn_drxn
        cr_ln=vec_strength




        binvls=[72,20]
        rangevls=[[0,2*np.pi],[0,1]]



        if cr_mn.size:

            heatmap, xedges, yedges = np.histogram2d(cr_mn,cr_ln, bins=binvls,range=rangevls)


            xbins=np.linspace(0,2*np.pi,binvls[0])
            ybins=np.linspace(0,1,binvls[1])
            theta,r=pylab.meshgrid(ybins,xbins)

            heatmap_norm=self.norm_heat_map(heatmap)


            heat_map_data=self.save_heat_map_data(heatmap_norm,xedges,yedges,r,theta)

            #plt.polar_heat_map(ax_traj_im,heat_map_data,aligned=0,plot_power_scale=5)
                    #self.datvls['sumstats']['exp']['heatmap_vls']={}




    def run_analysis(self):

        for crday in np.array(self.day_array):
            pylab.close('all')
            self.crday=pylab.num2date(crday)
            self.process_data_by_day()

    def process_data_by_day(self):

        self.datapath=fh.make_data_path(input_path,self.crday)


        dirlist=fh.get_directories(self.datapath)

        if self.exp_type =='open':
            path_end='/oloop/'
        else:
            path_end='/cloop/'

        if dirlist:
            if not self.fly_limiter:
                for crdir in dirlist:
                    self.in_data_path=self.datapath+crdir+path_end
                    self.process_data_by_fly()
            else:
                crdir=self.fly_limiter
                self.in_data_path=self.datapath+crdir+path_end
                self.process_data_by_fly()

    def process_data_by_fly(self):
        self.datvls={}
        files = sorted(os.listdir(self.in_data_path))
        self.exec_file=os.path.abspath(__file__)
        self.txtfiles = sorted([f for f in files if f[-3:]=='txt'])
        self.crdt={}

        for self.crfile in self.txtfiles:
            self.get_data_by_file()
            
            if self.params:
                try:
                    self.datpath=self.raw_save_file_name





                    self.crdt=fh.open_pickle(self.datpath)


                    self.pckname=self.in_data_path+str(self.fname)+'_rawdata.pck'
                except:

                    self.pckname=[]

                    print 'no pickle'
            else:
                #pdb.set_trace()
                self.crdt={}
                self.pckname=[]
                print 'no pickle'


            self.make_example_figure()
            #except:
            #    print 'no file'
                #ps.make_plots()

                #self.make_plots_by_file()

    def get_data_by_file(self):
        mod_file_name=os.path.join(self.in_data_path,'modified_raw_data.pck')
        self.fname = os.path.join(self.in_data_path,self.crfile)
        pck_fname=self.fname.strip('.txt')+'.pck'
        try:
            self.params=fh.open_pickle(pck_fname)
        except:
            print('open pickle error')
        #self.test_if_variable_gain()

        #self.add_params()

        try:
            if not os.path.isfile(mod_file_name):
                self.tmpdat=(np.array(np.genfromtxt(self.fname)))
            else:
                self.tmpdat=fh.open_pickle(mod_file_name)
        except:
            print('file name error')
            self.tmpdat=[]


        self.prep_data_fxn()

        if self.crdata:
            #if self.exp_type=='closed':

            self.anal_closed_stroke_data()
            if self.anal_type=='paired':
                self.datvls[self.pairind]={}

                self.datvls[self.pairind]['sumstats']=self.sumstats
                self.datvls[self.pairind]['splitstats']=self.splitstats
            else:
                self.datvls['PAIR_FLAG']=0
                self.datvls['sumstats']=self.sumstats
                self.datvls['splitstats']=self.splitstats

                self.datvls['params']=self.params

            # elif self.exp_type=='open':
            #     self.datvls['dat']=self.indt
            #     self.datvls['params']=self.params
            #     self.datvls['bndtimes']=self.bndtimes
            #     self.raw_save_file_name=self.fname.strip('.txt')+'_' + 'rawdata.pck'

            #     fh.save_to_file(self.exec_file,self.raw_save_file_name,self.datvls)
            #    self.datvls=self.open_sumdt

    def test_if_variable_gain(self):
        #notion is that variable_gain_experiment should have

        if 'variable_gain_duration' in self.params.keys():
            if 'variable_closed_loop_duration_sec' in self.params.keys():
                if self.params['variable_closed_loop_duration_sec']:
                    if 'pre_params' in self.params.keys():
                        self.variable_gain_flag=0
                    else:
                        self.variable_gain_flag=1
                else:
                    self.variable_gain_flag=0
            else:
                self.variable_gain_flag=0
        else:
            self.variable_gain_flag=0

    def test_if_ersatz_sun(self):
        if 'ersatz_sun' in self.params.keys():
            if self.params['ersatz_sun']:
                self.ersatz_sun_flag=1
            else:
                self.ersatz_sun_flag=0

    def add_params(self):

        self.params['filter_len_in_sec']=FILTER_LEN_IN_SEC
        self.params['step_len_in_sec']=STEP_LEN_IN_SEC
        self.params['vec_strength_threshold']=VEC_STRENGTH_THRESHOLD
    #this does some basic analysis of data
    #
    def anal_closed_stroke_data(self):

        self.outdict={}
        self.outdict['exp']={}
        self.outdict['ctrl']={}

        self.outdict['mot_rad']=self.crdata['adjusted_raw_motor_position_in_rad']

        self.outdict['mot_deg']=calc.rad_to_deg(self.crdata['adjusted_raw_motor_position_in_rad'])
        self.outdict['exp']['sumdeg']=self.outdict['mot_deg'].tolist()
        self.outdict['float_time']=self.crdata['time']
        self.outdict['raw_time']=self.crdata['net_time']
        self.outdict['period']=self.crdata['period']
        #self.outdict['corr_dict']=self.crdata['corr_dict']

        self.analyze_degree_data()


    def analyze_degree_data(self):
        self.sumstats={}
        self.splitstats={}
        deg_per_bin=10
        self.degbins = np.arange(0,370,deg_per_bin)

        #bnd_times=self.crdata['boundary_times']
        timevls=np.array(self.crdata['net_time'])
        #first find which boundary refers to 'positive' run
        #for i,crexp in enumerate( self.params['closed_loop_experiment_list']):
         #   if crexp == 'positive':
          #      self.closed_exp_index=i

        if self.bndtimes:
            if self.crday>datetime.datetime(2017,4,26,tzinfo=pytz.utc):
                self.assign_sumstats_modified(deg_per_bin,timevls)
            else:
                self.assign_sumstats_original(deg_per_bin,timevls)
        else:
            tst=1

    def assign_sumstats_crux(self,deg_per_bin,timevls):



        #nonaninds=~np.isnan(timevls)
        self.analinds=np.intersect1d(np.where(timevls>self.bndtimes[self.crkey][0])[0],np.where(timevls< self.bndtimes[self.crkey][1])[0])
        #self.analinds=nonaninds[tmpanalinds]

        crdeg=np.array(self.outdict['exp']['sumdeg'])[self.analinds]


        hstout_dt=calc.make_hist_calculations(crdeg,self.degbins)
        hstout_dt['start_rad']=calc.rad_to_deg(circ.circmean(calc.deg_to_rad(crdeg[0:2])))

        self.make_sumstats_values(hstout_dt)

        self.sumstats[self.crkey]['steps_per_rotation']=self.outdict['period']
        #TKTK save something here indicating index where experiment switched
        self.sumstats[self.crkey]['mean_deg']=calc.rad_to_deg(self.sumstats[self.crkey]['mnrad'])
        self.sumstats[self.crkey]['fname']=self.fname
        self.rad_per_bin=deg_per_bin*np.pi/180
        self.sumstats[self.crkey]['rad_per_bin']=deg_per_bin*np.pi/180
        self.xrad=calc.deg_to_rad(self.degbins)
        self.sumstats[self.crkey]['xrad']=calc.deg_to_rad(self.degbins)
        self.sumstats[self.crkey]['raw_time']=np.array(self.crdata['time'])[self.analinds]
        try:
            self.sumstats[self.crkey]['time_in_min']=np.array(self.crdata['net_time'])[self.analinds]/60.0
        except:
            print('time error')
            #pdb.set_trace()
        #self.sumstats[self.skey]['mot_vel']=np.array(self.crdata['mot_vel'])[self.analinds]
        self.sumstats[self.crkey]['lftwng']=np.array(self.crdata['lftwng'])[self.analinds]
        self.sumstats[self.crkey]['rtwng']=np.array(self.crdata['rtwng'])[self.analinds]

        #self.sumstats[self.skey]['pol_sensor']=np.array(self.crdata['pol_sensor'])[self.analinds]
        self.sumstats[self.crkey]['mot_position_in_rad']=np.array(self.crdata['adjusted_raw_motor_position_in_rad'])[self.analinds]

        if 'gain_coefficient' in self.crdata.keys():

            self.sumstats[self.crkey]['gain_coefficient']=np.array(self.crdata['gain_coefficient'])[self.analinds]

        self.sumstats[self.crkey]['mot_deg']=self.outdict['mot_deg'][self.analinds]

        self.sumstats[self.crkey]['mot_rad']=np.array(self.outdict['mot_rad'])[self.analinds]
        try:

            self.sumstats[self.crkey]['cloop_duration_in_min']=np.ceil(self.sumstats[self.skey]['time_in_min'][-1]-self.sumstats[self.skey]['time_in_min'][0])


                #pdb.set_trace()
            #this is hack specific for variable_gain_exps with 10/15
            #training_time=10
            #try:
             #   self.sumstats['exp']['switch_time']=self.bndtimes[0][0]/60.+10
            #except:
             #   tst=1
            #self.sumstats['exp']['switch_time']=self.sumstats['exp']['time_in_min'][analbnds[1][0]]
            self.calculate_continuous_heading()
            self.sumstats[self.crkey]['displacement_traj']=self.calc_displacement_trajectory()

            if 'len_vector_lst' in self.sumstats[self.crkey].keys():
                cr_ln=np.array(self.sumstats[self.crkey]['len_vector_lst'])

                cr_mn=np.array(self.sumstats[self.crkey]['mn_vector_lst'])

                threshinds=np.where(np.array(cr_ln)>VEC_STRENGTH_THRESHOLD)[0]

                out_dt=calc.make_hist_calculations(calc.rad_to_deg(cr_mn[threshinds]),self.degbins)

                self.sumstats[self.crkey]['vecmn']=out_dt['mnrad']
                self.sumstats[self.crkey]['vechst']=out_dt['normhst']
                self.sumstats[sef.crkey]['vecvar']=out_dt['circvar']
        except:
            tst=1



    def assign_sumstats_original(self,deg_per_bin,timevls):
        for expkey in ['stripe','sun']:
            self.sumstats[expkey]=[]


            self.sumstats[expkey].append({})


            stripe_inds=[0]
            sun_inds=[1]


        if len(self.bndtimes)>1:

            for ind,crbin in enumerate(self.bndtimes[0:2]):

                self.assign_sumstats_crux(deg_per_bin,timevls)
            

            if TRUNC_DATA_FLAG==0:
                self.raw_save_file_name=self.fname.strip('.txt')+'_' + 'rawdata.pck'
            else:
                self.raw_save_file_name=self.fname.strip('.txt')+'_' + 'TRUNC'+str(TRUNC_DATA_TIME_IN_SEC)+'rawdata.pck'


            fh.save_to_file(self.raw_save_file_name,self.sumstats,save_executable=self.exec_file)


        else:
            
            self.raw_save_file_name=[]
            self.PLOTFLAG=False
            return











    def assign_sumstats_modified(self,deg_per_bin,timevls):
        self.sumstats={}

        try:
            exp_type=self.params['exp_type']
        except:
            exp_type='random'




        for self.crkey,crbin in enumerate(self.bndtimes):

            self.assign_sumstats_crux(deg_per_bin,timevls)
            if TRUNC_DATA_FLAG==0:
                self.raw_save_file_name=self.fname.strip('.txt')+'_' + 'rawdata.pck'
            else:
                self.raw_save_file_name=self.fname.strip('.txt')+'_' + 'TRUNC'+str(TRUNC_DATA_TIME_IN_SEC)+'rawdata.pck'

           
            fh.save_to_file(self.raw_save_file_name,self.sumstats,save_executable=self.exec_file)


       
    def make_sumstats_values(self,indt):
        self.sumstats[self.crkey]={}
        for crkey in ['normhst','mnrad','circvar','circvar_mod360','mnrad_360','mnrad_max180','start_rad']:
            
            self.sumstats[self.crkey][crkey]=indt[crkey]
        
        self.sumstats[self.crkey]['entropy']=calc.entropy(self.sumstats[self.crkey]['normhst'])

    def calculate_continuous_heading(self):
        mn_lst=[]
        lenlst=[]
        mn_dot_product=[]
        vec_time_list=[]

        crdt=self.sumstats[self.skey]
        #nalinds=np.where(crdt['time_in_min']>=rdt['switch_time'])[0]

        analinds=np.where(crdt['time_in_min']>=crdt['time_in_min'][0])[0]
        mot_dt=crdt['mot_rad'][analinds]
        time_dt=crdt['time_in_min'][analinds]

        #FILTER_LEN_IN_PTS=FILTER_LEN_IN_SEC*SAMPLERATE

        STEP_LEN_IN_PTS=STEP_LEN_IN_SEC*SAMPLERATE

        if analinds.size:
            lenlst,mn_lst,vec_time_list=calc.get_dynamic_vec_strength(mot_dt,time_dt,'convolution','360',FILTER_LEN_IN_SEC)
        else:
            return


        try:
            mn_array=np.array(mn_lst)
        #calculate the overall vector strength
        except:
            print('mean calculation error')

        mn_array[mn_array<0]=mn_array[mn_array<0]+2*np.pi
        self.sumstats[self.skey]['mn_vector_lst']=mn_array
        #doubled_mn_array=np.mod(2*heading_data,2*np.pi)
        #self.sumstats[self.skey]['overall_vec_strength']=1-circ.circvar(doubled_mn_array)

        self.sumstats[self.skey]['len_vector_lst']=lenlst

        self.sumstats[self.skey]['mn_dot_product_lst']=mn_dot_product

        self.sumstats[self.skey]['vec_time_lst']=vec_time_list
        #self.sumstats[self.skey]['num_vec_overlap']=FILTER_LEN_IN_PTS/STEP_LEN_IN_PTS


    def make_plots_by_file(self):

        if self.PLOTFLAG:
            self.plot_combdata()
            self.figname=self.pckname.split('/')[-1].strip('_sumdata.pck')
            self.save_file_name=self.fname.strip('.txt')+'_' + self.anal_type + '_sumdata.pck'

            if os.path.isfile(self.save_file_name):

                os.remove(self.save_file_name)
            if os.path.isfile(self.fname.strip('.txt')+'_' + 'UNPAIRED_sumdata.pck'):
                os.remove(self.fname.strip('.txt')+'_' + 'UNPAIRED_sumdata.pck')
            
            fh.save_to_file(self.exec_file,self.save_file_name,self.datvls,fig=self.fig,figname=self.figname)

    #creates self.allpairs, which is a list of lists with each pair of data as a f
    def get_pairs(self):
        times=[]
        self.allpairs=[]
        for crfile in self.txtfiles:
            crtimestr=crfile.split('_')[2].strip('.txt')[0:4]

            times.append(calc.get_time(crtimestr))
        #clumsy way of doing this for each time in list that is not the last time #find all times that are greater than current time and within boundary.
        self.sorted_times=np.array(times)
        for i,crtime in enumerate(self.sorted_times):
            diffvls=self.sorted_times-crtime
            ind1=np.where(diffvls>0)[0]
            ind2=np.where(diffvls<PAIR_CRITERION_MINUTES)[0]
            finalind=np.intersect1d(ind1,ind2)

            if len(finalind):

                self.allpairs.append([i,finalind[0]])
#####
#####
#combines data into a global dict crdata
    def prep_data_fxn(self):

        alldata={}
        fname_list.append(self.fname)

        if self.params:
            self.PLOTFLAG=1
            self.parse_indata()

            self.get_exp_boundary_times()

            try:
                self.bndtimes
                BND_SET=1
            except:
                print('bound times error')
                BND_SET=0


            self.adjust_motor_data_to_correct_orientation()

            #self.indt['adjusted_raw_motor_position']=self.indt['raw_motor_position']
                    #self.indt['motor_position_radians']=2*np.pi*((np.mod(self.indt['adjusted_raw_motor_position'],self.calculated_motor_period))/self.calculated_motor_period)
            #time_list.append(int(self.fname[-10:-6]))
            if np.array(np.nonzero(np.isnan(self.indt['time']))).size:
                strt_vl=0
                end_vl=np.max(self.indt['time'])
            else:
                alldata['total_closed_loop_time']=0
            #returns crdata
            self.interpolate_times()

            self.crdata['period']=MAX_SUN_POSITION

        else:
            
            self.PLOTFLAG=False

#rewriting this so, whether closed loop or open loop,
    #self.bndtimes returns nx2 arraay wtih start and stop
    #called by prep_data_fxn
    def get_exp_boundary_times(self):
        #dealing with thorny issue of whether STOP was added to end of open loop data
        self.bndtimes=[]
        # NAN_AT_STOP=True
        nanvls=np.nonzero(np.isnan(self.indt['net_time']))[0]
        
        if len(nanvls)>0:

            self.bndtimes.append([0,self.indt['net_time'][nanvls[0]-1]])
            for crind in np.arange(len(nanvls)-1):

                self.bndtimes.append([self.indt['net_time'][nanvls[crind]+1],self.indt['net_time'][nanvls[crind+1]-1]])


        #self.bndtimes=[self.indt['net_time'][0],self.indt['net_time'][-1]]
    def parse_indata(self):
        #pdb.set_trace()

        self.indt={}
        try:
            self.indt['time']=self.tmpdat[:,0]
        except:
           print ('parse error')
        #self.indt['mot_vel']=self.tmpdat[:,1]
        self.indt['net_time']=(self.indt['time']-self.indt['time'][0])

        self.indt['raw_motor_position']=self.tmpdat[:,1]


        self.indt['lftwng']=self.tmpdat[:,2]*(180/np.pi)
        self.indt['rtwng']=self.tmpdat[:,3]*(180/np.pi)

        #if self.crday>datetime.datetime(2014,2,1,tzinfo=pytz.utc):
         #   self.indt['lftwng']=-self.tmpdat[:,3]*(180/np.pi)
          #  self.indt['rtwng']=-self.tmpdat[:,4]*(180/np.pi)
        #self.indt['pol_sensor']=self.tmpdat[:,5]
        #if np.shape(self.tmpdat)[1]>6:
         #   self.indt['gain_coefficient']=self.tmpdat[:,6]
        #else:
         #   self.indt['gain_coefficient']=[]

    #for stripe, 1790 is 0, 3844 is 180.
    #easiest conversion is to convert_to_radians immediately and then plot from there

    #for values greater or equal to ZERO_MOTOR_POSITON, map those to
    #mod ZERO_MOTOR_POSITION

    def adjust_motor_data_to_correct_orientation(self):
        invls=np.copy(self.indt['raw_motor_position'])
        invls_rad=(self.indt['raw_motor_position']/MAX_SUN_POSITION)*2*np.pi
        offset_in_rad=(ZERO_MOTOR_POSITION/MAX_SUN_POSITION)*2*np.pi
        corrected_rad=invls_rad-offset_in_rad
        try:
            neg_inds=np.where(corrected_rad<0)
            corrected_rad[neg_inds]=corrected_rad[neg_inds]+2*np.pi
        except:
            tst=1

        self.indt['adjusted_raw_motor_position_in_rad']=corrected_rad

        #self.indt['adjusted_raw_motor_position'][highinds]=np.mod(invls[highinds],zero_pos)
        #self.indt['adjusted_raw_motor_position'][lowinds]=invls[lowinds]+zero_pos



        #self.params['adjusted_start_angle']=self.params['start_pos_deg']+calc.rad_to_deg(offset_in_radians)
        #self.params['new_horiz_offset']=new_horiz_offset

        #self.params['fit_covmat']=covmat

    #for linear polarizer, this now returns offset of minimum from zero
    def make_fit_calculation(self,xvls,yvls,**kwargs):
        if 'filter_type' not in self.params:
            self.params['filter_type']='linear'
        if 'filter_type' in self.params:

            if 'ersatz_sun' in self.params:
                if 'ersatz_sun' in self.params:

                    sun_flag=1

                try:
                        #first take a histogram with 36 mean values, then calculate the peak using weighted mean
                    raw_yvls=kwargs['raw_yvls']
                    bins=np.linspace(0,2*np.pi,48)

                    if sun_flag:
                        goodinds= np.intersect1d(np.where(raw_yvls>0)[0],np.where(raw_yvls<12)[0])

                    fitvls=vmfit.make_fit(xvls,raw_yvls,plot_fit=True)
                    #fitvls=calc.vonmises_fit_new(xvls[goodinds],yvls[goodinds]-np.min(yvls[goodinds]))
                    fitout=fitvls[0][1]

                    if OVERWRITE_OFFSET_FLAG:
                        fitout=OVERWRITE_OFFSET_VL
                    if fitout>np.pi:
                        fit1=fitout-np.pi
                    else:
                        fit1=fitout+np.pi
                    fit0=1


                    polarizer_correction_value=0
                    
                except:
                    print ('fit error')
                    tst=1
            else:
                if self.params['filter_type']=='linear':
                    try:
                        fit0,fit1=calc.cos_ls_fit(xvls,yvls,num_cycles=2,plot_flag=1)
                    except:
                        fit0=0
                        fit1=0
                    offset_value_if_negative=np.pi/2
                    polarizer_correction_value=np.pi/2
                elif self.params['filter_type']=='circular':
                    fit0,fit1=calc.cos_ls_fit(xvls,yvls,num_cycles=1)
                    offset_value_if_negative=np.pi
                    polarizer_correction_value=0



                else:
                    fit0,fit1=calc.cos_ls_fit(xvls,yvls,num_cycles=2,plot_flag=1)
                    offset_value_if_negative=np.pi/2
                    polarizer_correction_value=np.pi/2


            if fit0<0:
                horiz_offset=fit1+offset_value_if_negative
            else:
                horiz_offset=fit1

            if 'intensity_cue' in self.params.keys():
                if self.params['intensity_cue']:
                    intensity_anal=1
                else:
                    intensity_anal=0
            else:
                intensity_anal=0

            if 'experiment_type' in self.params.keys():
                if self.params['experiment_type']=='wedge_45':
                    wedge_anal=1
                else:
                    wedge_anal=0
            else:
                wedge_anal=0



            if intensity_anal or wedge_anal:

                new_horiz_offset=polarizer_correction_value+self.adjust_motor_for_intensity_cues(xvls,yvls,horiz_offset)

            else:

                new_horiz_offset=polarizer_correction_value+horiz_offset

        return new_horiz_offset

    def adjust_motor_for_intensity_cues(self,in_xvls,in_yvls,horiz_offset):
        radian_bnds=[]
        crmed=[]
        adj_vls=(in_xvls-horiz_offset)+2*np.pi

        new_adj_vls=np.mod(adj_vls,2*np.pi)

        if self.params['filter_type']=='linear':
            #these are the bounds over which to search

            radian_bnds.append([np.pi/2-.2, np.pi/2+.2])
            radian_bnds.append([3*np.pi/2-.2 ,3*np.pi/2+.2])
            for crbnds in radian_bnds:
                crinds=np.intersect1d(np.where(new_adj_vls>crbnds[0])[0],np.where(new_adj_vls<crbnds[1])[0])
                crmed.append(np.median(in_yvls[crinds]))
            #make the 270 degree value the lower value
            if crmed[1]>crmed[0]:
                new_horiz_offset=horiz_offset+np.pi
            else:
                new_horiz_offset=horiz_offset
        elif self.params['filter_type']=='circular':
           new_horiz_offset=horiz_offset
        return new_horiz_offset
        return new_minimum_value





    def calculate_gain_minimum(self):
        #find inds where self.indt['gain_oeffienct is minimum
        #what is the min and max, of motor_position_radians(mod180) of those values
        #what is midpoint of min and max- that is the minimum

        minvl=self.params['minimum_gain']
        mininds=np.where(self.indt['net_time']<self.params['variable_gain_duration']) and np.where(self.indt['gain_coefficient']==minvl)
        mod_motor_position=np.mod(self.indt['motor_position_radians'],np.pi)
        motpos=np.unique(mod_motor_position[mininds[0]])
        return 0.5*circ.circmean(2*motpos)


    #called by prep_data_fxn
    #rewritten so that times between nan are not interpolated
    def interpolate_times(self):
        self.crdata={}
        keys=['net_time','adjusted_raw_motor_position_in_rad','lftwng','rtwng','time']
        try:
            if self.indt['gain_coefficient'].size:
                keys.append('gain_coefficient')
        except:
            tst=1
        for crkey in keys:
            self.crdata[crkey]=[]

        corr_time=self.indt['net_time']
        interp_bnds=np.nonzero(np.isnan(self.indt['net_time']))[0]

        for crkey in keys:
            #for i,crbnd in enumerate(interp_bnds):
            init_time=np.min(corr_time[~np.isnan(corr_time)])
            end_time=np.max(corr_time[~np.isnan(corr_time)])
            out_time=np.arange(init_time,end_time,1/SAMPLERATE)
            try:

                self.crdata[crkey]=self.crdata[crkey]+np.interp(out_time,corr_time,self.indt[crkey]).tolist()
            except:
                
                print ('no data')

    #takes crdata and transforms it into open_dt
    #open_sumdata[0-n] represents data from different conditions with positive and negative
    #open_sumdata[0]['pos'],open_sumdt[0]['neg']....
    def arrange_open_loop_data(self):
        self.open_sumdt=[]
        exparray=np.array(self.params['open_loop_exp_list'])
        for crlistind,cr_pos_gain_value in enumerate(np.unique(abs(exparray))):
            self.POSCT=0
            self.NEGCT=0
            self.open_sumdt.append({'pos':{},'neg':{}})
            pos_expinds=np.where(exparray==cr_pos_gain_value)[0]
            neg_expinds=np.where(exparray==-cr_pos_gain_value)[0]
            for self.crexpind in pos_expinds:
                #crbnds=get_current_bnds(self.bnds, crind,self.params)
                self.cut_open_loop_data(crlistind,'pos')
                self.POSCT=self.POSCT+1
            for self.crexpind in neg_expinds:
                #crbnds=get_current_bnds(self.bnds, crind,self.params)
                self.cut_open_loop_data(crlistind,'neg')
                self.NEGCT=self.NEGCT+1

    #called by arrange_open_loop_data
    def cut_open_loop_data(self,crlistind,open_drxn):
        time_bnds=self.get_index_bnds(self.crdata['net_time'],self.crexpind)

        sensor_bnds=np.where(np.array(self.crdata['pol_sensor'])[time_bnds]>POL_SENSOR_MINIMUM)
        try:
            crbnds=time_bnds[sensor_bnds]
        except:
            pdb.set_trace()
        tmpdt={}

        tmpdt['open_loop_gain']=np.array(self.params['open_loop_exp_list'])[self.crexpind]
        tmpdt['lftwng']=np.array(self.crdata['lftwng'])[crbnds]
        tmpdt['rtwng']=np.array(self.crdata['rtwng'])[crbnds]
        tmpdt['pol_sensor']=np.array(self.crdata['pol_sensor'])[crbnds]
        tmpdt['raw_motor_position']=np.array(self.crdata['adjusted_raw_motor_position'])[crbnds]
        tmpdt['motor_position_radians']=np.array(self.crdata['motor_position_radians'])[crbnds]

        tmpdt['net_time']=np.array(self.crdata['net_time'])[crbnds]-self.crdata['net_time'][crbnds[0]]
        tmpdt['wingdiff']=tmpdt['lftwng']-tmpdt['rtwng']
        if open_drxn=='pos':

            self.open_sumdt[crlistind]['pos'][self.POSCT]=tmpdt
        elif open_drxn=='neg':
            self.open_sumdt[crlistind]['neg'][self.NEGCT]=tmpdt

    #called by cut_open_loop_data
    def get_index_bnds(self,timevls,**kwargs):
        timeinds=[]
        plot_multiple=0
        if 'get_multiple_bnds' in kwargs:
            if kwargs['get_multiple_bnds']:
                plot_multiple=1
        if plot_multiple:

           for crbndind in np.arange(len(self.bndtimes)):
                crbnds=self.bndtimes[crbndind]
                bnds1=np.where(np.array(timevls)>crbnds[0])[0]
                bnds2=np.where(np.array(timevls)<crbnds[1])[0]
                crtminds=np.intersect1d(bnds1,bnds2)
                timeinds.append(crtminds)

                self.analinds=timeinds
                return(self.analinds)
        else:
            try:
                bnds1=np.where(np.array(timevls)>self.bndtimes[0][0])[0]
                bnds2=np.where(np.array(timevls)<self.bndtimes[0][1])[0]
                self.analinds=np.intersect1d(bnds1,bnds2)
                return self.analinds
            except:
                print 'bounds problem'

###
###
#functions for deatling with open loop data
###
###
    #this function will take every open loop experiment and ma


#called by process_data_by_fly


    def make_offset_corrections(self,gainind,crkey,crexpnum):
        KEYLIST=['pol_sensor','wingdiff']
        for current_data in KEYLIST:
            if self.keyind==0:
                crvls=self.crdict[current_data]
            else:
                crvls=self.crdict[current_data][::-1]
            if(self.offset<0):
                    corr_vls=crvls[-self.offset:]
                    #data needs to be moved to right
            elif(self.offset>0):
                    nanarray=np.empty(self.offset)
                    nanarray[:]=np.nan
                    corr_vls=np.concatenate((nanarray,crvls))
            else:
                corr_vls=crvls
            self.open_sumdt[gainind][crkey][crexpnum][current_data+'_offcorr']=corr_vls
            self.timediff=self.open_sumdt[0]['pos'][0]['net_time'][1]-self.open_sumdt[0]['pos'][0]['net_time'][0]
            corr_time=np.linspace(0,len(corr_vls)*self.timediff,num=len(corr_vls))
            self.open_sumdt[gainind][crkey][crexpnum]['net_time_offcorr']=corr_time

            if(current_data==KEYLIST[0]):
                self.lenlist.append(len(corr_vls))

    def reduce_to_single_cycle(self):
        NUM_FULL_CYCLES=1
        for gainind,crdt in enumerate(self.open_sumdt):
            cr_drxn_key='pos'

            mnpol=crdt[cr_drxn_key]['mn_pol_sensor']
            #L8=int(len(mnpol)/16)
            #calc_pol=mnpol[L8:-L8]-np.mean(mnpol[L8:-L8])
            calc_pol=mnpol-np.mean(mnpol)
            output_vls=np.correlate(calc_pol,calc_pol,'same')
            #need to determine period of mean value and then relate this to ave_time which gives time difference
            self.make_cut_inds(output_vls)
            self.rearrange_data_by_cut_points()

    def make_cut_inds(self,output_vls):
        import peakdetect as peakdetect
        self.cut_inds=[]

        peak_data=peakdetect.peakdetect(output_vls)
        if (len(peak_data[0])==3):
            full_cycle_length=peak_data[0][2][0]-peak_data[0][0][0]
        #use minima
        elif len(peak_data[1])==2:

            full_cycle_length=2*peak_data[1][1][0]-peak_data[1][0][0]
        elif len(peak_data[0])==2:
            full_cycle_length=2*peak_data[0][1][0]-peak_data[0][0][0]
        number_of_full_cycles=np.floor(len(self.open_sumdt[0]['ave_time'])/full_cycle_length)
        start_cut_ind=0
        for cycle_number in np.arange(0,number_of_full_cycles):
            self.cut_inds.append([start_cut_ind, start_cut_ind+full_cycle_length])
            start_cut_ind=start_cut_ind+full_cycle_length


    #start by just cutting the means
    def rearrange_data_by_cut_points(self):
        self.CUT_INITIAL_TIME=2
        for gainind,crdt in enumerate(self.open_sumdt):
            for cr_drxn_key in ['pos','neg']:

                for crdatakey in ['mn_wingdiff','mn_pol_sensor']:
                    tmplist=[]
                    cutkey=crdatakey+'_cut'

                    self.open_sumdt[gainind][cr_drxn_key][cutkey]={}
                    for cut_ind_num in np.arange(0,len(self.cut_inds)):
                        crinds=self.cut_inds[cut_ind_num]
                        tmpvl=crdt[cr_drxn_key][crdatakey][crinds[0]:crinds[1]]

                        self.open_sumdt[gainind][cr_drxn_key][crdatakey+'_cut'][cut_ind_num]=tmpvl
                        tmplist.append(tmpvl)
                    crmn,stdvl=calc.mnstd(tmplist)

                    self.open_sumdt[gainind][cr_drxn_key][crdatakey+'_cut']['mean']=crmn

    def calc_wingdiff_mn_std(self):
        MNCT=0
        mnvl=[]
        DATLEN=np.floor(np.median(self.lenlist))
        self.ave_time=np.linspace(0,DATLEN*self.timediff,num=DATLEN)

        for gainind,crgaindata in enumerate(self.open_sumdt):
            for cr_drxn_key in ['pos','neg']:
                for cr_datatype_key in ['pol_sensor_offcorr','wingdiff_offcorr']:
                    tmplist=[]
                    crdt=crgaindata[cr_drxn_key]
                    indlist=[x for x in crdt.keys() if type(x)  is int]
                    for crexpnum in indlist:
                        tmpvl=crdt[crexpnum][cr_datatype_key]
                        if len(tmpvl)>=DATLEN:
                            tmplist.append(tmpvl[0:DATLEN])
                        else:
                            difflen=DATLEN-len(tmpvl)
                            tmparray=np.empty(difflen)
                            tmparray[:]=np.nan
                            tmplist.append(np.concatenate((tmpvl,tmparray)))
                    crmn,stdvl=calc.mnstd(tmplist)
                    if cr_datatype_key=='wingdiff_offcorr':
                        mnvl.append(crmn)

                    self.open_sumdt[gainind][cr_drxn_key]['mn_' + cr_datatype_key.split('_offcorr')[0]]=crmn
                    MNCT=MNCT+1

                    self.open_sumdt[gainind][cr_drxn_key]['ster_wingdiff']=stdvl/np.sqrt(len(tmplist))

            self.open_sumdt[gainind]['summed_wing_diff']=np.sum(mnvl,axis=0)
            self.open_sumdt[gainind]['ave_time']=self.ave_time

    def get_gain(self):
        try:
            gains=self.params['gain_list'].tolist()
        except:
            gains=[self.params['closed_loop_gain']]
        if 'crgain' in self.params:
            self.crdata['crgain']=self.params['crgain']
        else:
            if self.params['closed_loop_experiment_list'][0]=='zero':
                self.crdata['crgain']=[0]+gains
            else:
                self.crdata['crgain']=gains+[0]

###
###
#functions for plotting
###
###

#called by make_plots_by_file
#crdata comes from prep_data_fxn




    def make_example_figure(self):

        if self.PLOTFLAG:
            fig=pylab.figure()
            axsubpol=[]
            self.sumax=[]
            axpol_hist={}
            ax_traj_plot=[]
            ax_displacement=[]
            gs = gridspec.GridSpec(19, 13)
            axmotor=[]
            axhist=[]
            axwing=[]
            #for key in ['stripe','sun']:

                #axpol_hist[key]=[]

            axwing.append(pylab.subplot(gs[0:4, 0:5]))
            axwing.append(pylab.subplot(gs[0:4, 6:11]))
            ax_mean_coherence=pylab.subplot(gs[5:6, 0:8])
            ax_vec_strength=pylab.subplot(gs[7:8,0:8])
            #axcolorbar=pylab.subplot(gs[5:6,9:11])
            #for key in ['stripe','exp']:
            #axmotor['stripe'].append(pylab.subplot(gs[14:18,0:5]))
            #axmotor['stripe'].append(pylab.subplot(gs[14:18,6:11]))

            #axmotor['sun'].append(pylab.subplot(gs[8:12,0:5],sharex=axwing[0]))
            #axmotor['sun'].append(pylab.subplot(gs[8:12,6:11],sharex=axwing[1]))
            axmotor=pylab.subplot(gs[14:18,0:5])

            try:
                for cr_fltnum in self.crdt.keys():

                    axhist.append(pylab.subplot(gs[14:18,5+cr_fltnum]))
            except:
                print('example error')

            fpb.make_raw_plot(self.crdt,axmotor, axhist)
            #for crind in np.arange(len(self.datvls['splitstats'])):
             #   axsubpol.append(pylab.subplot(gs[14:18,crind*2:crind*2+2], polar=True))

            #axgain=pylab.subplot(gs[14:18,5:7])
            axpos=[]
            #axpos.append(pylab.subplot(gs[14:18,7:9]))
            #axpos.append(pylab.subplot(gs[14:18,9:11]))
            axpolarizer=pylab.subplot(gs[0:4,9:11])

            #axpol_hist['stripe'].append(pylab.subplot(gs[14:18,5]))
            #axpol_hist['stripe'].append(pylab.subplot(gs[14:18,11]))
            #axpol_hist['sun'].append(pylab.subplot(gs[8:12,5]))
            #axpol_hist['sun'].append(pylab.subplot(gs[8:12,11]))


            fig_traj=pylab.figure()
            gs2 = gridspec.GridSpec(5, 2)

            ax_traj_plot.append(pylab.subplot(gs2[2,1],polar=True))

            ax_traj_plot.append(pylab.subplot(gs2[3,1],polar=True))


            try:

                crdt=self.crdt
            except:
                crdt=[]

                print ('no self.crdt!')
                return

            #plotsum=plot_sum16.Plot_Sum(example_dt=crdt,splitdt=self.datvls['splitstats'],params=self.params)
            #crdt['tim']

            if crdt:


                if TRUNC_DATA_FLAG==0:
                    self.save_file_name=self.fname.strip('.txt')+'_' + self.anal_type + '_sumdata.pck'

                if os.path.isfile(self.save_file_name):

                    os.remove(self.save_file_name)
                #if os.path.isfile(self.fname.strip('.txt')+'_' + 'UNPAIRED_sumdata.pck'):
                #   os.remove(self.fname.strip('.txt')+'_' + 'UNPAIRED_sumdata.pck')
                try:
                    self.figname=self.pckname.split('/')[-1].strip('_sumdata.pck')
                    self.fig_trajname=self.pckname.split('/')[-1].strip('_sumdata.pck')+'traj'
                except:
                    tst=1

                try:
                    
                    fh.save_to_file(self.save_file_name,self.datvls,fig=fig,fig_traj=fig_traj,figname=self.figname,fig_trajname=self.fig_trajname,save_executable=self.exec_file)
                except:
                    tst=1


            fig.savefig(self.datapath+'gains'+self.crfile.split('_')[-1].split('.txt.')[0]+'.pdf')
        else:
            tst=1



    def calc_displacement_trajectory(self):
        traj_position={}

        for i, calc_mode in enumerate (['raw']):
            traj_position[calc_mode]=[]


            crdt=self.sumstats['sun'][self.expkey]

                    #self.vec_inds=np.intersect1d(np.where(np.array(crdt['vec_time_lst'])>self.bndtimes[0][0]/60.)[0], np.where(np.array(crdt['vec_time_lst'])<self.bndtimes[0][1]/60.)[0])
            heading_lst=crdt['mn_vector_lst']
            traj_position[calc_mode]=calc.make_forward_trajectory(heading_lst,calc_type=calc_mode)


        return traj_position
    def plot_displacement_trajectory(self,axin):

        MODE='polar'
        #forward_position_list=calc.make_forward_trajectory(heading_lst,calc_mod_90=1)

        for i,calc_mode in enumerate(['raw', 'doubled']):
            ax=axin[i]
            try:

                crdt=self.sumstats['exp']['displacement_traj'][calc_mode]
            except:
                return
            if MODE is 'linear':
                #ax=pylab.subplot(2,1,i)
                crax.plot(crdt['x'],crdt['y'])
                ax.set_ylim([-500,500])
                ax.set_xlim([-500,500])
                ax.set_aspect('equal')
            else:
                forward_position_polar=calc.linear_to_polar(crdt)
                mid=np.floor(len(forward_position_polar['theta'])/2.0)

                ax.plot(forward_position_polar['theta'][0:mid],forward_position_polar['len'][0:mid],'r')
                #ax.plot(forward_position_polar['theta'][mid],forward_position_polar['len'][mid],'o','color','c')
                ax.plot(forward_position_polar['theta'][mid:],forward_position_polar['len'][mid:],'b')
                if calc_mode is 'mod_180':
                    ax.set_ylim([0,900])
                    ax.get_xaxis().set_ticks([0,np.pi/2.,np.pi,3.*(np.pi/2.)])
                    ax.get_xaxis().set_ticklabels(['0','45','90','135'],fontsize=8)
                    ax.get_yaxis().set_ticks([200,400,600,800])
                    ax.get_yaxis().set_ticklabels(['200','400','600','800'],fontsize=8)
                if calc_mode is 'raw':
                    ax.set_ylim([0,300])
                    ax.get_yaxis().set_ticks([100,200])
                    ax.get_yaxis().set_ticklabels(['100','200'],fontsize=8)

    def get_heading_vls(self):

        inds=np.intersect1d(np.where(np.array(self.crdt['vec_time_lst'])>self.bndtimes[0][0]/60.)[0], np.where(np.array(self.crdt['vec_time_lst'])<self.bndtimes[0][1]/60.)[0])
        return np.array(self.crdt['mn_vector_lst'])[inds],np.array(self.crdt['vec_time_lst'])[inds]
    def make_suptitle(self):
        tit_str=self.fname
        return tit_str

    def plot_gain_value(self,ax,axpos,**kwargs):
        gainfit=[]
        CANONICAL_PANEL_GAIN=40.0

        if 'plot_extra_fig' in kwargs:
            if kwargs['plot_extra_fig']:
                PLOT_EXTRA_FLAG=1
            else:
                PLOT_EXTRA_FLAG=0
        else:
            PLOT_EXTRA_FLAG=0

        wd_deg=self.crdt['lftwng']-self.crdt['rtwng']
        motvel_steps=self.crdt['mot_vel']
        steps_per_rotation=self.params['steps_per_rotation']
       #call standard panel gain 40.0 pixels/radian
        panel_gain_in_degrees_per_radian_wing=(CANONICAL_PANEL_GAIN/64)*360
        degrees_per_radian=(360/(2*np.pi))
        plt_panel_gain=panel_gain_in_degrees_per_radian_wing/degrees_per_radian
        try:
            motor_gain_degrees_per_radian=self.params['closed_loop_gain'][0]
            plt_motor_gain=-motor_gain_degrees_per_radian/degrees_per_radian
            mot_positions=self.crdt['mot_deg']
            self.make_position_plot(mot_positions,wd_deg,motvel_steps,steps_per_rotation,axpos)
        except:
            print('no closed_loop_gain')

        if 'cr_coeff_list' in self.params:
            motor_gain_degrees_per_radian=-150.0
            plt_motor_gain=-motor_gain_degrees_per_radian/degrees_per_radian
            #coeff=self.params['cr_coeff_list'][0]
            coeff=10.0
            #asy_vl=self.params['asy_vl_list'][0]

            asy_vl=5000.0
            linear_gains_to_plot=[plt_panel_gain/360.0,plt_motor_gain/360.0]

            fit_slope=plt.plot_wing_diff_rotation_speed(ax,wd_deg,-motvel_steps,steps_per_rotation,downsample=40,linear_gain=linear_gains_to_plot,nonlinear_gain=[coeff,asy_vl],calc_fit=1,wd_cutoff=10)
        else:
            try:
                motor_gain_degrees_per_radian=self.params['closed_loop_gain'][0]
                plt_motor_gain=-motor_gain_degrees_per_radian/degrees_per_radian
                linear_gains_to_plot=[plt_panel_gain/360.0,plt_motor_gain/360.0]

                fit_slope=plt.plot_wing_diff_rotation_speed(ax,wd_deg,-motvel_steps,steps_per_rotation,downsample=2,linear_gain=linear_gains_to_plot,calc_fit=1,wd_cutoff=10)
            except:
                fit_slope=0
                print('no closed_loop_gain')
        if PLOT_EXTRA_FLAG:
            fig=pylab.figure()
            ax=[]


            #positions are 45 degrees, 90 degrees, 135 degrees
            positions_to_plot=[[[40,50],[220,230]],[[87.5,92.5],[267.5,272.5]],[[60,65],[240,245]],[[310,320],[130,140]],[[77.5,80],[255.5,260]]]
            ax.append(fig.add_subplot(2,1,1))
            colvls=['r','c','k','g','m']
            for ind,crpos_list_to_plot in enumerate(positions_to_plot):
                crposinds=[]
                for crpos_to_plot in crpos_list_to_plot:

                    cr_mot_positions=mot_positions[np.arange(0,len(mot_positions)/2)]
                    crposinds.append(np.intersect1d(np.where(cr_mot_positions>crpos_to_plot[0])[0],np.where(cr_mot_positions<crpos_to_plot[1])[0]))
                cr_comb_inds=np.concatenate((crposinds[0],crposinds[1]),axis=2)

                try:
                    plt.plot_wing_diff_rotation_speed(ax[-1],wd_deg[cr_comb_inds],-motvel_steps[cr_comb_inds],steps_per_rotation,linear_gain=linear_gains_to_plot,wd_cutoff=20,color_flag=colvls[ind])
                except:
                    print('no gain plot')
                ax[-1].set_aspect(50)
            fig.savefig(self.datapath+'gains'+self.crfile.split('_')[-1].split('.txt.')[0]+'.pdf')

        return fit_slope


    def make_position_plot(self,mot_deg,wd_deg,motvel_steps,steps_per_rotation,axpos):
        motvel_rotations=(motvel_steps/steps_per_rotation)
        sy, bin_edges = np.histogram(mot_deg, bins=100, weights=np.abs(wd_deg))
        n, bin_edges = np.histogram(mot_deg, bins=100)
        sy2, bin_edges = np.histogram(mot_deg, bins=100, weights=np.abs(wd_deg)*np.abs(wd_deg))
        mean = sy / n
        ster = (np.sqrt(sy2/n - mean*mean))/np.sqrt(n)

        axpos[-2].errorbar((bin_edges[1:] + bin_edges[:-1])/2, mean, yerr=ster, fmt='r-')
        #axpos[-2].plot(self.crdt['mot_deg'],np.abs(wd_deg),'k.')
        #axpos[-1].plot(self.crdt['mot_deg'],np.abs(motvel_rotations),'k.')
    def plot_coherence_with_mean(self,ax):
        ax.plot(self.crdt['vec_time_lst'],self.crdt['mn_dot_product_lst'])


    def plot_raw_vector_strength(self,ax):

        ax.plot(self.crdt['vec_time_lst']-self.crdt['time_in_min'][0],self.crdt['len_vector_lst'])

        ax.set_xlim([0,15])
        ax.set_ylim([-.1,1.1])
        ax.get_yaxis().set_ticks([0,0.25,0.5,0.75,1.0])
        ax.get_yaxis().set_ticklabels(['0','0.25','0.75','1.0'],fontsize=6)

    def make_split_data(self,cr_mn,cr_ln,bin_vls,range_vls):
        #HARD CODE FOR 7.5 MINUTES
        if self.variable_gain_flag:
            try:
                HALF_SECONDS=self.params['variable_closed_loop_duration_sec']
            except:
                HALF_SECONDS=300.
            half_point=HALF_SECONDS/STEP_LEN_IN_SEC
            first_half_inds=np.arange(int(half_point))
        else:
            HALF_SECONDS=7.5*60
            half_point=HALF_SECONDS/STEP_LEN_IN_SEC
            first_half_inds=np.arange(int(half_point))


        if self.variable_gain_flag:
            splitmn=[]
            splitln=[]
            for crinds in [0,1]:
                splitmn.append(np.array(self.crdt['mn_vector_lst']))
                splitln.append(np.array(self.crdt['len_vector_lst']))
            return splitmn, splitln
        else:
            if len(cr_mn)>first_half_inds[-1]:

                second_half_inds=np.arange(first_half_inds[-1],len(cr_mn),1)
                splitmn=[]
                splitln=[]

                for crinds in [first_half_inds,second_half_inds]:
                    splitmn.append(np.array(self.crdt['mn_vector_lst'])[crinds])
                    splitln.append(np.array(self.crdt['len_vector_lst'])[crinds])

                return splitmn,splitln
            else:
                splitmn=[]
                splitln=[]
                return splitmn,splitln
        #for crind in [0,1]:
         #   if crind is 0:

          #  else:





    def calc_vec_strength_threshold_mean(self,input_ln,input_mn):

        hstvls={}
            #ax_traj.plot(cr_mn,cr_ln)

        threshinds=np.where(np.array(input_ln)>VEC_STRENGTH_THRESHOLD)[0]
        out_dt=calc.make_hist_calculations(calc.rad_to_deg(input_mn[threshinds]),self.degbins)
        for key in ['normhst','mnrad','circvar']:
            hstvls[key]=out_dt[key]

            #hstvls['xrad']=self.xrad
            #hstvls['rad_per_bin']=self.rad_per_bin
        return hstvls

    def norm_heat_map(self,heatmap):

        summed_heatmap=sum(sum(heatmap))
        heatmap_norm=heatmap/summed_heatmap
        return heatmap_norm

    def save_heat_map_data(self,heatmap_norm,xedges,yedges,r,theta,**kwargs):


        self.datvls['sumstats']['exp']['heat_map']={}
        self.datvls['sumstats']['exp']['heat_map']['norm_heat_map_vls']=heatmap_norm
        self.datvls['sumstats']['exp']['heat_map']['redges']=xedges
        self.datvls['sumstats']['exp']['heat_map']['thetaedges']=yedges
        self.datvls['sumstats']['exp']['heat_map']['r']=r
        self.datvls['sumstats']['exp']['heat_map']['theta']=theta
        self.datvls['sumstats']['exp']['heat_map']['realigned_norm_heat_map_vls']=realigned_heatmap_norm
        if 'vector_strength_mean' in kwargs:
            for ind in [0,1]:
                self.datvls['splitstats'][ind]['vector_strength_mean']=[]
        for ind in [0,1]:
            if split_heat_map_norm:
                if 'vector_strength_mean' in kwargs:

                    self.datvls['splitstats'][ind]['vector_strength_mean'].append(kwargs['vector_strength_mean'][ind])

                self.datvls['splitstats'][ind]['norm_heat_map_vls']=split_heat_map_norm[ind]

                self.datvls['splitstats'][ind]['realigned_norm_heat_map_vls']=realigned_split_heat_map_norm[ind]

        return self.datvls['sumstats']['exp']['heat_map']
    def plot_combdata(self):
        #make wingax,polarizerax,and sumax
        self.axhandles=[]
        if self.exp_type=='closed':
            if len(self.datvls.keys())>1:
                self.plot_closed_stroke_data()
                
                self.plot_closed_sumdata()
        elif self.exp_type=='open':
            self.plot_open_stroke_data()
            #self.plot_cut_data()

    #rewritten to use data from self.datvls
    #quickly hacked to automatically plot paired data
    #will need to modify in future to plot single
    def plot_closed_stroke_data(self):
        nrow=7
        ncol=1
        OFFSET_TIME=0
        self.axhandles,self.fig=set_up_closed_loop_figure(nrow,ncol)
        self.fig.suptitle(self.fname)
        if self.anal_type=='paired':
            for self.pairind in [0,1]:
                self.make_plots(OFFSET_TIME)
                OFFSET_TIME=OFFSET_TIME+20.0
        else:
            self.make_plots(OFFSET_TIME)

    def make_plots(self,OFFSET_TIME):
        MAXDEGPLOT=358

        PLOTFLAG=1
        MINDEGPLOT=2
        hstvls={}
        CONTROL_LAG_TIME=30
        if self.anal_type=='paired':
            crdt=self.datvls[self.pairind]['sumstats']['exp']
        else:
            try:
                crdt=self.datvls['sumstats']['exp']
            except:
                PLOTFLAG=0
        if PLOTFLAG:
            plot_time=OFFSET_TIME+crdt['time_in_min']
            plotinds1=np.where(crdt['mot_deg']<MAXDEGPLOT)
            plotinds2=np.where(crdt['mot_deg']>MINDEGPLOT)
            allinds=np.intersect1d(np.array(plotinds1[0]),np.array(plotinds2[0]))
            lftwng=crdt['lftwng']
            rtwng=crdt['rtwng']
            pol_sensor=crdt['pol_sensor']
            minbnds=np.where(pol_sensor>POL_SENSOR_MINIMUM)
            srtind=np.argsort(crdt['mot_deg'][minbnds])
            wingdiff=np.array(lftwng)-np.array(rtwng)

            MAXTIME= np.max(plot_time)
            self.axhandles[0].plot(plot_time,np.array(lftwng),'r',linewidth=0.4)
            self.axhandles[0].plot(plot_time,np.array(rtwng),'c',linewidth=0.4)
            self.axhandles[2].plot(crdt['vec_time_lst'],crdt['len_vector_lst'],'k')
            self.axhandles[3].plot(plot_time,pol_sensor,'k')


            #should change this to be vec_axh=self.axhandles[4]
            #fig2=pylab.figure()
            self.vecaxh=fig2.add_subplot(2,1,1,polar=True)
            self.vecaxh.plot(cr_mn,cr_ln)
            threshinds=np.where(np.array(cr_ln)>THRESHOLD)[0]

            tmphst,mnrad,tmpvar=calc.make_hist_calculations(calc.rad_to_deg(cr_mn[threshinds]),self.degbins)
            self.vechist_axh=fig2.add_subplot(2,1,2,polar=True)
            hstvls['normhst']=tmphst
            hstvls['mnrad']=mnrad
            hstvls['xrad']=self.xrad
            hstvls['rad_per_bin']=self.rad_per_bin
            
            plt.polar_plot(self.vechist_axh,hstvls)
            if self.anal_type=='paired':
                self.datvls[self.pairind]['sumstats']['exp']['vecmn']=mnrad
                self.datvls[self.pairind]['sumstats']['exp']['vecvar']=tmpvar
            else:
                self.datvls['sumstats']['exp']['vecmn']=mnrad
                self.datvls['sumstats']['exp']['vechst']=tmphst
                self.datvls['sumstats']['exp']['vecvar']=tmpvar


            splitinds=np.array_split(threshinds,np.array(np.where(np.diff(threshinds)!=1))[0]+1)
            for rotnum,crsplitinds in enumerate(splitinds):
                if np.size(crsplitinds):
                    self.vecaxh.plot(cr_mn[crsplitinds],cr_ln[crsplitinds],'r')



            self.vecaxh.set_ylim([0,1.2])

            if PLOT_SCATTER_POL:
                self.axhandles[4].plot(crdt['mot_deg'][minbnds][srtind],pol_sensor[minbnds][srtind])

                titlestr='period=%f/nentropy=%f'%(self.calculated_motor_period,crdt['entropy'])
                self.axhandles[4].set_title(titlestr,fontsize=8)
            splitinds=np.array_split(allinds,np.array(np.where(np.diff(allinds)!=1))[0]+1)
            for rotnum,crsplitinds in enumerate(splitinds):
                if np.size(crsplitinds):
                    self.axhandles[1].plot(plot_time[crsplitinds],crdt['mot_deg'][crsplitinds],'k')
            self.plot_boundary_lines()



#called by plot_combdata
    #plots polardata
    def plot_closed_sumdata(self):
        self.startrow=5
        self.numcol=6
        self.startcol=1
        if self.anal_type=='paired':
            for self.pairind in [0,1]:
                self.make_sum_plots()
        else:
            self.make_sum_plots()

    def make_sum_plots(self):
        PLOTFLAG=1
        if self.anal_type=='paired':
            crsumdt=self.datvls[self.pairind]['sumstats']['exp']
            crsplitdt=self.datvls[self.pairind]['splitstats']
        else:
            try:
                crsumdt=self.datvls['sumstats']['exp']
                crsplitdt=self.datvls['splitstats']
            except:
                PLOTFLAG=0

        if PLOTFLAG:
            self.vecaxh.plot([crsumdt['vecmn'],crsumdt['vecmn']+np.pi],[1.1,1.1],'k',linewidth=2)
            self.vecaxh.plot([crsumdt['mnrad'],crsumdt['mnrad']+np.pi],[1.1,1.1],'c',linewidth=2)
            self.axhandles.append(self.fig.add_subplot(self.startrow+1,self.numcol,self.startrow*self.numcol+self.startcol,polar=True))
            plt.polar_plot(self.axhandles[-1],crsumdt)
            self.startcol=self.startcol+1
            for crkey in crsplitdt:
                plot_dt=crsplitdt[crkey]
                self.axhandles.append(self.fig.add_subplot(self.startrow+1,self.numcol,self.startrow*self.numcol+self.startcol,polar=True))
                plt.polar_plot(self.axhandles[-1],plot_dt,add_text=1,title='SUB'+str(crkey))
                self.startcol=self.startcol+1

#first pass, just plot all the data
    def plot_open_stroke_data(self):
        plotkeys=['uncorrected']
        nrows=5
        ncol=2
        self.fig=[]
        for self.plot_type in plotkeys:
            for self.crgainvl,cr_gain_data in enumerate(self.open_sumdt):
                axh,fig=set_up_open_loop_figure(nrows,ncol)
                fig.suptitle(self.fname+'\n'+self.plot_type)
                craxvl=0
                keylist=['pos','neg']
                for crkey in keylist:
                    self.crdt=cr_gain_data[crkey]

                    cr_gain_vl=self.crdt[0]['open_loop_gain']
                    self.crax=axh[craxvl:craxvl+nrows]
                    self.crax[0].set_title('gain is %f'%cr_gain_vl)
                    self.raw_height=0
                    craxvl=craxvl+nrows
                    indlist=[x for x in self.crdt.keys() if type(x)  is int]
                    for self.crexpnum in indlist:
                        self.make_open_loop_plot()
                        self.raw_height=self.raw_height+50
                    self.calc_mean_wingdiff(crkey)
                    self.plot_mean_vls(crkey)
                self.calc_mean_wingdiff(keylist)
                self.plot_mean_vls(keylist)
                self.fig.append(fig)

    def calc_mean_wingdiff(self,inkey):
        tmplist=[]
        if type(inkey) is str:

            indlist=[x for x in self.crdt.keys() if type(x)  is int]
            for crnum in indlist:
                tmplist.append(self.crdt[crnum]['interp_wingdiff'])
            mdat = np.ma.masked_array(tmplist,np.isnan(tmplist))
            mm = np.mean(mdat,axis=0)
            stm=np.std(mdat,axis=0)
            self.open_sumdt[self.crgainvl][inkey]['wingdiff_mean']=mm
            self.open_sumdt[self.crgainvl][inkey]['wingdiff_std']=stm
        elif type(inkey) is list:
            for key in inkey:
                tmplist.append(self.open_sumdt[self.crgainvl][key]['wingdiff_mean'])

            self.open_sumdt[self.crgainvl]['summed_wingdiff']=np.sum(tmplist,axis=0)

    def plot_cut_data(self):
        for gainind,cr_gain_data in enumerate(self.open_sumdt):
            AXIND=0
            axh,fig=set_up_open_loop_figure(3,2,ax_order='normal')
            mnvl=[]
            for cr_drxn_key in ['pos','neg']:
                crdt=self.open_sumdt[gainind][cr_drxn_key]
                for crdatkey in ['mn_pol_sensor_cut', 'mn_wingdiff_cut']:
                    for crind in crdt[crdatkey].keys():
                        crvls=crdt[crdatkey][crind]
                        if crind=='mean':
                            col='r'
                            if crdatkey=='mn_wingdiff_cut':
                                mnvl.append(crvls)
                        else:
                            col='k'

                        axh[AXIND].plot(crvls,color=col)
                    AXIND=AXIND+1
            sumout=np.sum(mnvl,axis=0)
            for crmn in mnvl:
                axh[AXIND].plot(crmn,color='k')
            axh[AXIND].plot(sumout,color='r')

    #called by plot_open_stroke_data
    def make_open_loop_plot(self):
        crexp=self.crdt[self.crexpnum]
        timevls=crexp['net_time']
        #if self.plot_type == 'corrected':
        #    suffix='_offcorr'
        #elif self.plot_type == 'uncorrected':
        #    suffix=''
        self.crax[0].plot(timevls,crexp['rtwng']+self.raw_height,'r')
        #if not len(timevls) == len(crexp['wingdiff'+suffix]):
        self.crax[0].plot(timevls,crexp['lftwng']+self.raw_height,'c')
        self.crax[1].plot(timevls,crexp['wingdiff'],'k')
        self.crax[2].plot(timevls,crexp['pol_sensor'])

        mot=crexp['motor_position_radians']
        self.mot_binvls=np.linspace(0,2*np.pi,num=100)
        srt_ind=np.argsort(mot)
        interp_pol_sensor=np.interp(self.mot_binvls,mot[srt_ind],crexp['pol_sensor'][srt_ind])
        self.crax[3].plot(self.mot_binvls,interp_pol_sensor)
        interp_wing_diff=np.interp(self.mot_binvls,mot[srt_ind],crexp['wingdiff'][srt_ind])
        self.crdt[self.crexpnum]['interp_wingdiff']=interp_wing_diff

        self.crax[4].plot(self.mot_binvls,interp_wing_diff,'k')
        #self.crax[3].plot(uncorr_timevls,crexp['raw_motor_position'],'k')
        #fourth plot is radians vs. polarizer value
        #fifth plot is radians vs. wingbeat difference

    def plot_mean_vls(self,inkey):

        if type(inkey) is str:
            crdt=self.open_sumdt[self.crgainvl][inkey]['wingdiff_mean']

            self.crax[4].plot(self.mot_binvls,crdt,'r')
        elif type(inkey) is list:
            self.crax[4].plot(self.mot_binvls,self.open_sumdt[self.crgainvl]['summed_wingdiff'],'c',linewidth=2)

    #returns Hour-Minute int for time
    def get_time(self):
        return()

    def plot_boundary_lines(self):
        axh=self.axhandles[1]
        yvls=[-10,370]
        for i,vls in enumerate(self.bndtimes):
            if self.exp_type=='closed':
                for edges in vls:
                    axh.plot([edges/60,edges/60],[yvls[0], yvls[1]],'r--')
            #add text with current gain

                    axh.text(np.mean(vls)/60,yvls[1],self.params['gain_list'])
            else:
                axh.plot([vls/60,vls/60],[yvls[0], yvls[1]],'r--')
    #make new file (that can be deleted)

    def check_calibration(self):
        SAMPLE_SPACING=10
        import peakdetect
        period_data={}
        #calc_calibration_file=[x for x in os.listdir(datapath) if 'period.pck' in x]

        #filevl=[x for x in os.listdir(datapath) if 'calibration.pck' in x]
        calibrate_file=[x for x in os.listdir(self.datapath) if 'calculated_period.pck' in x]
        try:
            #
            period_data=fh.open_pickle(self.datapath+calibrate_file[0])

            self.calculated_motor_period=period_data['calculated_motor_period']
        except:
            crbnds=self.get_index_bnds(self.indt['net_time'],0)
            mot=self.indt['raw_motor_position'][crbnds]
            mot_x=np.arange(np.min(mot),np.max(mot),SAMPLE_SPACING)
            pol=self.indt['pol_sensor'][crbnds]
            output_vls,period,peakdata=calc.find_period(mot_x,mot,pol)
            self.calculated_motor_period=SAMPLE_SPACING*period
            period_data['calculated_motor_period']=SAMPLE_SPACING*period
            #pdb.set_trace()
            fh.save_to_pickle(self.datapath+'calculated_period.pck',period_data)
################



def check_if_period_calculated(fname):
    #check if file exists
    pck_fname=fname.strip('.txt')+'_sumdata.pck'
    if os.path.isfile(pck_fname):
        params=fh.open_pickle(pck_fname)
        try:
            return params['steps_per_rotation'],params['corr_dict']
        except:
            return
    else:
        return

def set_up_open_loop_figure(nrows,ncol,**kwargs):
    if 'ax_order' in kwargs:
        ax_order=kwargs['ax_order']
    else:
        ax_order='by_column'

    fig=pylab.figure()

    ax=[]
    totalaxes=nrows*ncol
    if ax_order=='by_column':
        startvls=[1,2]
        for cr_start_vl in startvls:
            for sub_axis_num in np.arange(0,totalaxes/2):
                print cr_start_vl+sub_axis_num*2
                ax.append(fig.add_subplot(nrows,ncol,cr_start_vl+sub_axis_num*2))
    else:
        for sub_axis_num in np.arange(1,totalaxes+1):
            ax.append(fig.add_subplot(nrows,ncol,sub_axis_num))

    
    return ax,fig

def set_up_closed_loop_figure(nrow,ncol):
    fig=pylab.figure()
    ax=[]
    totalaxes=nrow*ncol

    for craxnum in np.arange(1,totalaxes+1):
        if craxnum<(totalaxes-1):
            ax.append(fig.add_subplot(nrow,ncol,craxnum))
        #else:
            #ax.append(fig.add_subplot(nrow,ncol,craxnum,polar=True))
    return ax,fig


if __name__ == '__main__':
    import sys
    sun_anal = Sun_Anal(sys.argv)
    params=sun_anal.run_analysis()


#################
#################
