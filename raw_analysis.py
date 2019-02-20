#!/usr/bin/python

#raw_analysis.py performs initial analysis of raw
#data collected in flight behavior rig during sun experiments

###USAGE 
#python raw_analysis.py raw_date_range.txt

###Dependencies
#various libraries from https://github.com/timothylwarren/py_utilities

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


input_path='/users/tim/Dropbox/flystuff/uo_data/'
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
          

    def get_data_by_file(self):
        mod_file_name=os.path.join(self.in_data_path,'modified_raw_data.pck')
        self.fname = os.path.join(self.in_data_path,self.crfile)
        pck_fname=self.fname.strip('.txt')+'.pck'
        try:
            self.params=fh.open_pickle(pck_fname)
        except:
            print('open pickle error')
       
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

           

    def add_params(self):

        self.params['filter_len_in_sec']=FILTER_LEN_IN_SEC
        self.params['step_len_in_sec']=STEP_LEN_IN_SEC
        self.params['vec_strength_threshold']=VEC_STRENGTH_THRESHOLD
    
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
           
        self.sumstats[self.crkey]['lftwng']=np.array(self.crdata['lftwng'])[self.analinds]
        self.sumstats[self.crkey]['rtwng']=np.array(self.crdata['rtwng'])[self.analinds]

        self.sumstats[self.crkey]['mot_position_in_rad']=np.array(self.crdata['adjusted_raw_motor_position_in_rad'])[self.analinds]

        if 'gain_coefficient' in self.crdata.keys():

            self.sumstats[self.crkey]['gain_coefficient']=np.array(self.crdata['gain_coefficient'])[self.analinds]

        self.sumstats[self.crkey]['mot_deg']=self.outdict['mot_deg'][self.analinds]

        self.sumstats[self.crkey]['mot_rad']=np.array(self.outdict['mot_rad'])[self.analinds]
        try:

            self.sumstats[self.crkey]['cloop_duration_in_min']=np.ceil(self.sumstats[self.skey]['time_in_min'][-1]-self.sumstats[self.skey]['time_in_min'][0])


           
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
#functions for plotting
###
###


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
           

            axwing.append(pylab.subplot(gs[0:4, 0:5]))
            axwing.append(pylab.subplot(gs[0:4, 6:11]))
            ax_mean_coherence=pylab.subplot(gs[5:6, 0:8])
            ax_vec_strength=pylab.subplot(gs[7:8,0:8])
           
            axmotor=pylab.subplot(gs[14:18,0:5])

            try:
                for cr_fltnum in self.crdt.keys():

                    axhist.append(pylab.subplot(gs[14:18,5+cr_fltnum]))
            except:
                print('example error')

            fpb.make_raw_plot(self.crdt,axmotor, axhist)
           
            axpos=[]
            
            axpolarizer=pylab.subplot(gs[0:4,9:11])

           

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

           

            if crdt:


                if TRUNC_DATA_FLAG==0:
                    self.save_file_name=self.fname.strip('.txt')+'_' + self.anal_type + '_sumdata.pck'

                if os.path.isfile(self.save_file_name):

                    os.remove(self.save_file_name)
               
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
