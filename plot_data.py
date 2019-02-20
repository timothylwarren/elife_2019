#plot_data.py makes plots that comprise Figures 7, Figure 7 supplementary,
#and Figure 8 supplementary of Sullivan, Warren, Doe 2019, submitted to Elife

#USAGE 
#python plot_data.py 

#dependencies: figurefirst, 
#            various libraries from https://github.com/timothylwarren/py_utilities

##FLY_TYPE_TO_PLOT sets which data will be plotted. Only one FLY_TYPE_TO_PLOT value
#below should be uncommented. Other two commented out.

#Fig. 7 supplement A-D 
# ss00096Gal4/UAS-Kir2.1,empty split Gal4/UAS-Kir2.1
#FLY_TYPE_TO_PLOT=['ss96_uaskirj','es_uaskirj']

#Fig. 8 supplement 
# empty split GAl4/UAS Toy RNAi, ss00096Galr/UAS Toy RNAi
#FLY_TYPE_TO_PLOT=['estoy_lc_flies','ss96toy_lc_flies']

#Fig. 7 and Fig. 7 supplement F-I
# mcherryRNAI with heat shock, eyeless RNAi no heat shock, eyeless RNAi heat shock
FLY_TYPE_TO_PLOT=['mcherryhs_flies','ey_rnai_con_flies','eyrnai_exp_flies']


import numpy as np
import os
import pylab
import pdb
import pickle
from py_utilities import circ
import shutil
import scipy.stats.morestats as st
import scipy.stats as st_orig
from py_utilities import tw_plot_library3 as twplt
import datetime
from py_utilities import tw_filehandling as fh
from py_utilities import fp_library2 as fpl
from py_utilities import tw_calc_library as calc
from py_utilities import  fly_plot_basics as fpb
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt_ax
from scipy import stats
#from py_utilities import svg_to_axes as fifi
import figurefirst as fifi
mpl.rcParams['svg.fonttype']='none'
mpl.rcParams['font.family']='sans-serif'
mpl.rcParams['font.sans-serif']='Arial'

PLOT_DISPLACEMENT_FLAG=True
PLOT_NET_DISPLACEMENT_FLAG=True
PLOT_ALL_RAW=False
EXAMPLE_FLIES_TO_PLOT=[1,1,1]
EXAMPLE_FLIGHTS_TO_PLOT=[1,1,0]
PLOT_CORRELATION=True
PLOT_EX_TEXT_STR=True
#COMPARE_VECTOR_STRENGTH=True
CORR_PLOT_COLS=[0,1]
center_on_zero_flag=True

WITHHOLD_LEGEND_POLAR=True
COMBINE_FIRST_AND_SECOND_DISPLACEMENT_FLAG=True
PLOT_SCATTER=False
PLOT_VEC_STRENGTH=True
PLOT_HEADING_HIST=False
PLOT_MEAN_DISPLACEMENT=True
PERMUTE_BY_INDIVIDUALS_NOT_FLIGHTS=True
#2017-6-23 fly6
EXFLY=68
NROW=8
MAXROW=13
NUM_EXAMPLE_FIGS=1
NHISTCOL=3
CORR_PLOT_COLS=[0,1]
TIME_GAP=5
FRONTAL_THRESH=np.pi/4

dropbox_path='/users/tim/Dropbox/'
pylab.ion()

sum_layout_svg=dropbox_path+'generic_summary_layout.svg'
raw_plot_svg=dropbox_path+'generic_example_layoutb.svg'
input_data_path=dropbox_path+'flystuff/uo_data/'
sum_data_path=dropbox_path+'flystuff/summary_data/'
save_fig_folder='/raw_figs/'



COLVLS=['k','c','r']
VEC_STRENGTH_THRESH=0.2
PERMUTE_FLAG=True

#FLIES_TO_PLOT=[0,1,2]
class Example_Plot():
    #calls prep_data_fxn for each file
    def __init__ (self,**kwargs):
        self.datapath=sum_data_path
        self.pck_open_name='vecstrength_unpaired_sumsundata.pck'

        self.sumdt=fh.open_pickle(self.datapath+self.pck_open_name)['sumdt']
        
        self.adt=fh.open_pickle(self.datapath+self.pck_open_name)['adt']
        pdb.set_trace()



    def make_plots(self):
        self.ordered_fly_indices={}
        self.axraw={}
        svg_fig=raw_plot_svg

        self.layout={}
        for crflytype in FLY_TYPE_TO_PLOT:
            self.layout[crflytype]=[]
            for crind in np.arange(NUM_EXAMPLE_FIGS):
                self.layout[crflytype].append(fifi.FigureLayout(svg_fig))
                self.layout[crflytype][crind].make_mplfigures()
        for crflytype in FLY_TYPE_TO_PLOT:
            self.axraw[crflytype]=self.set_up_plot_ax(crflytype)
        self.sum_layout={}
        self.axsum={}
        for crflytype in FLY_TYPE_TO_PLOT:
            self.sum_layout[crflytype]=fifi.FigureLayout(sum_layout_svg)
            self.sum_layout[crflytype].make_mplfigures()
            self.axsum[crflytype]=self.set_up_sum_plot_ax(self.sum_layout[crflytype])

        if PLOT_ALL_RAW:
            for crtype in FLY_TYPE_TO_PLOT:

                self.plot_all_raw(crtype)
                self.ordered_fly_indices[crtype]=self.ordered_indlist

        #self.sumfig=pylab.figure()



        #if PLOT_EXAMPLE:
         #   self.plot_example(ax)
        self.plot_summary()
        if PLOT_SCATTER:
            self.axscatter['gal_tnt_inact_flies']=ax['tnt_inact_scatter']
            self.axscatter['wedge_kir_flies']=ax['wedge_scatter']
            for crtype in ['gal_tnt_inact_flies','wedge_kir_flies']:
                if crtype is 'gal_tnt_inact_flies':
                    crax=ax['tnt_inact_scatter']
                else:
                    crax=ax['wedge_scatter']

                self.make_plot_of_scatter_headings(crax,crtype)
        #self.sumfig.savefig(self.datapath+'sumfig.pdf')

        if PLOT_CORRELATION:
            ax=self.axsum
            self.plot_correlation(ax)



        self.savepath=input_data_path

        for crflytype in FLY_TYPE_TO_PLOT:
            for crind in np.arange(NUM_EXAMPLE_FIGS):
                add_str='sun_dataraw_%s%s.svg'%(str(crflytype),str(crind))
                save_svg=self.savepath+add_str
            #self.layout.set_layer_visability('Layer 1',False)
                #self.layout[crind].insert_figures()
                self.layout[crflytype][crind].save(save_svg)
        for crflytype in FLY_TYPE_TO_PLOT:
            save_sum_svg=self.savepath+str(crflytype)+'sun_sumdata.svg'
            #self.sum_layout[crfltype].insert_figures()
            self.sum_layout[crflytype].save(save_sum_svg)

        self.vecfig.savefig(dropbox_path+'vec_strength.eps')


    def plot_example(self,ax):
        self.get_data(EXFLY,0)
        axmot={}


        if self.params['stripe_first']:

            axmot['stripe']=[ax['ex1'],ax['ex3']]
            axmot['sun']=[ax['ex2'],ax['ex4']]

        else:
            axmot['stripe']=[ax['ex2'],ax['ex4']]
            axmot['sun']=[ax['ex1'],ax['ex3']]
        left_axis_flag=True
        for crkey in ['stripe','sun']:
            for fltnum in [0,1]:
                crax=axmot[crkey][fltnum]
                fpb.plot_motor(self.crdt[crkey][fltnum],crax,plot_vector=False,plot_split=1,plot_start_angle=0,xlim=[0,5],xticks=[0,5],xticklabels=['0','5'],plot_left_axis=left_axis_flag,withhold_bottom_axis=True, linewidth=0.3,subtract_zero_time=False)
                

                if left_axis_flag:
                    crax.text(-4,100,'Intrinsic Gal4 \n TNT-inactive',fontsize=7)
                    crax.set_ylabel('stimulus ($^\circ$)',fontsize=6)

                left_axis_flag=False


    def plot_summary(self):


        
        if PLOT_VEC_STRENGTH:
            self.vecfig=pylab.figure()
            self.ax_vecstrength=[]
            self.ax_hist_dir=[]
            self.ax_nonfrontal_vec=[]
            self.ax_vecstrength.append(pylab.subplot(341))
            self.ax_vecstrength.append(pylab.subplot(342))
            self.ax_vecstrength.append(pylab.subplot(343))
            self.ax_vecstrength.append(pylab.subplot(344))
            self.ax_nonfrontal_vec.append(pylab.subplot(345))
            self.ax_nonfrontal_vec.append(pylab.subplot(346))
            self.ax_nonfrontal_vec.append(pylab.subplot(347))
            self.ax_nonfrontal_vec.append(pylab.subplot(3,4,8))
            self.ax_hist_dir.append(pylab.subplot(3,4,9))
            self.ax_hist_dir.append(pylab.subplot(3,4,10))
            self.ax_hist_dir.append(pylab.subplot(3,4,11))
            self.ax_hist_dir.append(pylab.subplot(3,4,12))

            #self.plot_hist_dir(self.ax_hist_dir)
            self.plot_vec_strength_and_dir(self.ax_vecstrength,self.ax_hist_dir,self.ax_nonfrontal_vec)
        
        self.plot_net_displacement(self.axsum)

    def plot_correlation(self,ax):

        for flytype in FLY_TYPE_TO_PLOT:
            rad_list=[]
            vec_list=[]
            vec_plot_list=[]
            crax=ax[flytype]['corr_scatter']
            craxhist=ax[flytype]['corr_hist']
            for ind in CORR_PLOT_COLS:


                cr_rad=self.adt[flytype][ind]['mnrad']
                cr_vec=self.adt[flytype][ind]['vec_strength']


                vec_strength=np.reshape(cr_vec,len(cr_vec))
                radvls=np.reshape(cr_rad,len(cr_rad))
                #plt_rad = radvls[~np.isnan(radvls)]
                plt_rad=radvls
                #plt_vec=vec_strength[~np.isnan(radvls)]
                plt_vec=vec_strength
                if center_on_zero_flag:
                    highinds=np.where(plt_rad>np.pi)[0]
                    plt_rad[highinds]=-(2*np.pi-plt_rad[highinds])
                rad_list.append(plt_rad)
                vec_list.append(plt_vec)
            threshinds=np.intersect1d(np.where(vec_list[0]>VEC_STRENGTH_THRESH)[0],np.where(vec_list[1]>VEC_STRENGTH_THRESH)[0])

            vec_plot_list.append(vec_list[0][threshinds])
            vec_plot_list.append(vec_list[1][threshinds])
            #twplt.scatterplot(crax,rad_list[0][threshinds],rad_list[1][threshinds],plot_error_bar=True,dynamic_sizes=vec_plot_list,error_scale_factor=1)
            #twplt.scatterplot(crax,rad_list[0][threshinds],rad_list[1][threshinds]+2*np.pi,plot_error_bar=True,dynamic_sizes=vec_plot_list,error_scale_factor=1)
            if PERMUTE_FLAG:
                for crtype in ['thresh','all']:
                    if crtype is 'thresh':
                        crinds=threshinds
                        crcol='c'
                        ht=.02
                        markht=.01
                    else:
                        crinds=np.arange(len(rad_list[0]))
                        crcol='k'
                        ht=.04
                        markht=0

                    try:
                        return_dict=calc.permute_diffs(np.column_stack((rad_list[0][crinds],rad_list[1][crinds])),max_diff=np.pi)
                        actual_diff=calc.calc_heading_diff(np.column_stack((rad_list[0][crinds],rad_list[1][crinds])),calc_abs_diff=True,max_diff=np.pi)

                        self.plot_permute_hist(craxhist,return_dict['permuted_dist'],color=crcol)

                        craxhist.plot(np.nanmean(actual_diff),markht,marker='v',color=crcol,clip_on=False)

                        try:
                            prctile=stats.percentileofscore(return_dict['permuted_dist'],np.nanmean(actual_diff))
                        except:
                            print('percentile error')
                        craxhist.text(np.pi/6, ht,str(format(prctile/100.,'.2f')),fontsize=6,color=crcol)
                    except:
                        print('hist plot error')
            fpl.adjust_spines(crax,['left','bottom'])
            crax.get_yaxis().set_ticks([-np.pi,0,np.pi,2*np.pi, 3*np.pi])
            crax.get_xaxis().set_ticklabels(['-180','0','180'],fontsize=6)
            crax.get_xaxis().set_ticks([-np.pi,0,np.pi])
            crax.get_yaxis().set_ticklabels(['-180','0','180','360','180'],fontsize=6)
            crax.set_xlabel('first heading',fontsize=6)
            crax.set_ylabel('second heading',fontsize=6)
            crax.set_xlim([-np.pi-.5,np.pi+.5])
            crax.set_ylim([-np.pi-.5,3*np.pi+.5])
            crax.plot([-np.pi,np.pi],[-np.pi,np.pi],'b',linewidth=0.5)
            crax.plot([-np.pi,np.pi],[np.pi,3*np.pi],'b',linewidth=0.5)
            crax.set_aspect(1)

                            #crax.plot([plt_rad, plt_rad],[np.zeros(len(plt_vec)), plt_vec],linewidth=0.2,color=[0.5,0.5,0.5],zorder=1)

                #crax.scatter(plt_rad,plt_vec, color='0.2', s=6,alpha=0.5,zorder=2)


    def plot_permute_hist(self,ax, permuted_dist,**kwargs):
        try:
            crcol=kwargs['color']
        except:
            crcol='k'
        twplt.plot_hist(ax,permuted_dist,norm=True,hst_bnds=[0,np.pi],col=crcol,num_bins=200,linewidth=0.5)

        fpl.adjust_spines(ax,['left','bottom'])
        #ax.set_aspect('equal')
        #ax.plot([0,0.8],[0,0.8],'k--')
        ax.set_ylim([0, 0.075])

        ax.set_xlim([calc.deg_to_rad(0), calc.deg_to_rad(180)])
        xticks=[calc.deg_to_rad(0),calc.deg_to_rad(90),calc.deg_to_rad(180)]
        xticklabels=[0,90,180]
        yticks=[0,0.075]
        yticklabels=[0,0.075]
        ax.get_xaxis().set_ticks(xticks)
        ax.get_xaxis().set_ticklabels(xticklabels,fontsize=6)
        ax.get_yaxis().set_ticks(yticks)
        ax.get_yaxis().set_ticklabels(yticklabels,fontsize=6)
        ax.set_xlabel('heading difference\n($^\circ$)', fontsize=6)
        ax.set_ylabel('probability', fontsize=6)
        ax.xaxis.labelpad = 0
        ax.yaxis.labelpad=-12


    def plot_vec_strength_and_dir(self,axvec,axhist,axvec_nonfrontal):

        hst_bnds=[0,1]
        num_bins=10
        COLCT=0
        concat_dt_vec={}
        dt_vec_in_cols={}
        concat_dt_vec_front={}
        concat_dt_head={}
        self.concat_dt_head_in_cols={}
        dir_dt={}
        dir_dt_all={}
        all_dt_dict={}
        self.comb_head_dt={}
        for fly_type_ind,flytype in enumerate(FLY_TYPE_TO_PLOT):
            horiz_pos=0.2
            vert_pos=0.2

            for stim_type in ['sun']:
                frontal_dt={}
                for indnum in [0,1]:
                    #plot vec strength
                    indt=np.array(self.adt[flytype][indnum]['vec_strength'])
                    if indnum==1:
                        if flytype is 'tb_flies':
                            PLOTFLAG=False
                        else:
                            PLOTFLAG=True
                    else:
                        PLOTFLAG=True
                    if PLOTFLAG:
                       vert_offset=fly_type_ind
                       crcol=COLVLS[COLCT]
                       ax=axvec[indnum]
                       self.make_hist_cumulative(ax,indt,horiz_pos,vert_pos,vert_offset,crcol,hst_bnds)

                        #axvec[indnum].text(0.2,0.2+fly_type_ind*.2,str(len(indt)),color=COLVLS[COLCT])
                        #axvec[indnum].plot(np.mean(indt),0.1,'v',color=COLVLS[COLCT])
                        #cr_str="%.2f"%np.mean(indt)
                        #axvec[indnum].text(np.mean(indt),0.02*fly_type_ind*.02,cr_str,color=COLVLS[COLCT],fontsize=6)
                        #twplt.plot_hist(axvec[indnum],indt,hst_bnds=hst_bnds,num_bins=100,plt_type='cumulative',col=COLVLS[COLCT])



                    #plot position
                    indt=np.array(self.adt[flytype][indnum]['mnrad'])
                    thresh_inds=np.where(np.array(self.adt[flytype][indnum]['vec_strength'])>VEC_STRENGTH_THRESH)[0]
                    low_vec_inds=np.where(np.array(self.adt[flytype][indnum]['vec_strength'])<VEC_STRENGTH_THRESH)[0]
                    #map to 180 degrees
                    inds=np.where(indt[thresh_inds]>np.pi)[0]
                    dir_dt[indnum]=indt[thresh_inds]
                    dir_dt[indnum][inds]=2*np.pi-indt[thresh_inds][inds]
                    dir_dt_all[indnum]=np.array(self.adt[flytype][indnum]['mnrad'])
                    alldt=np.array(self.adt[flytype][indnum]['mnrad'])
                    
                    dt_to_add=np.copy(alldt)
                    dt_to_add[low_vec_inds]=np.nan
                    
                    all_dt_dict[indnum]=dt_to_add

                    alldt_inds=np.where(alldt>np.pi)[0]
                    
                    dir_dt_all[indnum][alldt_inds]=2*np.pi-alldt[alldt_inds]

                    hst_bnds=[0,np.pi]
                    ax=axhist[indnum]
                    self.make_hist_cumulative(ax,dir_dt[indnum],horiz_pos,vert_pos,vert_offset,crcol,hst_bnds)

                    #axhist[indnum].text(0.2,0.2+fly_type_ind*.2,str(np.count_nonzero(~np.isnan(thresh_inds))),color=COLVLS[COLCT])
                    #cr_str="%.2f"%np.mean(dir_dt[indnum])
                    #axhist[indnum].text(np.mean(dir_dt[indnum]),0.02*fly_type_ind*.02,cr_str,color=COLVLS[COLCT],fontsize=6)
                    #axhist[indnum].plot(np.mean(dir_dt[indnum]),0.1,'v',color=COLVLS[COLCT])
                    #twplt.plot_hist(axhist[indnum],dir_dt[indnum],hst_bnds=[0,np.pi],plt_type='cumulative',num_bins=100,col=COLVLS[COLCT])


                    #find position where angle is not frontal according to threshold.
                    ax=axvec_nonfrontal[indnum]
                    nonfrontal_inds_high=np.where(indt>FRONTAL_THRESH)[0]
                    nonfrontal_inds_low=np.where(indt<(2*np.pi-FRONTAL_THRESH))[0]
                    nonfrontal_inds=np.intersect1d(nonfrontal_inds_low,nonfrontal_inds_high)
                    indt=np.array(self.adt[flytype][indnum]['vec_strength'])[nonfrontal_inds]
                    frontal_dt[indnum]=np.reshape(indt,len(indt))
                    hst_bnds=[0,1]
                    self.make_hist_cumulative(ax,indt,horiz_pos,vert_pos,vert_offset,crcol,hst_bnds)
                #combined 1st and 2nd trials
                indnum=2
                #vec strength
                crdt=np.concatenate((np.array(self.adt[flytype][0]['vec_strength']),np.array(self.adt[flytype][1]['vec_strength'])),axis=0)
                crdt_in_cols=np.concatenate((np.array(self.adt[flytype][0]['vec_strength']),np.array(self.adt[flytype][1]['vec_strength'])),axis=1)
                crax=axvec[indnum]
                self.make_hist_cumulative(crax,crdt,horiz_pos,vert_pos,vert_offset,crcol,hst_bnds)
                concat_dt_vec[flytype]=np.reshape(crdt,len(crdt))
                dt_vec_in_cols[flytype]=crdt_in_cols
                crax=axvec_nonfrontal[indnum]

                front_dt_combined=np.concatenate((frontal_dt[0],frontal_dt[1]))
                concat_dt_vec_front[flytype]=front_dt_combined
                self.make_hist_cumulative(crax,front_dt_combined,horiz_pos,vert_pos,vert_offset,crcol,hst_bnds)

                #axvec[indnum].text(0.2,fly_type_ind*.1,str(np.count_nonzero(~np.isnan(crdt))),color=COLVLS[COLCT])

                #axvec[indnum].plot(np.mean(crdt),0.1,'v',color=COLVLS[COLCT])
                #twplt.plot_hist(axvec[indnum],crdt,hst_bnds=hst_bnds,num_bins=100,plt_type='cumulative',col=COLVLS[COLCT])
                #concat_dt_all[flytype]=np.reshape(indt,len(indt))

                #position

                crdt=np.concatenate((dir_dt[0],dir_dt[1]))

                self.comb_head_dt[flytype]=np.concatenate((all_dt_dict[0],all_dt_dict[1]))
               
                cr_hd_dt_in_cols=np.concatenate((dir_dt_all[0],dir_dt_all[1]),axis=1)
               
                concat_dt_head[flytype]=crdt
                for colind in [0,1]:
                    low_vec_strength_inds=np.where(crdt_in_cols[:,colind]<VEC_STRENGTH_THRESH)[0]
                    cr_hd_dt_in_cols[low_vec_strength_inds,colind]=np.nan
                
                self.concat_dt_head_in_cols[flytype]=cr_hd_dt_in_cols


                axhist[indnum].text(0.2,fly_type_ind*.1,str(len(crdt)),color=COLVLS[COLCT])
                
                axhist[indnum].plot(np.mean(crdt),0.1,'v',color=COLVLS[COLCT])
                twplt.plot_hist(axhist[indnum],crdt,hst_bnds=[0,np.pi],plt_type='cumulative',num_bins=100,col=COLVLS[COLCT])

            for stim_type in ['stripe']:


                indt=np.array(self.adt[flytype][2]['vec_strength'])
                try:
                    axvec[3].text(0.2,fly_type_ind*.1,str(np.count_nonzero(~np.isnan(indt))),color=COLVLS[COLCT])

                    twplt.plot_hist(axvec[3],indt,hst_bnds=hst_bnds,num_bins=100,plt_type='cumulative',col=COLVLS[COLCT])

                    indt=np.array(self.adt[flytype][2]['mnrad'])
                    thresh_inds=np.where(np.array(self.adt[flytype][2]['vec_strength'])>VEC_STRENGTH_THRESH)[0]
                        #map to 180 degrees

                    inds=np.where(indt[thresh_inds]>np.pi)[0]
                    pltdt=indt[thresh_inds]
                    pltdt=2*np.pi-indt[thresh_inds][inds]
                    axhist[3].text(0.2,0.2+fly_type_ind*.2,str(len(thresh_inds)),color=COLVLS[COLCT])

                    twplt.plot_hist(axhist[3],pltdt,hst_bnds=[0,np.pi],plt_type='cumulative',num_bins=100,col=COLVLS[COLCT])
                except:
                    print('stripe error')
            COLCT=COLCT+1

        FLY_TYPE_TO_PLOT
        
        try:
            if PERMUTE_BY_INDIVIDUALS_NOT_FLIGHTS:
                
                #sd,pctile_all=calc.permute_test_for_mean_diff_sampling_from_rows(dt_vec_in_cols[FLY_TYPE_TO_PLOT[0]], dt_vec_in_cols[FLY_TYPE_TO_PLOT[2]], num_permutations=1000)
                sd,pctile_head=calc.permute_test_for_mean_diff_sampling_from_rows(self.concat_dt_head_in_cols[FLY_TYPE_TO_PLOT[0]], self.concat_dt_head_in_cols[FLY_TYPE_TO_PLOT[1]], num_permutations=1000)
            else:
                 sd,pctile_front=calc.permute_test_for_mean_diff_between_two_groups(concat_dt_vec_front[FLY_TYPE_TO_PLOT[0]], concat_dt_vec_front[FLY_TYPE_TO_PLOT[1]], num_permutations=1000)
                 sd,pctile_head=calc.permute_test_for_mean_diff_between_two_groups(concat_dt_head[FLY_TYPE_TO_PLOT[0]], concat_dt_head[FLY_TYPE_TO_PLOT[1]], num_permutations=1000)
            
            
        except:
            print ('permutation error')
        for ind in [0,1,2,3]:
            fpl.adjust_spines(axvec[ind],['left','bottom'])
            fpl.adjust_spines(axhist[ind],['left','bottom'])
            fpl.adjust_spines(axvec_nonfrontal[ind],['left','bottom'])
            #crax.set_ylim([-calc.deg_to_rad(20.),calc.deg_to_rad(380.0)])
                #crax.plot(0.22,mnvl_in_rad,'r<')
            for crax in [axvec[ind],axvec_nonfrontal[ind]]:
                crax.set_xlim([0,1.0])
                crax.set_ylim([0,1])
                crax.get_xaxis().set_ticks([0,0.25,0.5,0.75,1.0])
                crax.get_xaxis().set_ticklabels(['0','0.25','0.5','0.75'],fontsize=6)
                crax.get_yaxis().set_ticks([0,0.5,1])
                crax.get_yaxis().set_ticklabels(['0','0.5','1'],fontsize=6)
                crax.set_xlabel('vector strength',fontsize=6)
                crax.set_ylabel('probability',fontsize=6)
                crax.set_aspect(1)

            axhist[ind].set_xlim([0,np.pi])
            axhist[ind].set_ylim([0,1])
            axhist[ind].get_xaxis().set_ticks([0,np.pi/2,np.pi])
            axhist[ind].get_xaxis().set_ticklabels(['0','90','180'],fontsize=6)
            axhist[ind].get_yaxis().set_ticks([0,0.5,1])
            axhist[ind].get_yaxis().set_ticklabels(['0','0.5','1'],fontsize=6)
            axhist[ind].set_xlabel('mean heading (deg)',fontsize=6)
            axhist[ind].set_ylabel('probability',fontsize=6)
            #axvec[ind].set_aspect(1)


        # for ind in [3]:
        #     fpl.adjust_spines(ax[ind],['left','bottom'])

        #     #crax.set_ylim([-calc.deg_to_rad(20.),calc.deg_to_rad(380.0)])
        #         #crax.plot(0.22,mnvl_in_rad,'r<')
        #     axvec[ind].set_xlim([0,1])
        #     axvec[ind].set_ylim([0,1])
        #     axvec[ind].get_xaxis().set_ticks([0,0.25,0.5,0.75,1])
        #     axvec[ind].get_xaxis().set_ticklabels(['0','0.25','0.5','0.75','1'],fontsize=6)
        #     axvec[ind].get_yaxis().set_ticks([0,0.5,1])
        #     axvec[ind].get_yaxis().set_ticklabels(['0','0.5','1'],fontsize=6)
        #     axvec[ind].set_xlabel('vector strength',fontsize=6)
        #     axvec[ind].set_ylabel('probability',fontsize=6)
        #     axvec[ind].set_aspect(1)


    def make_hist_cumulative(self,ax,indt,horiz_pos,vert_pos,vert_offset,crcol,hst_bnds):
        nvl=len(indt)
        mnvl=np.mean(indt)
        ax.text(horiz_pos,vert_pos+vert_offset*.2,nvl,color=crcol)
        
        ax.plot(mnvl,0.1,'v',color=crcol)
        cr_str="%.2f"%mnvl
        ax.text(np.mean(indt),vert_pos+vert_offset*.02,cr_str,color=crcol,fontsize=6)
        twplt.plot_hist(ax,indt,hst_bnds=hst_bnds,num_bins=100,plt_type='cumulative',col=crcol)

    def plot_hist_dir(self,ax):
        COLCT=0
        for fly_type_ind,flytype in enumerate(FLY_TYPE_TO_PLOT):
            for indnum in [0,1,2]:

                indt=np.array(self.adt[flytype][indnum]['mnrad'])

                thresh_inds=np.where(np.array(self.adt[flytype][indnum]['vec_strength'])>VEC_STRENGTH_THRESH)[0]
                crdt=np.cos(indt[thresh_inds])

                ax[indnum].text(0.2,0.2+fly_type_ind*.2,str(np.count_nonzero(~np.isnan(thresh_inds))),color=COLVLS[COLCT])

                twplt.plot_hist(ax[indnum],crdt,hst_bnds=[-1,1],plt_type='cumulative',num_bins=100,col=COLVLS[COLCT])
            COLCT=COLCT+1


        xvl=np.linspace(0,2*np.pi,10000)
        yvl=np.cos(xvl)
        for indnum in [0,1,2]:
            twplt.plot_hist(ax[indnum],yvl,hst_bnds=[-1,1],num_bins=100,plt_type='cumulative',col=[0.5,0.5,0.5])

        for indnum in [0,1,2]:
            crax=ax[indnum]
            fpl.adjust_spines(crax,['left','bottom'])

            #crax.set_ylim([-calc.deg_to_rad(20.),calc.deg_to_rad(380.0)])
                #crax.plot(0.22,mnvl_in_rad,'r<')
            #crax.set_xlim([0,0.2])
            crax.set_xlim([-1,1])
            crax.set_ylim([0,1])
            crax.get_xaxis().set_ticks([-1,0,1])
            crax.get_xaxis().set_ticklabels(['-1','0','1'],fontsize=6)
            crax.get_yaxis().set_ticks([0,0.5,1])
            crax.get_yaxis().set_ticklabels(['0','0.5',1],fontsize=6)
            crax.set_xlabel('cosine of heading',fontsize=6)
            crax.set_ylabel('cumulative probability',fontsize=6)
            crax.set_aspect(2)






    def plot_heading_hist(self,ax):
        COLCT=0
        for fly_type_ind,flytype in enumerate(FLY_TYPE_TO_PLOT):
            for stim_type in ['sun','stripe']:
                indt=np.array(self.adt[flytype][stim_type]['mnrad'])

                thresh_inds=np.where(np.array(self.adt[flytype][stim_type]['vec_strength'])>VEC_STRENGTH_THRESH)[0]
                crdt=indt[thresh_inds]

                highinds=np.where(crdt>np.pi)
                crdt[highinds]=-(2*np.pi-crdt[highinds])

                twplt.plot_hist(ax,crdt,hst_bnds=[-np.pi-np.pi/2,np.pi+np.pi/2],num_bins=12,col=COLVLS[COLCT])

                COLCT=COLCT+1

        fpl.adjust_spines(ax,['left','bottom'])

        #crax.set_ylim([-calc.deg_to_rad(20.),calc.deg_to_rad(380.0)])
            #crax.plot(0.22,mnvl_in_rad,'r<')
        #crax.set_xlim([0,0.2])
        ax.set_xlim([-np.pi,np.pi])
        ax.get_xaxis().set_ticks([-np.pi,0,np.pi])
        ax.get_xaxis().set_ticklabels(['-180','0','180'],fontsize=6)
        ax.get_yaxis().set_ticks([0,0.7])
        ax.get_yaxis().set_ticklabels(['0','0.7'],fontsize=6)
        ax.set_xlabel('mean heading',fontsize=6)
        ax.set_ylabel('probability',fontsize=6)
        ax.set_aspect(7)


    def make_plot_of_scatter_headings(self,ax,crtype):
        thresh=0.2
        vecvls=np.array(self.adt[crtype]['sun']['vec_strength'])
        inds1=np.where(vecvls[:,0]>thresh)
        inds2=np.where(vecvls[:,1]>thresh)

        comb_inds=np.intersect1d(inds1[0],inds2[0])

        pltvls=np.array(self.adt[crtype]['sun']['mnrad'])[comb_inds,:]








        twplt.scatterplot(ax,pltvls[:,0],pltvls[:,1],sizefactor=0.3,double_horizontal_ax=False)
        
        twplt.scatterplot(ax,pltvls[:,0],pltvls[:,1]+360,sizefactor=0.3,double_horizontal_ax=False)



    def plot_all_raw(self,crtype):

        #self.make_axes_for_raw(**kwargs)

        self.plot_all_raw_dt(crtype)

    def make_axes_for_raw(self,kwargs):
        #first determine how many distinct svg files we will need to open
        num_svg_files
        self.raw_layout=[]
        self.totalaxct=0
        for cr_svg_file_num in np.arange(num_svg_files):
            cr_svg=raw_plot_svg
            self.raw_layout.append(fifi.FigureLayout(svg_fig))
            self.raw_layout[-1].make_mplfigures()
            self.add_cr_layout_to_master_axes()
            self.totalaxct=self.totalaxct+num_row_per_raw_svg


    def add_cr_layout_to_master_axes():


        #gs = gridspec.GridSpec(NROW,NCOL)
        keylists=self.make_keylists()
        for crkeynum, cr_time_key in enumerate(timedt_keylist):
            ax['timedt'][crkey+self.totalaxct]=self.raw_layout[-1].axes[cr_timekey]
            cr_stripe_key=stripe_keylist[crkeynum]
            cr_sun_key_a=suna_keylist[crkeynum]
            cr_sun_key_b=sunb_keylist[crkeynum]
            ax['stripehst'][crkey+self.totalaxct]=self.raw_layout[-1].axes[cr_timekey]
            ax['sunhsta'][crkey+self.totalaxct]=self.raw_layout[-1].axes[cr_sun_key_a]
            ax['sunhstb'][crkey+self.totalaxct]=self.raw_layout[-1].axes[cr_sun_key_b]

        return ax




    def plot_all_raw_dt(self,exp_type_to_plot):
        cols_for_array=3
        gridpec_cols=21

        flies_to_plot=self.sumdt[exp_type_to_plot]
        self.fig={}
        gs={}
        axmst={}
        axmotor={}
        axhist={}
        for key in ['stripe','sun']:
            axmotor[key]=[]
            axhist[key]=[]

        axtext=[]
        axwing=[]
        SEC_FIG_INITIALIZED=False

        figct=0



        self.crplotnum=0


        #create order to plot
        vec_dict={}

        for ind in [0,1]:
            vec_strength_list=self.adt[exp_type_to_plot][ind]['vec_strength']
            try:
                vec_dict[ind]=1-np.reshape(vec_strength_list,len(vec_strength_list))
            except:
                print('reshape error')

        try:
            mnvls=np.nanmax([vec_dict[0],vec_dict[1]],axis=0)
        except:
            print('nanmax error')
        indlist=np.argsort(mnvls)

        for flyindnum,crfly in enumerate(np.array(self.adt[exp_type_to_plot][0]['flynum'])[indlist]):

            self.get_data(crfly)


            self.make_raw_plot(flyindnum,exp_type_to_plot)

        self.ordered_indlist=indlist
        #for crfig in np.arange(num_figs):
         #   self.fig[crfig].savefig(self.datapath+'example'+str(crfig)+'.pdf')

    def make_raw_plot(self,flyindnum,exp_type_to_plot):
        ADD_VEC_TEXT=True
        #COLNUM=-1

        if self.crdt:


            for cr_fltnum in self.crdt.keys():

                if self.crdt[cr_fltnum]:




                    mnvl_in_rad=self.crdt[cr_fltnum]['mnrad_360']
                    if mnvl_in_rad>np.pi:
                        mnvl_in_rad=-(2*np.pi-mnvl_in_rad)
                    halt_flag=False

                    offset_time=0
                    if cr_fltnum==1:
                        offset_time=self.crdt[cr_fltnum-1]['time_in_min'][-1]
                    elif cr_fltnum>1:
                        offset_time=self.crdt[cr_fltnum-1]['time_in_min'][-1]-TIME_GAP
                    if flyindnum<len(self.axraw[exp_type_to_plot]['tmdt']):
                        try:
                            fpb.plot_motor(self.crdt[cr_fltnum],self.axraw[exp_type_to_plot]['tmdt'][self.crplotnum],plot_vector=False,plot_split=1,plot_start_angle=0,subtract_zero_time=True,offset_time=offset_time,plot_vert_line_at_end=True, halt_flag=halt_flag,center_on_zero_flag=True,withhold_bottom_axis=True)
                        except:
                            print('plot_motor error')
                        self.axraw[exp_type_to_plot]['tmdt'][self.crplotnum].set_xlim([0,15.5])


                        if ADD_VEC_TEXT:
                            crvec=1-self.crdt[cr_fltnum]['circvar_mod360']

                            self.axraw[exp_type_to_plot]['tmdt'][self.crplotnum].text(1.5*5*cr_fltnum,30,str(crvec))
                        #if COLNUM:

                         #   axmotor[crkey][flyindnum][COLNUM].axis('off')
                          #  axhist[crkey][flyindnum][COLNUM].axis('off')
                        mot_deg=calc.center_deg_on_zero(calc.rad_to_deg(self.crdt[cr_fltnum]['mot_rad']))
                        degbins=np.arange(-180,190,10)
                        tst=calc.make_hist_calculations(mot_deg,degbins)
                        deg_per_bin=10

                        tst['rad_per_bin']=deg_per_bin*np.pi/180

                        tst['xrad']=calc.deg_to_rad(degbins)
                        try:
                            crax=self.axraw[exp_type_to_plot]['hst'][self.crplotnum][cr_fltnum]
                        except:
                            print('ax assignment error')

                        crax.step(tst['normhst'],tst['xrad'][0:-1]+tst['rad_per_bin']/2,color='k',linewidth=0.5)
                        if PLOT_EX_TEXT_STR:
                            textax=self.axraw[exp_type_to_plot]['txt'][self.crplotnum]
                            txtfile=self.crdt[cr_fltnum]['fname'].split('/')[-1]
                            textax.text(0,0,txtfile,fontsize=4, rotation='vertical')
                            textax.set_ylim([-12,2])

                        #crax.step(self.crdt[cr_fltnum]['normhst'],self.crdt[cr_fltnum]['xrad'][0:-1]+self.crdt[cr_fltnum]['rad_per_bin']/2,'k',linewidth=0.5)
                        #self.col_num[crkey]=self.col_num[crkey]+1
                        fpl.adjust_spines(crax,[])
                        crax.set_ylim([-calc.deg_to_rad(180.),calc.deg_to_rad(180.0)])
                        crax.plot(0.21,mnvl_in_rad,'r<',clip_on=False)
                        crax.set_xlim([0,0.24])
            self.crplotnum=self.crplotnum+1




    def add_text_name(self,ax):
        fname=self.crdt['sun'][0]['fname']
        cutname=fname.split('/cloop')
        #tst=['', 'home', 'timothy', 'data', '20170505', 'fly1-rnach']
        tst=cutname[0].split('/')
        #printstr='20170505/fly1-rnach'
        printstr=tst[-2]+'/'+tst[-1]
        ax.text(0,0,printstr,fontsize=5)
        ax.axis('off')

    def plot_net_displacement(self,ax, **kwargs):

        offset_angle=np.pi/2
        mn_rad_lst=[]
        vec_strength_lst=[]
        sumdt={}
        COLCT=0
        for fly_type_ind,flytype in enumerate(FLY_TYPE_TO_PLOT):

            if 'inact' in flytype:
                plot_key='tnt_inactive'
                text_key='intrinsic \n tnt'

            elif flytype is 'gal_tnt_flies':
                plot_key='tnt_active'
                text_key='intrinsic \n tnt-inactive'
            elif flytype is 'wedge096_kir_flies':
                plot_key='wedge_kir_flies'
                text_key='wedge Kir'
            elif flytype is 'wedge_kir_flies':
                plot_key='wedge_kir_flies'
                text_key='wedge Kir'
            elif flytype is 'tb_flies':
                plot_key='tb_flies'
                text_key='TB'
            elif flytype is 'split_vector_flies':
                plot_key=flytype
                text_key='empty_Kir'
            plot_left_axis=True


            for ind in [0,1,2]:

                crax=ax[flytype][ind]
                cr_rad=self.adt[flytype][ind]['mnrad']
                cr_vec=self.adt[flytype][ind]['vec_strength']
                if COMBINE_FIRST_AND_SECOND_DISPLACEMENT_FLAG:
                    if ind==0:
                        cr_vec_comb=np.concatenate((np.array(self.adt[flytype][0]['vec_strength']),np.array(self.adt[flytype][1]['vec_strength'])),axis=0)
                        cr_rad_comb=np.concatenate((np.array(self.adt[flytype][0]['mnrad']),np.array(self.adt[flytype][1]['mnrad'])),axis=0)
                else:
                    cr_vec=np.array(self.adt[flytype][ind]['vec_strength'])
                    cr_rad=np.array(self.adt[flytype][ind]['mnrad'])

                try:
                    mn_rad=np.reshape(cr_rad,len(cr_rad))+offset_angle
                except:
                    print('reshape error')
                try:
                    mn_rad_no_offset=np.reshape(cr_rad,len(cr_rad))
                except:
                    print('reshape error')
                try:
                    vec_strength=np.reshape(cr_vec,len(cr_vec))
                except:
                    print('reshape error')
                plt_rad = mn_rad[~np.isnan(mn_rad)]

                
                plt_vec=vec_strength[~np.isnan(vec_strength)]
                if PLOT_MEAN_DISPLACEMENT:
                    od=self.summary_stats(mn_rad_no_offset[~np.isnan(mn_rad_no_offset)],plt_vec)
                #twplt.scatter_polar(crax,plt_rad,plt_vec,withhold_direction_ticks=WITHHOLD_LEGEND_POLAR,plot_mean=PLOT_MEAN_DISPLACEMENT,sumstats=od)
                if COMBINE_FIRST_AND_SECOND_DISPLACEMENT_FLAG:
                    
                        
                    
                    if ind==0:
                        comb_pctiles=calc.resample_for_confidence_intervals(self.comb_head_dt[flytype],rows_flag=True)
                        twplt.scatter_polar(ax[flytype][3],np.reshape(cr_rad_comb,len(cr_rad_comb))+offset_angle,np.reshape(cr_vec_comb,len(cr_vec_comb)),withhold_direction_ticks=WITHHOLD_LEGEND_POLAR,plot_mean=PLOT_MEAN_DISPLACEMENT,sumstats=od)
                        ax[flytype][3].plot(np.linspace(comb_pctiles[0]+np.pi/2,comb_pctiles[2]+np.pi/2,100),.98*np.ones(100),'c',clip_on=False)
                        ax[flytype][3].plot([comb_pctiles[1]+np.pi/2,comb_pctiles[1]+np.pi/2],[0,1],'c')
                    #if ind==1:
                        #crind=self.ordered_fly_indices[flytype][1]
                        
                        #crax.plot(plt_rad[crind],0.9,'r<','MarkerSize',1)
                    crax.plot([0,np.pi],[VEC_STRENGTH_THRESH, VEC_STRENGTH_THRESH],'c')
                if ind==2:
                    low_inds=np.where(np.array(cr_vec)<VEC_STRENGTH_THRESH)
                    rad_copy=np.copy(np.array(cr_rad))
                    rad_copy[low_inds]=np.nan
                    

                    stripe_pctiles=calc.resample_for_confidence_intervals(rad_copy,rows_flag=False)
                    
                    crax.plot(np.linspace(stripe_pctiles[0]+np.pi/2,stripe_pctiles[2]+np.pi/2,100),.98*np.ones(100),'c',clip_on=False)
                
                    crax.plot([stripe_pctiles[1]+np.pi/2,stripe_pctiles[1]+np.pi/2],[0,1],'c')
                twplt.scatter_polar(crax,plt_rad,plt_vec,withhold_direction_ticks=WITHHOLD_LEGEND_POLAR,plot_mean=PLOT_MEAN_DISPLACEMENT,sumstats=od)
                
               
                crax.text(3*np.pi/2,1.4,flytype,fontsize=8)
            

                #    if plot_left_axis:
                 #       crax.text(np.pi,1.7,text_key+'\nn='+str(len(plt_vec)),color=COLVLS[COLCT],fontsize=7)

            #COLCT=COLCT+1


    def set_up_plot_ax(self,crflytype):
        LAYOUTIND=0

        ax={}
        ax['tmdt']=[]
        ax['hst']=[]
        ax['txt']=[]
        for cr_row in np.arange(NROW):
            ax['hst'].append([])

            axrow=np.mod(cr_row,MAXROW)
            if cr_row>MAXROW-1:
                LAYOUTIND=int(np.floor(cr_row/MAXROW))
            if PLOT_EX_TEXT_STR:
                hst_3_str='txt%s'%(str(axrow))
                try:
                    ax['txt'].append(self.layout[crflytype][LAYOUTIND].axes[hst_3_str])
                except:
                    print('text append error')
                fpl.adjust_spines(ax['txt'][-1],['none'])

                #hst_4_str='hst%s4'%(str(axrow))
                #axtmp=self.layout[crflytype][LAYOUTIND].axes[hst_4_str]
                #fpl.adjust_spines(axtmp,['none'])
            for crcol in np.arange(NHISTCOL):

                crhststr='hst%s%s'%(str(axrow),str(crcol))
                try:
                    ax['hst'][cr_row].append(self.layout[crflytype][LAYOUTIND].axes[crhststr])
                except:
                    print('layout error')


            crtmstr='ex%s'%str(axrow)
            try:
                ax['tmdt'].append(self.layout[crflytype][LAYOUTIND].axes[crtmstr])
            except:

                print('layout error')

        return ax

    def set_up_sum_plot_ax(self,crlayout):
        #NUMDT=5
        ax={}

        for crind,crfield in enumerate(['sum0','sum1','sum2','sum3','sum4']):
            try:
                ax[crind]=crlayout.axes[crfield]
            except:
                pdb.set_trace()
        for corr_field in ['corr_scatter','corr_hist']:
            ax[corr_field]=crlayout.axes[corr_field]
        return ax



    def plot_displacement(self,ax, **kwargs):

        try:
            plot_net_displacement_flag=kwargs['plot_net_displacement']
        except:
            plot_net_displacement=False
        plot_positions={}
        for key in ['x','y']:
            plot_positions[key]=[]
        if not plot_net_displacement:
            for ind in np.arange(len(self.crdt)):
                try:
                    tst=self.crdt[ind]['displacement_traj']['raw']
                except:
                    pdb.set_trace()
                for key in ['x','y']:
                    if ind:

                        plot_positions[key].append(plot_positions[key][0][-1]+tst[key])
                    else:
                        plot_positions[key].append(tst[key])






        polar_positions=calc.linear_to_polar(plot_positions)
        col=np.random.permutation(['r','b','g','m','c','k'])[0]

        indvls=np.arange(len(polar_positions['theta']))

        if kwargs['init_flag']:
            col='b'
            axsp=kwargs['init_ax']
            for ind in indvls:
                axsp.plot(polar_positions['theta'][ind]+np.pi/2,polar_positions['len'][ind],color=col)
                axsp.plot(polar_positions['theta'][ind][-1]+np.pi/2,polar_positions['len'][ind][-1],'o',color='k')
            axsp.get_xaxis().set_ticks([0,np.pi/2.,np.pi,3.*(np.pi/2.)])
            axsp.get_xaxis().set_ticklabels(['0','90','180','270'],fontsize=8)
            axsp.get_yaxis().set_ticks([])
            axsp.get_yaxis().set_ticklabels([],fontsize=8)

        for ind in indvls:
            ax.plot(polar_positions['theta'][ind]+np.pi/2,polar_positions['len'][ind],color=col)
            ax.plot(polar_positions['theta'][ind][-1]+np.pi/2,polar_positions['len'][ind][-1],'o',color='k')

        ax.get_xaxis().set_ticks([0,np.pi/2.,np.pi,3.*(np.pi/2.)])
        ax.get_xaxis().set_ticklabels(['0','90','180','270'],fontsize=8)
        ax.get_yaxis().set_ticks([])
        ax.get_yaxis().set_ticklabels([],fontsize=8)

    def summary_stats(self,rad_lst,vec_strength_lst):
        out_dict={}
        try:

            out_dict['wtmn'],out_dict['wtvar']=calc.weighted_mean(rad_lst,vec_strength_lst,anal_180=False)
        except:
            pdb.set_trace()
        out_dict['mn'],out_dict['var']=calc.circ_mean(rad_lst,anal_180=False)
        
        inds=[]
        inds.append(len(np.intersect1d(np.where(rad_lst>np.pi/4)[0],np.where(rad_lst<(np.pi/4+np.pi/2))[0])))
        inds.append(len(np.intersect1d(np.where(rad_lst>0.75*np.pi)[0],np.where(rad_lst<(1.25*np.pi))[0])))
        inds.append(len(np.intersect1d(np.where(rad_lst>1.25*np.pi)[0],np.where(rad_lst<calc.deg_to_rad(315))[0])))
        initinds=[len(rad_lst)-sum(inds)]

        out_dict['dist']=initinds+inds

        return out_dict



    def get_data(self,crfly):
        crexp=0
        crind=0

        try:
            split_name=self.sumdt[crfly][crexp]['sumstats'][crind]['fname'].split('Dropbox/')
            DAT_FLAG=True
        except:
            DAT_FLAG=False

        if DAT_FLAG:
            if len(split_name)<2:
                split_name=self.sumdt[crfly][crexp]['sumstats']['sun'][0]['fname'].split('data/')

            datpath=dropbox_path+split_name[1].strip('.txt')+'_rawdata.pck'
            params_datapath=dropbox_path+split_name[1].strip('.txt')+'.pck'
            try:
                self.crdt=fh.open_pickle(datpath)

                self.params=fh.open_pickle(params_datapath)
            except:
                pdb.set_trace()
        else:
            self.crdt=[]



if __name__ == '__main__':
    import sys
    example_plot=Example_Plot(input_args=sys.argv)
    example_plot.make_plots()
