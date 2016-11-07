import os
import sys
import numpy as np
import random as rand
import math as m

Ra_aim = 
Dec_aim = 

optical_cat = ""
radio_cat = ""
xray_cat = ""

def where(array,op,ref):
	temp_sort_ind = np.argsort(array)
	temp_sorted = array[temp_sort_ind]
	if op == "eq":
		up = np.searchsorted(array,ref,side='r')
		down = np.searchsorted(array,ref,side='l')
		if up == down: 
			out_inds = np.zeros(0)
		else:
			out_inds = temp_sort_ind[down:up]
	if op == "ge":
		down = np.searchsorted(array,ref,side='l')
		if down > len(array): 
			out_inds = np.zeros(0)
		else:
			out_inds = temp_sort_ind[down:len(array)+1]
	if op == "gt":
		down = np.searchsorted(array,ref,side='r')
		if down > len(array): 
			out_inds = np.zeros(0)
		else:
			out_inds = temp_sort_ind[down:len(array)+1]
	if op == "le":
		down = np.searchsorted(array,ref,side='r')
		if down <= 0: 
			out_inds = np.zeros(0)
		else:
			out_inds = temp_sort_ind[0:down]
	if op == "lt":
		down = np.searchsorted(array,ref,side='l')
		if down <= 0: 
			out_inds = np.zeros(0)
		else:
			out_inds = temp_sort_ind[0:down]
	return out_inds



def SphDist(RAi,Deci,Raf,Decf):
	dist = 2.0*m.asin(m.sqrt((m.sin(0.5*m.pi*(Decf-Deci)/360))**2+m.cos(Decf)*m.cos(Deci)*(m.sin(0.5*m.pi*(RAf-RAi)/360))**2))
	return dist

def FindCloseSources(ra,dec,tol,ra_opt,dec_opt):
	ra_box_temp_ind1 = where(ra_opt,"ge",ra-tol)
	ra_box_temp1 = ra_opt[ra_box_temp_ind1]
	ra_box_temp_ind2 = where(ra_box_temp1,"le",ra+tol)
	ra_box_temp2 = ra_box_temp1[ra_box_temp_ind2]
	if len(ra_box_temp2) > 0:
		dec_box_temp1 = dec_opt[ra_box_temp_ind1[ra_box_temp_ind2]]
		dec_box_temp_ind2 = where(dec_box_temp1,"ge",dec-tol)
		dec_box_temp2 = dec_box_temp1[dec_box_temp_ind2]
		dec_box_temp_ind3 = where(dec_box_temp2,"le",dec+tol)
		dec_box_temp3 = dec_box_temp2[dec_box_temp_ind3]
		if len(dec_box_temp3) > 0:
			ra_box_temp3 = ra_box_temp2[dec_box_temp_ind2[dec_box_temp_ind3]]
			temp_dist_ar = np.zeros(len(dec_box_temp3))
			for i in range(0,len(temp_dist_ar)):
				temp_dist_ar[i] = SphDist(ra_box_temp3[i],dec_box_temp3[i],ra,dec)
			inside_tol_ind = where(temp_dist_ar,"le",tol)
			return ra_box_temp_ind1[ra_box_temp_ind2[dec_box_temp_ind2[dec_box_temp_ind3[inside_tol_ind]]]]
		else:
			return np.zeros(0)
	else:
		return np.zeros(0)

def OptMatch(raX,decX,d_xerr,ra_opt,dec_opt):
	close_ind = FindCloseSources(raX,decX,d_xerr)
	ra_close = ra_opt[close_ind]
	dec_close = ra_opt[close_ind]
	magI_close = magI_opt[close_ind]
	numdens_close = loc_num_dens_gtMag[close_ind]
	close_likelihoods = np.zeros(len(close_ind))
	close_dist = np.zeros(len(close_ind))
	for j in range(0,len(close_ind)):
		close_dist[j] = SphDist(raX,decX,ra_close[j],dec_close[j])
		close_likelihoods[j] = m.exp(-0.5*(close_dist[j]**2)/(d_xerr))/(d_xerr**2*numdens_close[j])
	return close_likelihoods,close_ind
	
def MonteCarlo(reps,ra_UB,ra_LB,dec_UB,dec_LB,err,ra_opt,dec_opt):
	rand.seed()
	output = np.zeros(0)
	for i in range(0,reps):
		rand.seed()
		rand1 = rand.random
		rand2 = rand.random
		ra_rand = (ra_UB-ra_LB)*rand1 + ra_LB
		dec_rand = (dec_UB-dec_LB)*rand2 + dec_LB
		likelihoods = OptMatch(ra_rand,dec_rand,err,ra_opt,dec_opt)
		output = np.append(output,likelihoods)
	return output
	
	

#X-Ray
cr_xray = read_file(xray_cat)
# openw, 1, 'Cl0023.xray_phot.soft_hard_full.dat'
# for i=0, n_elements(ra)-1 do printf,1, ra[i], dec[i], flux_softz[i], flux_hardz[i], flux_fullz[i],  netcnts_corr_softz[i], netcnts_corr_hardz[i], netcnts_corr_fullz[i], sig_soft[i], sig_hard[i], sig_full[i], wsig_soft[i], wsig_hard[i], wsig_full[i], wmask[i], wflag[i],format='(1(G," ",G," ",G," ",G," ",G," ",G," ",G," ",G," ",G," ",G," ",G," ",G," ",G," ",G," ",I," ",I," "))'
# close, 1
raX = get_colvals(cr_xray,'col1')
decX = get_colvals(cr_xray,'col2')
fluxX_soft = get_colvals(cr_xray,'col3')
fluxX_hard = get_colvals(cr_xray,'col4')
fluxX_full = get_colvals(cr_xray,'col5')
netcnts_corrX_soft = get_colvals(cr_xray,'col6')
netcnts_corrX_hard = get_colvals(cr_xray,'col7')
netcnts_corrX_full = get_colvals(cr_xray,'col8')
sigX_soft = get_colvals(cr_xray,'col9')
sigX_hard = get_colvals(cr_xray,'col10')
sigX_full = get_colvals(cr_xray,'col11')
wflagX = get_colvals(cr_xray,'col13')
off_axisX = np.zeros(len(raX))
ncnts_corrX = np.zeros((3,len(raX)))
ncnts_corrX[0][:] = netcnts_corrX_soft[:]
ncnts_corrX[1][:] = netcnts_corrX_hard[:]
ncnts_corrX[2][:] = netcnts_corrX_full[:]
xwnetcnts = np.zeros(len(raX))
for i in range(0,len(raX)):
	off_axisX[i] = SphDist(Ra_aim,Dec_aim,raX[i],decX[i])
	xwnetcnts[i] = ncnts_corrX[wflagx][i]

cr_opt = read_file(optical_cat)
ra_opt = get_colvals(cr_opt,'')
dec_opt = get_colvals(cr_opt,'')

# Calculate errors
d_xerr = np.zeros(len(raX))
for i in range(0,len(raX)):
	if xwnetcnts[i] <= 137.816: 10**(0.1145*off_axisX[i] - 0.4958*m.log10(xwnetcnts[i])+0.1932)
	if xwnetcnts[i] > 137.816: d_xerr[i] = 10**(((0.0968*off_axis[i] - 0.2064*m.log10(xwnetcnts[i])-0.4260)))
  	if off_axis[i] >= 15.0: d_xerr[i] = 60.0   # Upper limit is 1' 


#Get source densities
ra_opt_cen = 
dec_opt_cen = 
radius_opt_cen = 
area_opt_cen = m.pi*(radius_opt_cen**2)

d_fromoptcen = np.zeros(len(ra_opt))
for i in range(0,len(ra_opt)):
	d_fromoptcen = SphDist(ra_opt[i],dec_opt[i],ra_opt_cen,dec_opt_cen)

loc_dens_gtMag = np.zeros(len(ra_opt))
for i in range(0,len(ra_opt)):
	ind_local = where(d_fromoptcen,"le",8.0/60.0)
	magI_loc = magI_opt[ind_local]
	ind_lemagI = where(magI_loc,"le",magI_opt[i])
	num_gt_limit = len(ind_lemagI)
	loc_num_dens_gtMag[i] = num_gt_limit/area_opt_cen

#Optical Matching

matched_raX = np.zeros(0)
matched_decX = np.zeros(0)
matched_errX = np.zeros(0)
matched_num_matches = np.zeros(0)
matched_dec_opt1 = np.zeros(0)
matched_dec_opt2 = np.zeros(0)
matched_dec_opt3 = np.zeros(0)
matched_ra_opt1 = np.zeros(0)
matched_ra_opt2 = np.zeros(0)
matched_ra_opt3 = np.zeros(0)
matched_ind_opt1 = np.zeros(0)
matched_ind_opt2 = np.zeros(0)
matched_ind_opt3 = np.zeros(0)
matched_prob_opt1 = np.zeros(0)
matched_prob_opt2 = np.zeros(0)
matched_prob_opt3 = np.zeros(0)
matched_prob_none = np.zeros(0)
unmatched_raX = np.zeros(0)
unmatched_decX = np.zeros(0)
unmatched_errX = np.zeros(0)
unmatched_indX = np.zeros(0)
unmatched_prob = np.zeros(0)

for i in range(0,len(raX)):
	likelihoods,like_inds = OptMatch(raX[i],decX[i],d_xerr[i],ra_opt,dec_opt)
	num_matches = len(likelihoods)
	reliability_ij = np.zeros(len(likelihoods))
	rel_ref = MonteCarlo(reps,ra_UB,ra_LB,dec_UB,dec_LB,d_xerr[i],ra_opt,dec_opt)
	for j in range(0,len(likelihoods)):
		MC_temp = where(rel_ref,"ge",likelihoods[j])
		rel_cnt = len(MC_temp)
		reliability_ij[j] = 1 - rel_cnt/reps
	probability_ij = np.zeros(len(likelihoods))
	for j in range(0,len(likelihoods)):
		probability_ij[j] = reliability_ij[j]
		for k in range(0,len(likelihoods)):
			if k != j: probability_ij[j] *= (1-reliability_ij[k])
	probability_nonej = 1.0
	for j in range(0,len(likelihoods)): probability_nonej *= (1-reliability_ij[j])
	prob_sum = probability_nonej + np.sum(probability_ij)
	probability_ij /= prob_sum
	probability_nonej /= prob_sum
	if num_matches == 1:
		matched_ra_opt2 = np.append(matched_ra_opt2,-1)
		matched_dec_opt2 = np.append(matched_dec_opt2,-1)
		matched_ra_opt3 = np.append(matched_ra_opt3,-1)
		matched_dec_opt3 = np.append(matched_dec_opt3,-1)
		matched_raX = np.append(matched_raX,raX[i])
		matched_decX = np.append(matched_decX,decX[i])
		matched_errX = np.append(matched_errX,d_xerr[i])
		matched_indX = np.append(matched_indX,i)
		matched_ind_opt2 = np.append(matched_indX,-1)
		matched_ind_opt3 = np.append(matched_indX,-1)
		matched_prob_opt2 = np.append(matched_prob_opt2,-1)
		matched_prob_opt3 = np.append(matched_prob_opt3,-1)
		matched_prob_none = np.append(match_prob_none,probability_nonej[i])
		if probability_nonej <= 0.2: 
			matched_ra_opt1 = np.append(matched_ra_opt1,ra_opt[like_inds[0]])
			matched_dec_opt1 = np.append(matched_dec_opt1,dec_opt[like_inds[0]])
			matched_num_matches = np.append(matched_num_matches,1)
			matched_prob_opt1 = np.append(matched_prob_opt1,1-probability_nonej)
			matched_ind_opt1 = np.append(matched_indX,like_inds[0])
			
		else:
			unmatched_raX = np.append(unmatched_raX,raX[i])
			unmatched_decX = np.append(unmatched_decX,decX[i])
			unmatched_errX = np.append(unmatched_errX,d_xerr[i])
			unmatched_indX = np.append(unmatched_indX,i)
			unmatched_prob = np.append(unmatched_prob,probability_nonej[i])
			matched_ra_opt1 = np.append(matched_ra_opt1,-1)
			matched_dec_opt1 = np.append(matched_dec_opt1,-1)
			matched_num_matches = np.append(matched_num_matches,0)
			matched_prob_opt1 = np.append(matched_prob_opt1,-1)
			matched_ind_opt1 = np.append(matched_indX,-1)
	if num_matches > 1:
		if probability_nonej <= 0.2:
			prob_probe_temp = np.zeros(num_matches)
			for k in range(0,num_matches):
				prob_probe_temp = probability_ij[k]/(1 - probability_nonej - probability_ij[k])
			cand_temp_ind = where(prob_probe_temp,"ge",4.0)
			if len(cand_temp_ind) > 1: sys.exit("Probabilities don't sum to 1")
			if len(cand_temp_ind) == 1:
				matched_ra_opt2 = np.append(matched_ra_opt2,-1)
				matched_dec_opt2 = np.append(matched_dec_opt2,-1)
				matched_ra_opt3 = np.append(matched_ra_opt3,-1)
				matched_dec_opt3 = np.append(matched_dec_opt3,-1)
				matched_raX = np.append(matched_raX,raX[i])
				matched_decX = np.append(matched_decX,decX[i])
				matched_errX = np.append(matched_errX,d_xerr[i])
				matched_indX = np.append(matched_indX,i)
				matched_ind_opt2 = np.append(matched_indX,-1)
				matched_ind_opt3 = np.append(matched_indX,-1)
				matched_prob_opt2 = np.append(matched_prob_opt2,-1)
				matched_prob_opt3 = np.append(matched_prob_opt3,-1)
				matched_prob_none = np.append(match_prob_none,probability_nonej[i])
				matched_ra_opt1 = np.append(matched_ra_opt1,ra_opt[like_inds[cand_temp_ind]])
				matched_dec_opt1 = np.append(matched_dec_opt1,dec_opt[like_inds[cand_temp_ind]])
				matched_num_matches = np.append(matched_num_matches,1)
				matched_prob_opt1 = np.append(matched_prob_opt1,probability_ij[cand_temp_ind])
				matched_ind_opt1 = np.append(matched_indX,like_inds[cand_temp_ind])
			if len(cand_temp_ind) == 0:
				if num_matches < 1: sys.exit("Inconsistency: 001")
				test_sig_cands_ind = where(probability_ij,"ge",0.2)
				if len(test_sig_cands_ind) > 3: sys.exit("More than 3 matches")
				if len(test_sig_cands_ind) > 0:
					matched_prob_none = np.append(match_prob_none,probability_nonej[i])
					matched_ra_opt1 = np.append(matched_ra_opt1,ra_opt[like_inds[test_sig_cands_ind[0]]])
					matched_dec_opt1 = np.append(matched_dec_opt1,dec_opt[like_inds[test_sig_cands_ind[0]]])
					matched_prob_opt1 = np.append(matched_prob_opt1,probability_ij[test_sig_cands_ind[0]])
					matched_ind_opt1 = np.append(matched_indX,like_inds[test_sig_cands_ind[0]])
					matched_raX = np.append(matched_raX,raX[i])
					matched_decX = np.append(matched_decX,decX[i])
					matched_errX = np.append(matched_errX,d_xerr[i])
					matched_indX = np.append(matched_indX,i)
					if len(test_sig_cands_ind) == 1:
						matched_len(test_sig_cands_ind) = np.append(matched_len(test_sig_cands_ind),1)
						matched_ind_opt2 = np.append(matched_indX,-1)
						matched_ind_opt3 = np.append(matched_indX,-1)
						matched_prob_opt2 = np.append(matched_prob_opt2,-1)
						matched_prob_opt3 = np.append(matched_prob_opt3,-1)
						matched_ra_opt2 = np.append(matched_ra_opt2,-1)
						matched_dec_opt2 = np.append(matched_dec_opt2,-1)
						matched_ra_opt3 = np.append(matched_ra_opt3,-1)
						matched_dec_opt3 = np.append(matched_dec_opt3,-1)
					if len(test_sig_cands_ind) == 2:
						matched_len(test_sig_cands_ind) = np.append(matched_len(test_sig_cands_ind),2)
						matched_ra_opt2 = np.append(matched_ra_opt1,ra_opt[like_inds[test_sig_cands_ind[1]]])
						matched_dec_opt2 = np.append(matched_dec_opt1,dec_opt[like_inds[test_sig_cands_ind[1]]])
						matched_prob_opt2 = np.append(matched_prob_opt1,probability_ij[test_sig_cands_ind[1]])
						matched_ind_opt2 = np.append(matched_indX,like_inds[test_sig_cands_ind[1]])
						matched_ra_opt3 = np.append(matched_ra_opt3,-1)
						matched_dec_opt3 = np.append(matched_dec_opt3,-1)
						matched_prob_opt3 = np.append(matched_prob_opt3,-1)
						matched_ind_opt3 = np.append(matched_indX,-1)
					if len(test_sig_cands_ind) == 3:
						matched_len(test_sig_cands_ind) = np.append(matched_len(test_sig_cands_ind),3)
						matched_ra_opt3 = np.append(matched_ra_opt1,ra_opt[like_inds[test_sig_cands_ind[2]]])
						matched_dec_opt3 = np.append(matched_dec_opt1,dec_opt[like_inds[test_sig_cands_ind[2]]])
						matched_prob_opt3 = np.append(matched_prob_opt1,probability_ij[test_sig_cands_ind[2]])
						matched_ind_opt3 = np.append(matched_indX,like_inds[test_sig_cands_ind[2]])
				else:
					unmatched_raX = np.append(unmatched_raX,raX[i])
					unmatched_decX = np.append(unmatched_decX,decX[i])
					unmatched_errX = np.append(unmatched_errX,d_xerr[i])
					unmatched_indX = np.append(unmatched_indX,i)
					unmatched_prob = np.append(unmatched_prob,probability_nonej[i])
					matched_ra_opt1 = np.append(matched_ra_opt1,-1)
					matched_dec_opt1 = np.append(matched_dec_opt1,-1)
					matched_num_matches = np.append(matched_num_matches,0)
					matched_prob_opt1 = np.append(matched_prob_opt1,-1)
					matched_ind_opt1 = np.append(matched_indX,-1)
					matched_ra_opt2 = np.append(matched_ra_opt2,-1)
					matched_dec_opt2 = np.append(matched_dec_opt2,-1)
					matched_ra_opt3 = np.append(matched_ra_opt3,-1)
					matched_dec_opt3 = np.append(matched_dec_opt3,-1)
					matched_ind_opt2 = np.append(matched_indX,-1)
					matched_ind_opt3 = np.append(matched_indX,-1)
					matched_prob_opt2 = np.append(matched_prob_opt2,-1)
					matched_prob_opt3 = np.append(matched_prob_opt3,-1)
					matched_raX = np.append(matched_raX,raX[i])
					matched_decX = np.append(matched_decX,decX[i])
					matched_errX = np.append(matched_errX,d_xerr[i])
					matched_indX = np.append(matched_indX,i)
					matched_prob_none = np.append(match_prob_none,probability_nonej[i])
			
					
		else:
			unmatched_raX = np.append(unmatched_raX,raX[i])
			unmatched_decX = np.append(unmatched_decX,decX[i])
			unmatched_errX = np.append(unmatched_errX,d_xerr[i])
			unmatched_indX = np.append(unmatched_indX,i)
			unmatched_prob = np.append(unmatched_prob,probability_nonej[i])
			matched_ra_opt1 = np.append(matched_ra_opt1,-1)
			matched_dec_opt1 = np.append(matched_dec_opt1,-1)
			matched_num_matches = np.append(matched_num_matches,0)
			matched_prob_opt1 = np.append(matched_prob_opt1,-1)
			matched_ind_opt1 = np.append(matched_indX,-1)
			matched_ra_opt2 = np.append(matched_ra_opt2,-1)
			matched_dec_opt2 = np.append(matched_dec_opt2,-1)
			matched_ra_opt3 = np.append(matched_ra_opt3,-1)
			matched_dec_opt3 = np.append(matched_dec_opt3,-1)
			matched_ind_opt2 = np.append(matched_indX,-1)
			matched_ind_opt3 = np.append(matched_indX,-1)
			matched_prob_opt2 = np.append(matched_prob_opt2,-1)
			matched_prob_opt3 = np.append(matched_prob_opt3,-1)
		       	matched_raX = np.append(matched_raX,raX[i])
		       	matched_decX = np.append(matched_decX,decX[i])
		       	matched_errX = np.append(matched_errX,d_xerr[i])
		       	matched_indX = np.append(matched_indX,i)
			matched_prob_none = np.append(match_prob_none,probability_nonej[i])
			
