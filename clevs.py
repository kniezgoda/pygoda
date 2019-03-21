import numpy as np
import matplotlib.pyplot as plt
from pygoda import niceClev


def getClev(var, season, data = None):
	prettyClev = 19
	if data is not None:
		prettyClev = niceClev(data)
	clevs = { \
				"PRECT_H2O" : \
					{"ANN_" : np.linspace(0,16,17), \
					"ANN_NA" : np.linspace(0,9,19), \
					"ANN_MC" : np.linspace(0,9,19), \
					"ANN_IM" : np.linspace(0,9,19), \
					"ANN_EP" : np.linspace(0,9,19), \
					\
					"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJF_" : np.linspace(0,16,17), \
					"DJF_NA" : np.linspace(0,16,17), \
					"DJF_MC" : np.linspace(0,16,17), \
					"DJF_IM" : np.linspace(0,9,19), \
					"DJF_EP" : np.linspace(0,9,19), \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJA_" : np.linspace(0,30,21), \
					"JJA_NA" : np.linspace(0,30,21), \
					"JJA_MC" : np.linspace(0,30,21), \
					"JJA_IM" : np.linspace(0,9,19), \
					"JJA_EP" : np.linspace(0,9,19), \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \
				"PRECT" : \
					{"ANN_" : np.linspace(0,16,17), \
					"ANN_NA" : np.linspace(0,9,19), \
					"ANN_MC" : np.linspace(0,9,19), \
					"ANN_IM" : np.linspace(0,9,19), \
					"ANN_EP" : np.linspace(0,9,19), \
					\
					"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJF_" : np.linspace(0,16,17), \
					"DJF_NA" : np.linspace(0,16,17), \
					"DJF_MC" : np.linspace(0,16,17), \
					"DJF_IM" : np.linspace(0,9,19), \
					"DJF_EP" : np.linspace(0,9,19), \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJA_" : np.linspace(0,30,21), \
					"JJA_NA" : np.linspace(0,30,21), \
					"JJA_MC" : np.linspace(0,30,21), \
					"JJA_IM" : np.linspace(0,9,19), \
					"JJA_EP" : np.linspace(0,9,19), \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"fluxDelta" : \
					{"ANN_" : np.linspace(0,16,17), \
					"ANN_NA" : np.linspace(0,9,19), \
					"ANN_MC" : np.linspace(0,9,19), \
					"ANN_IM" : np.linspace(0,9,19), \
					"ANN_EP" : np.linspace(0,9,19), \
					\
					"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJF_" : np.linspace(0,16,17), \
					"DJF_NA" : np.linspace(0,16,17), \
					"DJF_MC" : np.linspace(0,16,17), \
					"DJF_IM" : np.linspace(0,9,19), \
					"DJF_EP" : np.linspace(0,9,19), \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJA_" : np.linspace(0,30,21), \
					"JJA_NA" : np.linspace(0,30,21), \
					"JJA_MC" : np.linspace(-30,30,21), \
					"JJA_IM" : np.linspace(0,9,19), \
					"JJA_EP" : np.linspace(0,9,19), \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : np.linspace(-30,30,21), \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"PRECT_d18O" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"OMEGA850" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"OMEGA500" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"OMEGA200" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"V850" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"V500" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"V200" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"VQ850" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"VQ500" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"VQ200" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"U850" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"U500" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"U200" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"UQ850" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"UQ500" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}, \

				"UQ200" : \
					{"ANNdiff_" : prettyClev, \
					"ANNdiff_NA" : prettyClev, \
					"ANNdiff_MC" : prettyClev, \
					"ANNdiff_IM" : prettyClev, \
					"ANNdiff_EP" : prettyClev, \
					\
					"DJFdiff_" : prettyClev, \
					"DJFdiff_NA" : prettyClev, \
					"DJFdiff_MC" : prettyClev, \
					"DJFdiff_IM" : prettyClev, \
					"DJFdiff_EP" : prettyClev, \
					\
					"JJAdiff_" : prettyClev, \
					"JJAdiff_NA" : prettyClev, \
					"JJAdiff_MC" : prettyClev, \
					"JJAdiff_IM" : prettyClev, \
					"JJAdiff_EP" : prettyClev, \
					}}
	return clevs.get(var, dict()).get(season, prettyClev)

def getCmap(var, which):
	cmap = { \
				"PRECT_H2O" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"PRECT" :  { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"PRECT_d18O" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"Q" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"V" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"VT" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"VQ" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"U" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"UT" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"UQ" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"T" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"OMEGA" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"Z3" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"dDV" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"d18OV" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}, \

				"dxsV" : { \
					"cmap" : plt.cm.coolwarm, \
					"diffcmap" : plt.cm.RdBu_r \
					}}
	return cmap.get(var, dict()).get(which, plt.cm.coolwarm)