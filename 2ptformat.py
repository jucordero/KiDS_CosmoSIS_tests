from astropy.io import fits
import numpy as np

def angbin(angle):
	switcher = {
		0.71336: 1,
		1.45210: 2,
		2.95582: 3,
		6.01675: 4,
		12.24745: 5,
		24.93039: 6,
		50.74725: 7,
		103.29898: 8,
		210.27107: 9,

		0.713365: 1,
		1.452096: 2,
		2.955825: 3,
		6.016752: 4,
		12.24745: 5,
		24.93039: 6,
		50.74726: 7,
		103.299: 8,
		210.271: 9,
	}
	return switcher.get(angle, 99999)

def bin_num(p,m):
	if p == 1:
		if m == 1:
			return 0
		if m == 2:
			return 1
		if m == 3:
			return 2
		if m == 4:
			return 3
	if p == 2:
		if m == 2:
			return 4
		if m == 3:
			return 5
		if m == 4:
			return 6
	if p == 3:
		if m == 3:
			return 7
		if m == 4:
			return 8
	if p == 4:
		if m == 4:
			return 9

	return 99999
		

#-----------------------------Covariance Matrix----------------------------------
cov = np.genfromtxt('./KiDS_Data/COV_MAT/Cov_mat_all_scales.txt')
COVMAT = np.zeros((180, 180))

for i in range(len(cov)):

	COVMAT[int(cov[i, 2]*90 + 10*(angbin(cov[i, 3])-1)+bin_num(cov[i, 0], cov[i, 1])), 
	int(cov[i, 6]*90 + 10*(angbin(cov[i, 7])-1)+bin_num(cov[i, 4], cov[i, 5]))] = cov[i, 8]+cov[i, 9]+cov[i, 10]

	COVMAT[int(cov[i, 6]*90 + 10*(angbin(cov[i, 7])-1)+bin_num(cov[i, 4], cov[i, 5])), 
	int(cov[i, 2]*90 + 10*(angbin(cov[i, 3])-1)+bin_num(cov[i, 0], cov[i, 1]))] = cov[i, 8]+cov[i, 9]+cov[i, 10]

hdu_covmat = fits.ImageHDU(COVMAT)
hdu_covmat.header['EXTNAME'] = 'COVMAT'
hdu_covmat.header['NAME_0'] = 'xi_plus'
hdu_covmat.header['NAME_1'] = 'xi_minus'
hdu_covmat.header['STRT_0'] = 0
hdu_covmat.header['STRT_1'] = 90
hdu_covmat.header['COVDATA'] = True
hdu_covmat.header['LEN_0'] = 90
hdu_covmat.header['LEN_1'] = 90


#--------------------------------------------------------------------------------
#-----------------------------Correlation Function-------------------------------
XI_PLUS = np.zeros((90,5))
XI_MINUS = np.zeros((90,5))

index = 0
for i in range(1,5):
	for j in range (i,5):
		currxi = np.genfromtxt('./CorrelationFunctions/xi_' + str(i) + '-' + str(j))
		for k in range(len(currxi)):			
			XI_PLUS[k*10+index,:] = i, j, angbin(currxi[k,0]), currxi[k,1], currxi[k,0]
			XI_MINUS[k*10+index,:] = i, j, angbin(currxi[k,0]), currxi[k,2], currxi[k,0]
		index += 1
		
col1_p = fits.Column(name = 'BIN1', format = 'K', array = XI_PLUS[:,0])
col2_p = fits.Column(name = 'BIN2', format = 'K', array = XI_PLUS[:,1])
col3_p = fits.Column(name = 'ANGBIN', format = 'K', array = XI_PLUS[:,2])
col4_p = fits.Column(name = 'VALUE', format = 'D', array = XI_PLUS[:,3])
col5_p = fits.Column(name = 'ANG', format = 'D', array = XI_PLUS[:,4])

col1_m = fits.Column(name = 'BIN1', format = 'K', array = XI_MINUS[:,0])
col2_m = fits.Column(name = 'BIN2', format = 'K', array = XI_MINUS[:,1])
col3_m = fits.Column(name = 'ANGBIN', format = 'K', array = XI_MINUS[:,2])
col4_m = fits.Column(name = 'VALUE', format = 'D', array = XI_MINUS[:,3])
col5_m = fits.Column(name = 'ANG', format = 'D', array = XI_MINUS[:,4])

cols_p = fits.ColDefs([col1_p,col2_p,col3_p,col4_p,col5_p])
cols_m = fits.ColDefs([col1_m,col2_m,col3_m,col4_m,col5_m])

hdu_xi_plus = fits.BinTableHDU.from_columns(cols_p)
hdu_xi_minus = fits.BinTableHDU.from_columns(cols_m)

hdu_xi_plus.header['EXTNAME'] = 'xi_plus'
hdu_xi_plus.header['QUANT1'] = 'G+R'
hdu_xi_plus.header['QUANT2'] = 'G+R'
hdu_xi_plus.header['KERNEL_1'] = 'NZ_CC'
hdu_xi_plus.header['KERNEL_2'] = 'NZ_CC'
hdu_xi_plus.header['2PTDATA'] = True
hdu_xi_plus.header['WINDOWS'] = 'SAMPLE'
hdu_xi_plus.header['N_ZBIN1'] = 4
hdu_xi_plus.header['N_ZBIN2'] = 4
hdu_xi_plus.header['N_ANG'] = 7
hdu_xi_plus.header['TUNIT5'] = 'arcmin'

hdu_xi_minus.header['EXTNAME'] = 'xi_minus'
hdu_xi_minus.header['QUANT1'] = 'G+R'
hdu_xi_minus.header['QUANT2'] = 'G+R'
hdu_xi_minus.header['KERNEL_1'] = 'NZ_CC'
hdu_xi_minus.header['KERNEL_2'] = 'NZ_CC'
hdu_xi_minus.header['2PTDATA'] = True
hdu_xi_minus.header['WINDOWS'] = 'SAMPLE'
hdu_xi_minus.header['N_ZBIN1'] = 4
hdu_xi_minus.header['N_ZBIN2'] = 4
hdu_xi_minus.header['N_ANG'] = 6
hdu_xi_minus.header['TUNIT5'] = 'arcmin'

#--------------------------------------------------------------------------------
#-----------------------------Redshift number density distribution---------------

nz_dir = np.zeros((70,7))
nz_cc = np.zeros((40,7))

for i in range(4):
	ccfile = np.genfromtxt('./KiDS_Data/Nz_CC/Nz_CC_z0.' + str(1+i*2) + 't0.' + str(3+i*2) + '.asc')
	dirfile = np.genfromtxt('./KiDS_Data/Nz_DIR/Nz_DIR_Mean/Nz_DIR_z0.' + str(1+i*2) + 't0.' + str(3+i*2) + '.asc')
	for j in range(len(ccfile)):
		nz_cc[j,3+i] = ccfile[j,1]
		if i == 3:
			nz_cc[j,0] = ccfile[j,0]
			if j < len(ccfile)-1:
				nz_cc[j,2] = ccfile[j+1,0]
			else:
				nz_cc[j,2] = ccfile[j,0]+0.05
			nz_cc[j,1] = 0.5*(nz_cc[j,0] + nz_cc[j,2])

	for j in range(len(dirfile)):
		nz_dir[j,3+i] = dirfile[j,1]
		if i == 3:
			nz_dir[j,0] = dirfile[j,0]
			if j < len(dirfile)-1:
				nz_dir[j,2] = dirfile[j+1,0]
			else:
				nz_dir[j,2] = dirfile[j,0]+0.05
			nz_dir[j,1] = 0.5*(nz_dir[j,0] + nz_dir[j,2])


col1_cc = fits.Column(name = 'Z_LOW', format = 'D', array = nz_cc[:,0])
col2_cc = fits.Column(name = 'Z_MID', format = 'D', array = nz_cc[:,1])
col3_cc = fits.Column(name = 'Z_HIGH', format = 'D', array = nz_cc[:,2])
col4_cc = fits.Column(name = 'BIN1', format = 'D', array = nz_cc[:,3])
col5_cc = fits.Column(name = 'BIN2', format = 'D', array = nz_cc[:,4])
col6_cc = fits.Column(name = 'BIN3', format = 'D', array = nz_cc[:,5])
col7_cc = fits.Column(name = 'BIN4', format = 'D', array = nz_cc[:,6])

col1_dir = fits.Column(name = 'Z_LOW', format = 'D', array = nz_dir[:,0])
col2_dir = fits.Column(name = 'Z_MID', format = 'D', array = nz_dir[:,1])
col3_dir = fits.Column(name = 'Z_HIGH', format = 'D', array = nz_dir[:,2])
col4_dir = fits.Column(name = 'BIN1', format = 'D', array = nz_dir[:,3])
col5_dir = fits.Column(name = 'BIN2', format = 'D', array = nz_dir[:,4])
col6_dir = fits.Column(name = 'BIN3', format = 'D', array = nz_dir[:,5])
col7_dir = fits.Column(name = 'BIN4', format = 'D', array = nz_dir[:,6])

cols_cc = fits.ColDefs([col1_cc, col2_cc, col3_cc, col4_cc, col5_cc, col6_cc, col7_cc])
cols_dir = fits.ColDefs([col1_dir, col2_dir, col3_dir, col4_dir, col5_dir, col6_dir, col7_dir])

hdu_nz_cc = fits.BinTableHDU.from_columns(cols_cc)
hdu_nz_dir = fits.BinTableHDU.from_columns(cols_dir)


hdu_nz_cc.header['EXTNAME'] = 'NZ_CC'
hdu_nz_cc.header['NZDATA'] = True
hdu_nz_cc.header['NBIN'] = 4
hdu_nz_cc.header['NZ'] = 40

hdu_nz_dir.header['EXTNAME'] ='NZ_DIR'
hdu_nz_dir.header['NZDATA'] = True
hdu_nz_dir.header['NBIN'] = 4
hdu_nz_dir.header['NZ'] = 70



#--------------------------------------------------------------------------------
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, hdu_covmat, hdu_xi_plus, hdu_xi_minus, hdu_nz_cc, hdu_nz_dir])
thdulist.writeto('KiDS_like.fits', clobber = True)





















