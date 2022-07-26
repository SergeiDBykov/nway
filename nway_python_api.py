print('This is NWAY Python API')


import sys
import numpy
from numpy import log10, pi, exp, logical_and
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import argparse
import nwaylib 
import nwaylib.progress as progress
print('nwaylib file', nwaylib.__file__)
import nwaylib.fastskymatch as match



def table_from_fits(fitsname, poserr_value=None, area=None, magnitude_columns=[]):
	fits_table = pyfits.open(fitsname)[1]
	table_name = fits_table.name
	ra = fits_table.data['RA']
	dec = fits_table.data['DEC']
	if 'pos_err' in fits_table.data.columns.names:
		poserr = fits_table.data['pos_err']
	else:
		assert poserr_value is not None, ('"pos_err" column not found in file "%s", and no poserr_value passed' % fitsname)
		poserr = poserr_value * numpy.ones(len(ra))
	if area is None:
		area = fits_table.header['SKYAREA'] * 1.0 
	
	# magnitude columns
	mags = []
	maghists = []
	magnames = []
	#for mag in magnitude_columns:
	for col_name, magfile in magnitude_columns:
		assert col_name in fits_table.data.dtype.names
		
		mag_all = fits_table.data[col_name]
		# mark -99 as undefined
		mag_all[mag_all == -99] = numpy.nan
		
		mags.append(mag_all)
		magnames.append(col_name)
		if magfile == 'auto':
			maghists.append(None)
		else:
			bins_lo, bins_hi, hist_sel, hist_all = numpy.loadtxt(magfile).transpose()
			maghists.append((bins_lo, bins_hi, hist_sel, hist_all))
	
	return dict(name=table_name, ra=ra, dec=dec, error=poserr, area=area, mags=mags, maghists=maghists, magnames=magnames)
	# area in square degrees
	# error in arcsec
	# ra/dec in degrees
	# mag: column of something
	# maghists: either (bin, sel, all) tuple or None (for auto)



def create_shifted_catalogue(inputfile: str, outputfile:str,
shift_dec: float = 0, shift_ra: float = 0,  radius: float = 0):


    filename = inputfile
    outfile = outputfile
    print('opening', filename)
    inputfitsfile = pyfits.open(filename)
    header_hdu = inputfitsfile[0]
    table = inputfitsfile[1].data

    if shift_ra==0 and shift_dec==0:
        raise ValueError('ERROR: You have to set either shift-ra or shift-dec to non-zero')

    ra_key  = match.get_tablekeys(table, 'RA')
    print('    using RA  column: %s' % ra_key)
    dec_key = match.get_tablekeys(table, 'DEC')
    print('    using DEC column: %s' % dec_key)

    ra_orig = table[ra_key]
    dec_orig = table[dec_key]

    ra = ra_orig + shift_ra / 60. / 60
    dec = dec_orig + shift_dec / 60. / 60

    # for each of them, check that there is no collision
    excluded = []

    pbar = progress.bar(ndigits=6)
    for i, (ra_i, dec_i) in pbar(list(enumerate(zip(ra, dec)))):
        d = match.dist((ra_i, dec_i), (ra_orig, dec_orig))
        excluded.append((d * 60 * 60 < radius).any())

    excluded = numpy.array(excluded)
    print('removed %d sources which collide with original positions' % (excluded.sum()))

    table[ra_key] = ra
    table[dec_key] = dec

    newcolumns = []
    for col in table.columns:
        newcolumns.append(pyfits.Column(name=col.name, format=col.format, array=col.array[~excluded]))

    tbhdu = match.fits_from_columns(newcolumns)
    print('writing "%s" (%d rows)' % (outfile, len(tbhdu.data)))
    for k, v in inputfitsfile[1].header.items():
        if k not in tbhdu.header:
            tbhdu.header[k] = v

    hdulist = pyfits.HDUList([header_hdu, tbhdu])
    hdulist.writeto(outfile, **progress.kwargs_overwrite_true)




def calibrate_cutoff(realfile_table, fakefile_table):


    data = realfile_table
    mask = data['ncat'] == 1
    p_any0 = data['prob_has_match'][mask]
    mask = data['match_flag'] == 1
    p_any = data['prob_has_match'][mask]
    p_i = data['prob_this_match'][mask]
    
    # mask = data['ncat'] == 1
    # p_any0 = data['p_any'][mask]
    # mask = data['match_flag'] == 1
    # p_any = data['p_any'][mask]
    # p_i = data['p_i'][mask]

    data = fakefile_table
    mask = data['ncat'] == 1
    p_any0_offset = data['prob_has_match'][mask]
    mask = data['match_flag'] == 1
    p_any_offset = data['prob_has_match'][mask]
    p_i_offset = data['prob_this_match'][mask]

    # mask = data['ncat'] == 1
    # p_any0_offset = data['p_any'][mask]
    # mask = data['match_flag'] == 1
    # p_any_offset = data['p_any'][mask]
    # p_i_offset = data['p_i'][mask]

    cutoffs = numpy.linspace(0, 1, 101)

    efficiency = [(p_any0 > cutoff).mean() for cutoff in cutoffs]
    error_rate = [(p_any0_offset > cutoff).mean() for cutoff in cutoffs]



    for rate in [0.01,0.03,0.05,0.1]:
        print()
        # find where error-rate is < rate
        mask = numpy.array(error_rate) < rate
        if not mask.any():
            print('A false detection rate of <%d%% is not possible.' % (rate*100))
        else:
            i = numpy.min(numpy.where(mask)[0])
            #print(i)
            print('For a false detection rate of <%d%%' % (rate*100))
            print('--> use only counterparts with p_any>%.2f (%.2f%% of matches)' % (cutoffs[i], efficiency[i]*100))



    plt.figure(figsize=(5,5))
    plt.plot(cutoffs, efficiency, '-', color='k', label='completeness (selection efficiency)', lw=2)
    plt.plot(cutoffs, purity:=1-numpy.asanyarray(error_rate), '--', color='r', label='purity (1-false selection rate)')
    #find the intersection of the two curves: efficiency and purity
    cutoff_intersection = numpy.argmin(numpy.abs(efficiency - purity))
    print('The efficiency is %.2f%%' % (tmp_val:=efficiency[cutoff_intersection]*100))
    print('The purity is  %.2f%%' % (purity[cutoff_intersection]*100))
    plt.axvline(cutoffs[cutoff_intersection], color='k', ls='--', label=f'purity=completeness={tmp_val:.2g}%')


    plt.ylabel('fraction(>p_any)')
    plt.xlabel('p_any')
    plt.legend(loc='lower left', prop=dict(size=10))
    #plt.savefig(realfile + '_p_any_cutoffquality.pdf', bbox_inches='tight')
    #plt.savefig(realfile + '_p_any_cutoffquality.png', bbox_inches='tight')
    plt.show()
    #print('created plot "%s_p_any_cutoffquality.pdf"' % realfile)



    plt.figure(figsize=(5,5))
    plt.plot(cutoffs, efficiency, '-', color='k', label='selection efficiency', lw=2)
    plt.plot(cutoffs, error_rate, '--', color='r', label='false selection rate')

    plt.ylabel('fraction(>p_any)')
    plt.xlabel('p_any')
    plt.legend(loc='lower left', prop=dict(size=10))
    #plt.savefig(realfile + '_p_any_cutoffquality.pdf', bbox_inches='tight')
    #plt.savefig(realfile + '_p_any_cutoffquality.png', bbox_inches='tight')
    plt.show()
    #print('created plot "%s_p_any_cutoffquality.pdf"' % realfile)

    plt.figure(figsize=(5,5))
    plt.plot(p_any, p_i, '. ', color='r', label='real')
    plt.plot(p_any_offset, p_i_offset, '. ', color='gray', label='offset')

    plt.ylabel('p_i')
    plt.xlabel('p_any')
    plt.legend(loc='lower left', prop=dict(size=10))
    #plt.savefig(realfile + '_p_any_p_i.pdf', bbox_inches='tight')
    #plt.savefig(realfile + '_p_any_p_i.png', bbox_inches='tight')
    plt.show()
    #print('created plot "%s_p_any_p_i.pdf"' % args.realfile)

    return cutoffs, efficiency, error_rate



# def explain(data, id, primary_id_col, cols_ra = ['RA'], cols_dec=['DEC'], cols_err = ['pos_err']):

#     mask = data[primary_id_col] == id

#     assert mask.sum() != 0, print('ERROR: ID not found. Was searching for %s == %s' % (primary_id_col, id))


#     center_ra = data[cols_ra[0]][mask][0]
#     center_dec = data[cols_dec[0]][mask][0]
#     p_any = data['p_any'][mask][0]




#     print('NWAY results for Source %s:' % id)
#     print()
#     if p_any > 0.8:
#         print('This source probably has a counterpart (p_any=%.2f)' % p_any)
#     elif p_any < 0.1:
#         print('This source probably does not a counterpart (p_any=%.2f)' % p_any)
#     else:
#         print('It is uncertain if this source has a counterpart (p_any=%.2f)' % p_any)
#     print()
#     print("Assuming it has a counterpart, we have the following possible associations:")
#     print()

#     j_option = 0
#     def print_option(name, i):
#         global j_option
#         j_option += 1
#         matchflag = data['match_flag'][mask][i]
#         matchflagstars = '**' if matchflag == 1 else ('*' if matchflag==2 else '')
#         if p_any < 0.1:
#             matchflagstars = ''
#         if matchflag == 0:
#             print('Association %d: probability p_i=%.2f ' % (j_option, data['p_i'][mask][i]))
#         else:
#             print('Association %d%s[match_flag==%d]: probability p_i=%.2f ' % (j_option, matchflagstars, matchflag, data['p_i'][mask][i]))
#         print('     Involved catalogues:  %s ' % (name))
#         for col in header['BIASING'].split(', '):
#             if col.strip() == '': continue
#             bias = data['bias_' + col][mask][i]
#             if bias >= 2:
#                 print('     prior %-15s increased the probability (bias_%s=%.2f)' % (col, col, bias))
#             elif bias <= 0.5:
#                 print('     prior %-15s decreased the probability (bias_%s=%.2f)' % (col, col, bias))
#         print()

#     def convx(ra):
#         return (ra - center_ra) * 3600
#     def convy(dec):
#         return (dec - center_dec) * 3600
#     def converr(err):
#         return err * 3600

#     markers = ['x', '+', '^', '<', '>', 'v', 'p'] * 10
#     colors = ['b', 'c', 'g', 'r', 'k', 'brown'] * 10
#     all_options = []
#     ii = numpy.argsort(data['p_i'][mask])[::-1]
#     for col_ra, col_dec, col_err, marker, color in zip(cols_ra, cols_dec, cols_err, markers, colors):
#         tblname, err = col_err.split('_', 1)
#         if err.startswith(':'):
#             err_components = err[1:].split(':')
#             if len(err_components) == 1:
#                 ra_err = data[tblname + '_' + err_components[0]][mask]
#                 dec_err = ra_err
#                 pa_err = numpy.zeros(ra_err.shape)
#             elif len(err_components) == 2:
#                 ra_err = data[tblname + '_' + err_components[0]][mask]
#                 dec_err = data[tblname + '_' + err_components[1]][mask]
#                 pa_err = numpy.zeros(ra_err.shape)
#             elif len(err_components) == 3:
#                 ra_err = data[tblname + '_' + err_components[0]][mask]
#                 dec_err = data[tblname + '_' + err_components[1]][mask]
#                 pa_err = data[tblname + '_' + err_components[2]][mask]
#             else:
#                 assert len(err_components) == 3, err_components
#         else:
#             ra_err = numpy.ones(mask.sum()) * float(err)
#             dec_err = numpy.ones(mask.sum()) * float(err)
#             pa_err = numpy.zeros(ra_err.shape)
#         pos = set(zip(data[col_ra][mask], data[col_dec][mask], ra_err, dec_err, pa_err))
#         ras = numpy.array([ra for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99])
#         decs = numpy.array([dec for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99])
#         ra_errs = numpy.array([ra_err for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99]) / 60. / 60.
#         dec_errs = numpy.array([dec_err for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99]) / 60. / 60.
#         pa_errs = numpy.array([pa_err for ra, dec, ra_err, dec_err, pa_err in pos if ra != -99])
#         r,  = plt.plot(convx(ras), convy(decs), marker=marker, mec=color, mfc='None', ms=8, mew=2, ls=' ', label='%s %s' % (col_ra, col_dec))
#         patches = [
#             Ellipse((convx(ra), convy(dec)),
#             2 * converr(ra_err),
#             2 * converr(dec_err),
#             angle=90 - pa_err)
#                 for ra, dec, ra_err, dec_err, pa_err in zip(ras, decs, ra_errs, dec_errs, pa_errs)]
#         p = PatchCollection(patches)
#         p.set_facecolor('None')
#         p.set_edgecolor(color)
#         plt.gca().add_collection(p)
#         options = [(-99, -99)]
#         for i in ii:
#             ra, dec = data[col_ra][mask][i], data[col_dec][mask][i]
#             if (ra, dec) not in options:
#                 options.append((ra, dec))
#         all_options.append(options)

#     def graph_make(all_options, highlight=None):
#         for j, options in enumerate(all_options):
#             if j != 0:
#                 plt.plot(j, 0, marker='x', color='gray')
#             for k in range(len(options)-1):
#                 plt.plot(j, -(k + 1), marker='o', color='k')

#     def graph_highlight(all_options, selected):
#         x = numpy.arange(len(all_options))
#         plt.plot(x, selected, '-', color='k')


#     out_options = []
#     outfilename = '%s_explain_%s_options.pdf' % (args.matchcatalogue, args.id)
#     pp = PdfPages(outfilename)
#     maxy = max([len(o)-1 for o in all_options])
#     maxx = len(all_options)-1
#     for i in ii:
#         plt.figure(figsize=(3+maxx, 3))
#         plt.axis('off')
#         graph_make(all_options)
#         j = []
#         name = []
#         for col_ra, col_dec, options, tablename in zip(cols_ra, cols_dec, all_options, tablenames):
#             radec = data[col_ra][mask][i], data[col_dec][mask][i]
#             plt.text(len(j), 0.1, col_ra + '\n' + col_dec, 
#                 rotation=90, size=6, ha='center', va='bottom')
#             k = options.index(radec)
#             j.append(-k)
#             if k == 0:
#                 name.append("")
#             elif k == 1:
#                 name.append("%s" % tablename)
#         print_option('-'.join(name), i)
#         graph_highlight(all_options, j)
#         plt.text(-0.1, -1, 'p_i=%.2f' % data['p_i'][mask][i], ha='right', va='center')
#         plt.text(maxx + 0.1, 0, '$\leftarrow$ absent', ha='left', va='center')
#         plt.ylim(-maxy-0.5, 0.5)
#         plt.xlim(-0.5, maxx+0.5)
#         plt.savefig(pp, format='pdf', bbox_inches='tight')
#         plt.close()
#     pp.close()
#     print('plotting to %s' % outfilename)

#     # go through each association and highlight
#     for j, i in enumerate(numpy.argsort(data['p_single'][mask])[::-1][:3]):
#         ras = []
#         decs = []
#         for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
#             if data[col_ra][mask][i] == -99:
#                 continue
#             ra = data[col_ra][mask][i]
#             dec = data[col_dec][mask][i]
#             ras.append(ra)
#             decs.append(dec)
        
#         plt.plot(convx(ras), convy(decs), '-', lw=(3-j), label='top %s by distance (p_single=%.2f, %d cat.)' % (j+1, data['p_single'][mask][i], data['ncat'][mask][i]), color='y')

#     mask2 = numpy.logical_and(mask, data['match_flag'] == 1)
#     ras = []
#     decs = []
#     first = True
#     for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
#         ra = float(data[col_ra][mask2])
#         dec = float(data[col_dec][mask2])
#         if ra == -99:
#             continue
#         if not first:
#             plt.text(convx(ra), convy(dec), ' 1', va='top', ha='left', alpha=0.5, size=16, fontweight='bold')
#         first = False
#         if ra == -99:
#             continue
#         ras.append(ra)
#         decs.append(dec)
#     plt.plot(convx(ras), convy(decs), '-', lw=1.7, label='p_i=%.2f (match_flag=1)' % (float(data['p_i'][mask2])), color='orange')

#     mask2 = numpy.logical_and(mask, data['match_flag'] == 2)
#     for i in numpy.where(mask2)[0]:
#         ras = []
#         decs = []
#         first = True
#         for col_ra, col_dec, marker in zip(cols_ra, cols_dec, markers):
#             ra = data[col_ra][i]
#             dec = data[col_dec][i]
#             if ra == -99: 
#                 continue
#             if not first:
#                 plt.text(convx(ra), convy(dec), 
#                     ' 2', va='top', ha='left', alpha=0.5, 
#                     size=16, fontweight='bold')
#             first = False
#             ras.append(ra)
#             decs.append(dec)
#         plt.plot(convx(ras), convy(decs), '-', lw=0.5, label='p_i=%.2f (match_flag=2)' % (data['p_i'][i]), color='yellow')

#     plt.xlabel('$\Delta$RA [arcsec]')
#     plt.ylabel('$\Delta$DEC [arcsec]')
#     plt.title('Source %s, p_any=%.2f' % (args.id, p_any))
#     xlo, xhi = plt.xlim()
#     ylo, yhi = plt.ylim()
#     hi = max(-xlo, xhi, -ylo, yhi)
#     plt.ylim(-hi, hi) # DEC
#     plt.xlim(hi, -hi) # RA goes the other way
#     plt.legend(loc='best', numpoints=1, prop=dict(size=8))
#     outfilename = '%s_explain_%s.pdf' % (args.matchcatalogue, args.id)
#     print('plotting to %s' % outfilename)
#     print()
#     print("Disclaimer: These results assume that the input (sky densities, positional errors, and priors) are correct.")
#     print()
#     plt.savefig(outfilename, bbox_inches='tight')
#     plt.close()
