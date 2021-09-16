import os

def writefile(protid, out_dir, pident, isoform, consequence, df, maptype, csv= False, hdf = False):
    if csv is True: 
        out_csv = os.path.join(out_dir, 'csv')
        with open(os.path.join(out_csv,(maptype + '_pident' + str(pident) + '_isoform_' +
                '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.csv')), 'a') as f:
           df.to_csv(f, sep=',', index=False,
                                          header=f.tell() == 0)
    if hdf is True: 
        out_hdf = os.path.join(out_dir, 'hdf5', maptype)
        if not os.path.exists(out_hdf):
            os.mkdir(out_hdf)
        import vaex
        fn = os.path.join(out_hdf,( maptype + '_pident' + str(pident) + '_isoform_' +
                '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '_' + protid + '.hdf5'))
        vaex_df = vaex.from_pandas(df, copy_index=False)
        vaex_df.export(fn)