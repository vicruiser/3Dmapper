# coding: utf-8
import os
import pandas as pd
from tabulate import tabulate
from .run_subprocess import call_subprocess


def stats(var_statsfile, int_statsfile, mapped_infofile, out_dir):
    # # PDBmapper stats
    # var_stats = pd.read_csv(var_statsfile, sep=" |\t", engine='python')
    # int_stats = pd.read_csv(int_statsfile, sep=" |\t", engine='python')

    # # number of mapped variants
    # n_variants_unique_cmd = "awk 'NR>1{{print $1}}' {} | uniq | wc -l"
    # n_variants_unique, err1 = call_subprocess(
    #     n_variants_unique_cmd.format(mapped_infofile))

    # n_variants_cmd = "awk 'NR>1{{print $1}}' {}  | wc -l"
    # n_variants, err1 = call_subprocess(
    #     n_variants_cmd.format(mapped_infofile))

    # # number of mapped genes
    # n_genes_cmd = "awk 'NR>1{{print $4}}' {} | uniq | wc -l"
    # n_genes, err2 = call_subprocess(n_genes_cmd.format(mapped_infofile))

    # # process results
    # n_variants, n_variants_unique, n_genes = n_variants.decode(
    #     'utf-8').rstrip(), n_variants_unique.decode(
    #     'utf-8').rstrip(), n_genes.decode('utf-8').rstrip()

    # # number of mapped interfaces, type of interfaces and consequences
    # interfaces_info = pd.read_csv(mapped_infofile, sep=' ', usecols=[
    #                               'interface_id', 'Consequence'])
    # interfaces_info[['pdb_id', 'ENSP', 'temp_chain', 'int_chain', 'type']
    #                 ] = interfaces_info.interface_id.str.split("_", expand=True)

    # # Consequence table counts
    # summary = interfaces_info.groupby(
    #     ['type', 'Consequence']).size().reset_index(name='Count')
    # summary['type'] = summary['type'].replace({'ligand': 'interface_with_ligand',
    #                                            'protein': 'interface_with_protein',
    #                                            'nucleic': 'interface_with_dna'})

    # summary = summary.pivot(index='Consequence',
    #                         columns='type', values='Count')
    # summary = summary.fillna(0)

    # summary.to_csv(os.path.join(out_dir, 'pdbmapper_consequences_stats.info'),
    #                sep=' ', encoding='utf-8', index=True)
    # stats_table = tabulate(
    #     summary, headers='keys', tablefmt='psql')

    # # number of mapped proteins
    # n_prot = interfaces_info['ENSP'].drop_duplicates().count()

    # # number of mapped interfaces
    # n_int = interfaces_info['interface_id'].drop_duplicates().count()

    # # number of interfaces ligand
    # n_int_ligand = interfaces_info.loc[
    #     interfaces_info['type'] == "ligand", 'interface_id'].drop_duplicates().count()

    # # number of interfaces protein
    # n_int_prot = interfaces_info.loc[
    #     interfaces_info['type'] == "protein", 'interface_id'].drop_duplicates().count()

    # # number of interfaces nucleic
    # n_int_dna = interfaces_info.loc[
    #     interfaces_info['type'] == "nucleic", 'interface_id'].drop_duplicates().count()

    # # stats table
    # stats_message = {'Variable': ['Variants ids', 'Unique variants ids',
    #                               'Interface ids', 'Interface ids (dna)',
    #                               'Interface ids (ligand)', 'Interface ids (protein)',
    #                               'Protein ids', 'Gene ids'],
    #                  'Total':
    #                  [int(n_variants), int(n_variants_unique),
    #                   n_int, n_int_dna,
    #                   n_int_ligand, n_int_prot,
    #                   n_prot, int(n_genes)],
    #                  '% (mapped / input)':
    #                  ['-',
    #                   round((int(n_variants_unique) /
    #                          var_stats.iloc[0]['n_variants']) * 100, 2),
    #                   round((int(n_int) /
    #                          int_stats.iloc[0]['n_interfaces']) * 100, 2),
    #                   round((int(n_int_dna) /
    #                          int_stats.iloc[0]['n_interfaces_nucleic']) * 100, 2),
    #                   round((int(n_int_ligand) /
    #                          int_stats.iloc[0]['n_interfaces_ligand']) * 100, 2),
    #                   round((int(n_int_prot) /
    #                          int_stats.iloc[0]['n_interfaces_protein']) * 100, 2),
    #                   round((int(n_prot) /
    #                          int_stats.iloc[0]['n_ENSP']) * 100, 2),
    #                   round((int(n_genes) /
    #                          var_stats.iloc[0]['n_genes']) * 100, 2)]}
    # # conver it to data frame and save it
    # stats_message = pd.DataFrame(stats_message)
    # stats_message.set_index('Variable')
    # stats_message2 = tabulate(
    #     stats_message, headers='keys', tablefmt='psql')
    # stats_message.to_csv(os.path.join(out_dir, 'pdbmapper_stats.info'),
    #                      sep='\t', encoding='utf-8', index=False)

    # return stats_message2, stats_table
    pass
