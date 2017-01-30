__author__ = 'maximkuleshov'

from collections import defaultdict
import json
import requests
import os.path
import pickle
from retrying import retry
from operator import itemgetter



ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/'


@retry
def get_enrichr_results(gene_set_library, genelist, description):
    addlist_url = ENRICHR_URL + 'addList'
    payload = {
        'list': (None, genelist),
        'description': (None, description)
    }

    response = requests.post(addlist_url, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    # sleep(1)
    data = json.loads(response.text)

    enrich_url = ENRICHR_URL + '/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    response = requests.get(enrich_url + query_string % (user_list_id, gene_set_library))
    # sleep(1)
    return json.loads(response.text)


def parse_gmt(gmt, direction):
    tfs = defaultdict(list)
    for line in gmt:
        term, desc, *genes = line.strip().split('\t')
        genes = [gene.split(',')[0] for gene in genes]
        if term.split(sep='-')[1] == direction:
            tf = term.split(sep='-')[0].upper()
            tfs[tf].append(genes)
    return tfs


def filter_library(ref, lib):
    lib_keys = set(line.split()[0].split(sep='_')[0] for line in lib)
    filtered_keys = set(ref.keys()).intersection(lib_keys)
    filtered_ref = []
    for key in sorted(filtered_keys):
        filtered_ref.extend([[key, ref[key][pos]] for pos in range(len(ref[key]))])
    return sorted(filtered_ref)


def map_tf(tf, res, ref):
    indices = [i for i, x in enumerate(res) if x == tf]
    for index in indices:
        ref[index] += 1
    return ref


def main():
    libraries = [['ChEA_2016', 645], ['ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 104]]
    dirs = ['up', 'dn']
    gmt_file = 'single_gene_perturbations-v1.0.gmt'

    for direction in dirs:
        for lib_data in libraries:
            library, lib_size = lib_data
            reference = parse_gmt(open(gmt_file, 'r').readlines(), direction)
            lib_file = open('%s.gmt' % library, 'r').readlines()
            reference = filter_library(reference, lib_file)

            # if os.path.isfile('%s_%s.pickle' % (library, direction)):
            #     jar = pickle.load(open('%s_%s.pickle' % (library, direction), 'rb'))
            #     start_pos, pval_hist, adj_pval_hist, old_pval_hist, old_adj_pval_hist = jar
            # else:
            start_pos = 0
            pval_hist = [0] * lib_size
            adj_pval_hist = [0] * lib_size
            old_pval_hist = [0] * lib_size
            old_adj_pval_hist = [0] * lib_size

            for pos, line in enumerate(reference[start_pos:]):
                key, genes = line
                if not genes:
                    continue
                results = []
                data = get_enrichr_results(library, '\n'.join(genes), '')
                for res in data[library]:
                    tf = res[1].split(sep='_')[0]
                    pval = res[2]
                    adj_pval = res[6]
                    old_pval = res[7]
                    old_adj_pval = res[8]
                    results.append([tf, pval, adj_pval, old_pval, old_adj_pval])
                s_pval = [line[0] for line in sorted(results, key=itemgetter(1))]
                pval_hist = map_tf(key, s_pval, pval_hist)

                s_adj_pval = [line[0] for line in sorted(results, key=itemgetter(2))]
                adj_pval_hist = map_tf(key, s_adj_pval, adj_pval_hist)

                s_old_pval = [line[0] for line in sorted(results, key=itemgetter(3))]
                old_pval_hist = map_tf(key, s_old_pval, old_pval_hist)

                s_old_adj_pval = [line[0] for line in sorted(results, key=itemgetter(4))]
                old_adj_pval_hist = map_tf(key, s_old_adj_pval, old_adj_pval_hist)

                status = [start_pos + pos + 1, pval_hist, adj_pval_hist, old_pval_hist, old_adj_pval_hist]
                pickle.dump(status, open('%s_%s.pickle' % (library, direction), 'wb'))
    return None

if __name__ == '__main__':
    main()
