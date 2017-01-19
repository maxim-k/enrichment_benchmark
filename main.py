__author__ = 'maximkuleshov'

from collections import defaultdict
import json
import requests
from time import sleep
from operator import itemgetter
# import matplotlib.pyplot as plt


ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/'
REF_SIZE = 104


def get_enrichr_results(gene_set_library, genelist, description):
    addlist_url = ENRICHR_URL + 'addList'
    payload = {
        'list': (None, genelist),
        'description': (None, description)
    }

    response = requests.post(addlist_url, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    sleep(1)
    data = json.loads(response.text)

    enrich_url = ENRICHR_URL + '/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    response = requests.get(enrich_url + query_string % (user_list_id, gene_set_library))
    sleep(1)
    return json.loads(response.text)


def parse_gmt(gmt):
    tfs = defaultdict(list)
    for line in gmt:
        term, desc, *genes = line.strip().split('\t')
        genes = [gene.split(',')[0] for gene in genes]
        tf = term.split(sep='_')[0]
        tfs[tf].append(genes)
    return tfs


def map_tf(tf, res, ref):
    indices = [i for i, x in enumerate(res) if x == tf]
    for index in indices:
        ref[index] += 1
    return ref


def main():
    library = 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X'
    gmt_file = 'ChEA_2016.gmt'
    chea2016 = parse_gmt(open(gmt_file, 'r').readlines())

    pval_hist = [0] * REF_SIZE
    adj_pval_hist = [0] * REF_SIZE
    old_pval_hist = [0] * REF_SIZE
    old_adj_pval_hist = [0] * REF_SIZE

    for key in chea2016.keys():
        for genes in chea2016[key]:
            data = get_enrichr_results(library, '\n'.join(genes), '')
            results = []
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
    return None


if __name__ == '__main__':
    main()
