__author__ = 'maximkuleshov'

from collections import defaultdict
import json
import requests
from time import sleep
from operator import itemgetter
import matplotlib.pyplot as plt


ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/'

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


def parse_gmt(gmt, dir):
    tfs = defaultdict(list)
    for line in gmt:
        term, desc, *genes = line.strip().split('\t')
        genes = [gene.split(',')[0] for gene in genes]
        if term.split(sep='-')[1] == dir:
            tf = term.split(sep='-')[0].upper()
            tfs[tf].append(genes)
    return tfs, len(gmt)


def map_tf(tf, res, ref):
    indices = [i for i, x in enumerate(res) if x == tf]
    for index in indices:
        ref[index] += 1
    return ref


def draw_hist_cmp(result1, result2, label1, label2, title, ref_size):
    plt.xlim(-5, 20)
    # plt.ylim(0, 200)
    plt.hist(result1, alpha=0.5, color='blue', bins=ref_size, label=label1)
    plt.hist(result2, alpha=0.5, color='green', bins=ref_size, label=label2)
    plt.xlabel('tf hits in library')
    plt.title(title)
    plt.xlabel('ranks')
    plt.ylabel('matches')
    plt.legend(loc='upper right')
    plt.show()
    return None


def main():
    libraries = ['ChEA_2016'] #, 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X']
    dirs = ['up'] #, 'down']
    gmt_file = 'single_gene_perturbations-v1.0.gmt'

    for direction in dirs:
        for library in libraries:
            reference, ref_size = parse_gmt(open(gmt_file, 'r').readlines(), direction)
            pval_hist = [0] * ref_size
            adj_pval_hist = [0] * ref_size
            old_pval_hist = [0] * ref_size
            old_adj_pval_hist = [0] * ref_size

            for key in list(reference.keys())[:100]:
                for genes in reference[key]:
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
            draw_hist_cmp(pval_hist, old_pval_hist, 'p-value', 'old p-value', '%s %s' % (library, direction), ref_size)
            draw_hist_cmp(adj_pval_hist, old_adj_pval_hist, 'adjusted p-value', 'old adjusted p-value', '%s %s' % (library, direction), ref_size)
    return None

if __name__ == '__main__':
    main()
