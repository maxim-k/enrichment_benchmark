__author__ = 'maximkuleshov'

from collections import defaultdict
import json
import requests
from time import sleep
from operator import itemgetter

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


def parse_gmt(gmt):
    tfs = defaultdict(list)
    for line in gmt:
        term, desc, *genes = line.strip().split('\t')
        genes = [gene.split(',')[0] for gene in genes]
        tf = term.split(sep='_')[0]
        tfs[tf].append(genes)
    return tfs


def main():
    library = 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X'
    gmt_file = 'ChEA_2016.gmt'
    chea2016 = parse_gmt(open(gmt_file, 'r').readlines())
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
            s_pval = sorted(results, key=itemgetter(1))
            s_adj_pval = sorted(results, key=itemgetter(2))
            s_old_pval = sorted(results, key=itemgetter(3))
            s_old_adj_pval = sorted(results, key=itemgetter(4))
            print()
    return None


if __name__ == '__main__':
    main()
