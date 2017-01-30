__author__ = 'maximkuleshov'

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pickle


def binify(res, size):
    bins_chain = []
    bins = {'label': [], 'sum': []}
    for i in range(0, len(res) - 1, size):
        bins_chain.append([i, i + size])
    bins_chain[-1][1] = len(res)

    for bin in bins_chain:
        sum = 0
        for i in range(bin[0], bin[1]):
            sum += res[i]
        bins['label'].append(bin[0])
        bins['sum'].append(sum)
    return bins


def histify(res):
    hist = []
    for pos, i in enumerate(res):
        hist.extend([pos]*i)
    return hist


def draw_hist_cmp(pv, opv, apv, aopv):
    bins = binify(pv, 10)
    ind = range(len(bins))  # the x locations for the groups
    width = 5  # the width of the bars
    # plt.xticks(bins + width / 2)
    plt.bar(bins['label'], bins['sum'], width)

    sns.set(style="white", palette="muted", color_codes=True)
    f, axes = plt.subplots(2, 2, figsize=(7, 7))

    # Plot a simple histogram with binsize determined automatically
    sns.distplot(histify(pv), color="b", ax=axes[0, 0], bins=65, axlabel='p-value')
    plt.setp(axes, yticks=[])
    plt.hist(histify(pv), bins=65)
    plt.tight_layout()
    sns.despine(left=True)
    plt.show()
    return None


def main():
    chea_up = pickle.load(open('ChEA_2016_up.pickle', 'rb'))
    chea_dn = pickle.load(open('ChEA_2016_dn.pickle', 'rb'))
    encode_up = pickle.load(open('ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_up.pickle', 'rb'))
    encode_dn = pickle.load(open('ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_dn.pickle', 'rb'))
    pval_hist, old_pval_hist, adj_pval_hist, old_adj_pval_hist = chea_up[1:]
    draw_hist_cmp(pval_hist, old_pval_hist, adj_pval_hist, old_adj_pval_hist)
    return None

if __name__ == '__main__':
    main()