__author__ = 'maximkuleshov'

import matplotlib.pyplot as plt
import seaborn as sns
import pickle


def draw_hist_cmp(pv, opv, apv, aopv):
    print(sns.axes_style())
    sns.set(style="ticks", color_codes=True)
    sns.set(style="ticks", palette="muted", color_codes=True)
    f, axes = plt.subplots(figsize=(7, 7), nrows=2, ncols=2, sharex='col', sharey='row')

    # sns.distplot(pv, ax=axes[0, 0], color='b')
    # axes[0, 0].set_title('p-value')
    # sns.distplot(apv, ax=axes[0, 1], color='b')
    # axes[0, 1].set_title('adjusted p-value')
    # sns.distplot(opv, ax=axes[1, 0], color='g')
    # axes[1, 0].set_title('old p-value')
    # sns.distplot(aopv, ax=axes[1, 1], color='g')
    # axes[1, 1].set_title('old adjusted p-value')

    plt.tight_layout()

    # plt.hist(result1, alpha=0.5, color='blue', linewidth=0.8, label=label1)
    # plt.hist(result2, alpha=0.5, color='green', linewidth=0.8, label=label2)
    #
    # yint = range(min(pv + apv), max(pv + apv) + 1, 5)
    # axes[0, 0].set_yticks(yint)
    # plt.yticks(yint)
    # plt.title(title.replace('_', ' '))
    # plt.xlabel('ranks')
    # plt.ylabel('matches')
    # plt.legend(loc='upper right')

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