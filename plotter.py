import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3



def get_csv(path):
    print(f'Reading Files: {path}')
    df = pd.read_csv(path, header=0)
    return df


def filter_csv(csv):
    print(f'Filtering: {csv}')
    filtered_csv = csv.loc[(csv['pvalue']<=0.05) & (csv['adj_pvalue']<= 0.05)].copy()
    filtered_csv = filtered_csv.dropna(subset=['feature_id'])
    return filtered_csv


def venn_dia(all_data, labels=('Aza vs Control', 'Gua vs Control', 'Aza vs Gua')):
    print('Plotting')
    sets = [
        set(data['feature_id'].unique())
        for data in all_data
    ]
    plt.figure(figsize=(9, 7))
    venn = venn3(
        subsets=sets,
        set_labels=labels
    )
    plt.title('Significant Differential Transcript Usage (DRIMSeq)')
    plt.tight_layout()
    plt.savefig(
        '../results/EPIC_CAM_SLR24/drimseq/venn_transcripts.png',
        dpi=300,
        bbox_inches='tight'
    )
    plt.close()

    all_sign = sets[0] & sets[1] & sets[2]
    aza_gua = sets[0] & sets[1]
    aza_coef = sets[0] & sets[2]
    gua_coef = sets[1] & sets[2]
    stats = {
    'category': [
        'Significant transcripts',
        'Significant transcripts',
        'Significant transcripts',
        'Shared transcripts',
        'Shared transcripts',
        'Shared transcripts',
        'Shared transcripts'
    ],
    'comparison': [
        labels[0],
        labels[1],
        labels[2],
        'All three',
        f'{labels[0]} & {labels[1]}',
        f'{labels[0]} & {labels[2]}',
        f'{labels[1]} & {labels[2]}'
    ],
    'n_transcripts': [
        len(sets[0]),
        len(sets[1]),
        len(sets[2]),
        len(all_sign),
        len(aza_gua),
        len(aza_coef),
        len(gua_coef)
    ]
    }

    csv_res = {
        'all': all_sign,
        'aza_control_gua_control': aza_gua,
        'aza_control_aza_gua': aza_coef,
        'gua_control_aza_gua': gua_coef
    }
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(
        '../results/EPIC_CAM_SLR24/drimseq/venn_statistics.csv',
        index=False
    )
    for nam, transcript_set in csv_res.items():
        df = pd.DataFrame(sorted(transcript_set), columns=['feature_id'])
        df.to_csv(f'../SLR24_results/drimseq/{nam}_sign.csv',index=False)
    return stats_df


if __name__=='__main__':
    f_path = ['../results/EPIC_CAM_SLR24/drimseq/res_aza_txp.csv',
              '../results/EPIC_CAM_SLR24/drimseq/res_gua_txp.csv',
              '../results/EPIC_CAM_SLR24/drimseq/res_coef_txp.csv']
    all_data_list = []
    for path in f_path:
        data = filter_csv(get_csv(path))
        all_data_list.append(data)
        name = os.path.splitext(path)[0] + '_filtered.csv'
        data.to_csv(name, index=False)
    venn_dia(all_data_list)
