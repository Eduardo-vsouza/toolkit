import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



def count(data):
    df = pd.read_csv(data, sep='\t')
    fixed = []
    peps = set(df["peptide"].tolist())
    for pep in peps:
        fpep = pep.replace("-", "").replace('.', '').replace("[UNIMOD:", '').replace("]", '')
        fpep = ''.join([i for i in fpep if not i.isdigit()])
        if fpep not in fixed:
            fixed.append(fpep)

    for pep in fixed:
        print(pep)

def quick_plot():
    # fig = plt.figure()
    data = {'seps': [49, 3, 1838, 33], 'origin': ['Human', 'Virus', 'Human', 'Virus'], 'category': ['Novel', 'Novel',
                                                                                                    'Annotated', 'Annotated']}
    # ax = fig.add_axes([0, 0, 1, 1])
    df = pd.DataFrame(data)

    students = [49, 3, 1838, 33]
    # plt.bar(langs, students, edgecolor='black')
    sns.barplot(x="category", y="seps", hue="origin", data=df)

    plt.show()

# c8 = '/media/eduardo/DATA/Eduardo/adenovirus/proteogenomics/nanopore_assemblies/human_virus/01_03_22/102_c8_hr/Transcriptome/post_perc/results/'

# total = '/media/eduardo/DATA/Eduardo/adenovirus/proteogenomics/nanopore_assemblies/human_virus/01_03_22/102_total_hr/Transcriptome/post_perc/results/'
# count(data='/media/eduardo/DATA/Eduardo/adenovirus/proteogenomics/nanopore_assemblies/human_virus/01_03_22/102_c8_hr/Transcriptome/post_perc/results/human_novel.xls')
# count('/media/eduardo/DATA/Eduardo/adenovirus/proteogenomics/nanopore_assemblies/human_virus/01_03_22/102_total_hr/Transcriptome/post_perc/results/human_novel.xls')
# count('/media/eduardo/DATA/Eduardo/adenovirus/proteogenomics/nanopore_assemblies/human_virus/01_03_22/102_c8_hr/Transcriptome/post_perc/results/virus_novel.xls')
# count('/media/eduardo/DATA/Eduardo/adenovirus/proteogenomics/nanopore_assemblies/human_virus/01_03_22/102_total_hr/Transcriptome/post_perc/results/virus_novel.xls')
# count(f'{c8}/human_all.txt')
# # count(f'{total}/human_all.txt')
# # count(f'{c8}/virus_all.txt')
# count(f'{total}/virus_all.txt')
quick_plot()