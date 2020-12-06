# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import pandas as pd
import sys
import numpy as np

sys.path.append("/usr/local/lib/python3.8/site-packages")
import RNA


def fold_rna(seq):
    # # Use a breakpoint in the code line below to debug your script.
    # p = Popen(path,stdin=PIPE,stdout=PIPE)
    # ans = p.communicate(seq.encode())
    # ans = str (ans[0].decode()).replace('\r','').split('\n')[1]
    # return ans[: ans.index(' ')]
    return RNA.fold(seq)


def draw_2nd_structure(srna):
    for i in range(3):
        name = srna.loc[i, "name"]
        sequence = srna.loc[i, "sequence"]
        folding = fold_rna(srna.loc[i, "sequence"])[0]
        print("##############", name, "###############")
        print(sequence)
        print(folding)
        print()

        fx_test = ">{0}\n{1}\n{2}\n".format(name, sequence, folding)
        textfile = open('./resources/' + srna.loc[i, "name"] + '.fx', "w")
        textfile.write(fx_test)
        textfile.close()
    for i in range(3):
        ##Print Structure
        plt.figure(figsize=(20, 20))
        cg = forgi.load_rna('./resources/' + srna.loc[i, "name"] + '.fx', allow_many=False)
        fvm.plot_rna(cg, text_kwargs={"fontweight": "black"}, lighten=0.7, backbone_kwargs={"linewidth": 3})
        plt.show()
        # plt.savefig(srna.loc[i, "name"]+'.png')
    return srna


def find_str(s, char):
    index = 0

    if char in s:
        c = char[0]
        for ch in s:
            if ch == c:
                if s[index:index + len(char)] == char:
                    return index

            index += 1

    return -1


def align(sequence, fragments):
    align = list()
    align.append(list())
    align[0].append(sequence)
    for frag in fragments:
        index = find_str(sequence, frag)
        spaces = ' ' * (index)
        frag = spaces + frag
        frag = frag + ' ' * (len(sequence) - len(frag))
        align[0].append(frag)
    return align


def count_Nuc(alignment):
    count_array = []
    for i in range(len(alignment[0][0])):
        counter = 0
        for j in range(1, len(alignment[0])):
            if alignment[0][j][i] == alignment[0][0][i]:
                counter += 1
        count_array.append(counter)
    return count_array


def count_hybrid(alignment, sRNA,mRNA, df):
    count_array = []
    for i in range(len(alignment[0][0])):
        counter = 0
        for j in range(1, len(alignment[0])):
            if alignment[0][j][i] == alignment[0][0][i]:
                counter += df.loc[df['sRNA'] == sRNA][df["sRNA fragment (5'-3')"] == alignment[0][j].strip()]
        count_array.append(counter)
    return count_array


def split(word):
    return [char for char in word]


def histogram(alignment, sRNA, unique, df=None):
    y_axis = []
    if unique:
        y_axis = count_Nuc(alignment)
    else:
        y_axis = count_hybrid(alignment, sRNA, df)
    x_axis = split(alignment[0][0])
    # style
    plt.style.use('seaborn-darkgrid')
    # create a color palette
    palette = plt.get_cmap('Set1')
    x = range(len(x_axis))
    plt.xticks(x, x_axis)
    # multiple line plot
    plt.plot(x, y_axis, marker='', color=palette(1), linewidth=1)
    plt.tick_params(axis='x', which='major', labelsize=5)
    plt.ylim(ymax=max(y_axis) + 1, ymin=min(y_axis) - 1)
    ymax = max(y_axis)
    xpos = y_axis.index(ymax)
    xpos2 = len(y_axis) - 1 - y_axis[::-1].index(ymax)
    xmax = x[xpos]
    xmax2 = x[xpos2]
    plt.annotate('most likeable interaction fragment', xy=(xmax, ymax), xytext=(xmax, ymax + 0.5),
                 arrowprops=dict(facecolor='black', shrink=0.001))
    plt.annotate(text="", xy=(xmax2, ymax), xytext=(xmax2, ymax + 0.5),
                 arrowprops=dict(facecolor='black', shrink=0.001))
    # plt.bar(x, y_axis, color=palette(1), linewidth=1, alpha=0.9)
    # Add legend
    # plt.legend(loc=2, ncol=2)
    # Add titles
    plt.title("interaction", loc='left', fontsize=12, fontweight=0, color='orange')
    plt.xlabel("Nucleotide")
    plt.ylabel("Count")
    # plt.show()
    plt.savefig(sRNA + '.png')
    plt.clf()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    rdlB = "rdlB"
    rybB = "rybB"
    ryhB = "ryhB"
    srna = pd.read_csv('./resources/srna.csv')
    rdlBseq = srna.loc[0, "sequence"]
    rybBseq = srna.loc[1, "sequence"]
    ryhBseq = srna.loc[2, "sequence"]

    rdlBfrag = srna.loc[srna['sRNA'] == rdlB][["sRNA fragment (5'-3')","mRNA"]]
    rybBfrag = srna.loc[srna['sRNA'] == rybB][["sRNA fragment (5'-3')","mRNA"]]
    ryhBfrag = srna.loc[srna['sRNA'] == ryhB][["sRNA fragment (5'-3')","mRNA"]]

    print("#############################    rdlB    #############################\n")
    rdlB_sys = align(rdlBseq, rdlBfrag["sRNA fragment (5'-3')"])

    print("\n#############################    rybB    #############################\n")
    rybB_sys = align(rybBseq, rybBfrag["sRNA fragment (5'-3')"])

    print("\n#############################    ryhB    #############################\n")
    ryhB_sys = align(ryhBseq, ryhBfrag["sRNA fragment (5'-3')"])

    # plot
    histogram(rdlB_sys, rdlB, True)
    # plot
    histogram(rybB_sys, rybB, True)
    # plot
    histogram(ryhB_sys, ryhB, True)
    # plot hybrid
    histogram(rdlB_sys, rdlB, False, srna)
    # plot hybrid
    histogram(rybB_sys, rybB, False, srna)
    # plot hybrid
    histogram(ryhB_sys, ryhB, False, srna)
