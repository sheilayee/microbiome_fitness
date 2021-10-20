import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd
from scipy.stats import sem
import matplotlib.patches as patches
from matplotlib.legend_handler import HandlerPatch
from matplotlib.offsetbox import OffsetImage,AnnotationBbox

def relative_abundance(rows):
    # find the mean number of CFUs for each bacteria for each of the 32 treatments
    LP_mean = mean_CFU_per_treatment(rows, 6)
    LB_mean = mean_CFU_per_treatment(rows, 7)
    AP_mean = mean_CFU_per_treatment(rows, 8)
    AT_mean = mean_CFU_per_treatment(rows, 9)
    AO_mean = mean_CFU_per_treatment(rows, 10)

    # find the mean number of total CFUs for each of the 32 treatments
    mean_CFUs = mean_CFU_per_treatment(rows, 11)

    # calculate the proportion of each bacteria for each of the 32 treatments
    # to calculate this, use the mean CFU count of each bacteria to determine the percentage of its presence out of the
    # mean of the total CFUS. Do this for every treatment
    # LP_ratios = ratios(LP_mean, mean_CFUs)
    # # print(LP_ratios)
    # LB_ratios = ratios(LB_mean, mean_CFUs)
    # # print(LB_ratios)
    # AP_ratios = ratios(AP_mean, mean_CFUs)
    # # print(AP_ratios)
    # AT_ratios = ratios(AT_mean, mean_CFUs)
    # # print(AT_ratios)
    # AO_ratios = ratios(AO_mean, mean_CFUs)
    # # print(AO_ratios)

    # find the SEM
    LP_SEM = SEM(rows, 6)
    LB_SEM = SEM(rows, 7)
    AP_SEM = SEM(rows, 8)
    AT_SEM = SEM(rows, 9)
    AO_SEM = SEM(rows, 10)
    total_CFUs_SEM = SEM(rows, 11)

    # concatenate all of the ratio arrays together into one array
    all_sums = np.zeros((5, 32))
    # append a row of ratios to ith row in all_ratios array
    all_sums[0, :] = LP_mean
    all_sums[1, :] = LB_mean
    all_sums[2, :] = AP_mean
    all_sums[3, :] = AT_mean
    all_sums[4, :] = AO_mean

    # transpose the all_ratios array
    all_sums = all_sums.T
    # print(all_sums.shape)

    d = np.zeros(all_sums.shape)
    for j, row in enumerate(all_sums):
        g = np.zeros(len(row) + 1)
        G = np.sum(row)
        g[1:] = np.cumsum(row)
        f = 10 ** (g / G * np.log10(G))
        f[0] = 0
        d[j, :] = np.diff(f)
    # print('Preserving linear ratios logarithmically')
    # print(d)
    # print(d.shape)

    df = pd.DataFrame(d)
    # print(df)

    ax = df.plot.bar(stacked=True, log=True, figsize=(12, 3), color=['orangered', 'gold', 'darkorchid', 'royalblue', 'limegreen'],
                 edgecolor='black', linewidth=.6, width=.8, legend=None)
    # df2.plot.bar(stacked=True, log=True, figsize=(10, 3), color=['orangered', 'gold', 'darkorchid', 'royalblue', 'limegreen'],
    #              edgecolor='black', linewidth=.6, width=.8, yerr=total_CFUs_SEM)
    plt.ylabel('Bacterial load(CFUs)')
    plt.ylim(1, 10**6)
    plt.yticks([10**1, 10**2, 10**3, 10**4, 10**5, 10**6])

    colors = ["r", "y", "p", "b", "g",
              "ry", "rp", "rb", "rg", "yp", "yb", "yg", "pb", "pg", "bg",
              "pry", "bry", "gry", "prb", "prg", "brg", "pyb", "pyg", "byg", "pgb",
              "pryb", "pryg", "bryg", "prgb", "pygb",
              "prygb",
              "w"]
    for index, col in enumerate(colors):
        # for x in range(0, 5):
        pie_AB = offset_pie_image(index, col, ax)
        # print(ab)
        ax.add_artist(pie_AB)

    # add error bars
    x = []
    for i in range(0, 32):
        x.append(i)
    ax.errorbar(x, mean_CFUs, yerr=total_CFUs_SEM, capsize=2, lw=.8, ms=8, zorder=3, ls='none', color='black')

    # make x axis pie charts
    plt.tick_params(axis='y', top=True, bottom=True, right='off', which='both', direction='in', labeltop=True)
    plt.subplots_adjust(top=.964, left=.086, bottom=.1, right=.83, hspace=.28, wspace=.268)

    # plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    # ax.tick_params(axis='x', which='minor', direction='in')

    # make legend denoting the bacterial species that is associated with each color
    circ1 = patches.Circle((.1, .1), 1, facecolor='orangered', edgecolor='black', linewidth=.3)
    circ2 = patches.Circle((0, 0), 1, facecolor='gold', edgecolor='black', linewidth=.3)
    circ3 = patches.Circle((0, 0), 1, facecolor='darkorchid', edgecolor='black', linewidth=.3)
    circ4 = patches.Circle((0, 0), 1, facecolor='royalblue', edgecolor='black', linewidth=.3)
    circ5 = patches.Circle((0, 0), 1, facecolor='limegreen', edgecolor='black', linewidth=.3)
    circ6 = patches.Circle((0, 0), 1, facecolor='whitesmoke', edgecolor='black', linewidth=.8)
    plt.legend((circ1, circ2, circ3, circ4, circ5, circ6),
               ('$\it{Lactobacillus\,plantarum}$',
                '$\it{Lactobacillus\,brevis}$',
                '$\it{Acetobacter\,pasteurianus}$',
                '$\it{Lactobacillus\,tropicalis}$',
                '$\it{Lactobacillus\,orientalis}$',
                'germ-free'),
               handler_map={patches.Circle: Circle(), },
               fontsize=7,
               loc='upper left',
               bbox_to_anchor=(1.01, .8),
               columnspacing=0.5,
               labelspacing=1,
               handletextpad=1,
               edgecolor='black')
    ax.set_xticklabels([])


    plt.show()

def get_pie_png(name):
    path_to_png = "/Users/Sheila/Documents/BF550/Homework/Project1/pies/{}.png".format(name)
    im = plt.imread(path_to_png)
    return im

def offset_pie_image(coordinate, name, ax):
    # print(coord)
    image = get_pie_png(name)
    im = OffsetImage(image, zoom=.01)
    im.image.axes = ax

    # note that the y axis is on a logarithmic scale. For the xy parameter, a y value of 1 is actually
    # the log of 1 which is 0...this is below the x-axis of the plot, hence why the png images of the 32
    # pie charts are plotted below the x-axis
    pie_AB = AnnotationBbox(im, xy=(coordinate, 1), xybox=(0, -10), frameon=False,
                        xycoords='data',  boxcoords="offset points", pad=0)
    return pie_AB
    # ax.annotate(ab)

def mean_CFU_per_treatment(rows, column):
    means = [] # create empty list to store means
    num_treatments = 32  # the total number of treatments
    start_flies_index = 13
    end_flies_index = 36

    # for replicates 2 and 3
    for treatment in range(num_treatments): # iterate through 32 treatments
        sum = 0
        for i in range(start_flies_index, end_flies_index+1):
            sum += float(rows[i][column])
            # sum += math.log10(float(rows[i][column])) # add CFU count of replicate 1 to the sum
            # sum += math.log10(float(rows[i+36][column])) # add CFU count of replicate 4 to the sum
            mean = sum/24 # find the mean

        means.append(mean)
        start_flies_index += 48
        end_flies_index += 48

    # for replicates 1 and 4
    # for treatment in range(num_treatments): # iterate through 32 treatments
    #     sum = 0
    #     for i in range(start_flies_index, end_flies_index+1):
    #         # print(str(i) + ': the sum is ' + str(sum))
    #         # print('...' + rows[i][column])
    #         print(i)
    #         # print(rows[i][column])
    #         # print(rows[i+36][column])
    #         print(means)
    #         sum += float(rows[i][column])
    #         sum += float(rows[i+36][column])
    #         # sum += math.log10(float(rows[i][column])) # add CFU count of replicate 1 to the sum
    #         # sum += math.log10(float(rows[i+36][column])) # add CFU count of replicate 4 to the sum
    #         mean = sum/24 # find the mean
    #     means.append(mean)
    #     start_flies_index += 48
    #     end_flies_index += 48

    return means

def SEM(rows, column):
    # calculate SEM for mean total CFU for each bacterial combination

    all_data = []  # create empty list to store means
    num_treatments = 32  # the total number of treatments
    start_flies_index = 13
    end_flies_index = 36

    # for replicates 2 and 3
    for treatment in range(num_treatments):  # iterate through 32 treatments
        data = []
        sum = 0
        for i in range(start_flies_index, end_flies_index + 1):
            x = float(rows[i][column])
            data.append(x)
        all_data.append(data)
        start_flies_index += 48
        end_flies_index += 48

    # print(all_data)
    # # x = all_data[0][0]

    SEMs = []
    for i in range(num_treatments):  # iterate through 32 treatments
        x = sem(all_data[i])
        SEMs.append(x)
    # print(SEMs)

    return SEMs

def ratios(bacterial_mean, CFU_mean):
    ratios = []
    for i in range(len(bacterial_mean)):
        # print(i)
        if CFU_mean[i] > 0:
            ratio = bacterial_mean[i]/CFU_mean[i]
            # print(ratio)
            ratios.append(ratio)
        else:
            ratios.append(0)
    # print(ratios)
    ratios_array = np.array(ratios)
    print(ratios_array)
    # ratios_array_resized = ratios_array[np.newaxis, :]
    # print(ratios_array_resized)
    return ratios_array

def pie_colors():
    # 100% pie graphs
    r = ['orangered']
    y = ['gold']
    p = ['darkorchid']
    b = ['royalblue']
    g = ['limegreen']

    # make 50%/50% pie graphs (2 treatment combinations)
    ry = ['orangered','gold']
    rp = ['orangered', 'darkorchid']
    rb = ['orangered', 'royalblue']
    rg = ['orangered', 'limegreen']
    yp = ['gold', 'darkorchid']
    yb = ['gold', 'royalblue']
    yg = ['gold', 'limegreen']
    pb = ['darkorchid', 'royalblue']
    pg = ['darkorchid', 'limegreen']
    bg = ['royalblue', 'limegreen']

    # make 33%/33%/33% graphs (3 treatment combinations)
    pry = ['darkorchid', 'gold', 'orangered']
    bry = ['royalblue', 'gold', 'orangered']
    gry = ['limegreen', 'gold', 'orangered']
    prb = ['darkorchid', 'royalblue', 'orangered']
    prg = ['darkorchid', 'limegreen', 'orangered']
    brg = ('royalblue', 'limegreen', 'orangered')
    pyb = ('darkorchid', 'royalblue', 'gold')
    pyg = ('darkorchid', 'limegreen', 'gold')
    byg = ('royalblue', 'limegreen', 'gold')
    pgb = ('darkorchid', 'royalblue', 'limegreen')

    # make 25%/25%/25%/25% graphs (4 treatment combinations)
    pryb = ['darkorchid', 'royalblue', 'gold', 'orangered']
    pryg = ['darkorchid', 'limegreen', 'gold', 'orangered']
    bryg = ['royalblue', 'limegreen', 'gold', 'orangered']
    prgb = ['darkorchid', 'royalblue', 'limegreen', 'orangered']
    pygb = ['darkorchid', 'royalblue', 'limegreen', 'gold']

    # make 20%/20%/20%/20%/20% graph (5 treatment combination)
    prygb = ['darkorchid', 'royalblue', 'limegreen', 'gold', 'orangered']

    colors = [r, y, p, b, g, ry, rp, rb, rg, yp, yb, yg, pb, pg, bg, pry, bry, gry, prb, prg, brg, pyb, pyg, byg, pgb,
               pryb, pryg, bryg, prgb, pygb, prygb]
    # print(len(colors))
    # print(colors)

    return colors

class Circle(HandlerPatch):
    # note, the function name HAS to be create_artists (it comes with a pre-defined set of parameters)
    # this function is for generating the circles in the legend

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        r = 5
        x = r + width//2
        y = height//2
        patch = patches.Circle(xy=(x, y), radius=r)
        self.update_prop(patch, orig_handle, legend)
        patch.set_transform(trans)
        return [patch]

def pie_graphs_100(color):
    # generate pie graphs for x-axis to represent 1 single bacterial treatment
    sizes = [100]
    colors = [color]
    plt.pie(sizes, colors=colors, shadow=False, startangle=90, wedgeprops={"edgecolor":"0",'linewidth':5,'antialiased': True})
    # plt.show()

def pie_graphs_50(color1, color2):
    # generate pie graphs for x-axis to represent 2 bacterial treatments
    sizes = [50, 50]
    colors = [color1, color2]
    plt.pie(sizes, colors=colors, shadow=False, startangle=90, wedgeprops={"edgecolor":"0",'linewidth':2,'antialiased': True})
    # plt.show()

def pie_graphs_33(color1, color2, color3):
    # generate pie graphs for x-axis to represent 3 bacterial treatments
    sizes = [33.3, 33.3, 33.3]
    colors = [color1, color2, color3]
    plt.pie(sizes, colors=colors, shadow=False, startangle=90, wedgeprops={"edgecolor":"0",'linewidth':2,'antialiased': True})
    # plt.show()

def pie_graphs_25(color1, color2, color3, color4):
    # generate pie graphs for x-axis to represent 4 bacterial treatments
    sizes = [25, 25, 25, 25]
    colors = [color1, color2, color3, color4]
    plt.pie(sizes, colors=colors, shadow=False, startangle=90, wedgeprops={"edgecolor":"0",'linewidth':2,'antialiased': True})
    # plt.show()

def pie_graphs_20():
    # generate pie graphs for x-axis to represent 5 bacterial treatments
    sizes = [20, 20, 20, 20, 20]
    colors = ['darkorchid', 'royalblue', 'limegreen', 'gold', 'orangered']
    plt.pie(sizes, colors=colors, shadow=False, startangle=90, wedgeprops={"edgecolor":"0",'linewidth':2,'antialiased': True})
    # plt.show()

def generate_pie_graphs():
    # make pie graphs to plot on the x axis underneath its respective bar

    # make 100% graphs (1 single treatment)
    pie_graphs_100('orangered')
    pie_graphs_100('gold')
    pie_graphs_100('darkorchid')
    pie_graphs_100('royalblue')
    pie_graphs_100('limegreen')

    # make 50%/50% pie graphs (2 treatment combinations)
    pie_graphs_50('orangered', 'gold')
    pie_graphs_50('orangered', 'darkorchid')
    pie_graphs_50('orangered', 'royalblue')
    pie_graphs_50('orangered', 'limegreen')
    pie_graphs_50('gold', 'darkorchid')
    pie_graphs_50('gold', 'royalblue')
    pie_graphs_50('gold', 'limegreen')
    pie_graphs_50('darkorchid', 'royalblue')
    pie_graphs_50('darkorchid', 'limegreen')
    pie_graphs_50('royalblue', 'limegreen')

    # make 33%/33%/33% graphs (3 treatment combinations)
    pie_graphs_33('darkorchid', 'gold', 'orangered')
    pie_graphs_33('royalblue', 'gold', 'orangered')
    pie_graphs_33('limegreen', 'gold', 'orangered')
    pie_graphs_33('darkorchid', 'royalblue', 'orangered')
    pie_graphs_33('darkorchid', 'limegreen', 'orangered')
    pie_graphs_33('dodgerblue', 'limegreen', 'orangered')
    pie_graphs_33('darkorchid', 'royalblue', 'gold')
    pie_graphs_33('darkorchid', 'limegreen', 'gold')
    pie_graphs_33('royalblue', 'limegreen', 'gold')
    pie_graphs_33('darkorchid', 'royalblue', 'limegreen')

    # make 25%/25%/25%/25% graphs (4 treatment combinations)
    pie_graphs_25('darkorchid', 'royalblue', 'gold', 'orangered')
    pie_graphs_25('darkorchid', 'limegreen', 'gold', 'orangered')
    pie_graphs_25('royalblue', 'limegreen', 'gold', 'orangered')
    pie_graphs_25('darkorchid', 'royalblue', 'limegreen', 'orangered')
    pie_graphs_25('darkorchid', 'royalblue', 'limegreen', 'gold')

    # make 20%/20%/20%/20%/20% graph (5 treatment combination)
    pie_graphs_20()

    # make 0% pie graph (0 treatments)
    pie_graphs_100('white')

def main():

    f = open('FlygutCFUsData.csv')
    csv_f = csv.reader(f)

    rows = []
    for row in csv_f:
        rows.append(row)

    relative_abundance(rows)
    # generate_pie_graphs()

main()

