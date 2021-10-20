import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import itertools
import matplotlib.patches as patches
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter
import matplotlib.image as img
from matplotlib.offsetbox import OffsetBox, AnnotationBbox
import matplotlib.patches as mpatch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerPatch
from matplotlib.offsetbox import OffsetImage,AnnotationBbox
from decimal import Decimal

def predictions(fitness_trait_mean, fitness_trait_SE, binary_IDS, plot_index, y_label, fig_label1, fig_label2, sci_notation=False, rotation=False):
    # prepare subplot
    ax = plt.subplot(4, 2, plot_index)
    plt.subplots_adjust(top=.899, left=.086, bottom=.071, right=.973, hspace=.28, wspace=.268)

    # get a list of colors for each of the pie charts for each bacterial combination
    col = pie_colors()

    # MEASUREMENTS (figures A, C, D, F)
    # make x values to plot on graph so that data points on the graph do not overlap:
    x_vals = np.zeros(32)
    lastnum = 1
    index = np.zeros(32)
    for i in range(0, 6):
        # find the number of different possible combinations for each number of bacteria that could be present
        number_combinations = scipy.special.binom(5, i)
        if number_combinations > 1:
            spacing = np.linspace(-.4, .4, int(number_combinations))
            x_vals[lastnum:lastnum+spacing.size] = spacing
            index[lastnum:lastnum+spacing.size] = i
            lastnum = lastnum + spacing.size
        elif number_combinations ==1:
            pass
    index[lastnum] = 5

    x_vals = np.add(x_vals, index)

    # plot the data consisting of the collected measurements
    meas = ax.errorbar(x_vals[1:], fitness_trait_mean[1:], [1.96 * x for x in fitness_trait_SE[1:]],
                       capsize=2, lw=.8, ms=8, label='measured traits', zorder=3, ls='none', barsabove=False,
                       color='black')

    # for figure G only, transform the y-axis tick mark labels to be written in scientific notation
    if sci_notation is True:
        plt.tick_params(axis='both', which='major', labelsize=8)
        f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
        g = lambda x, pos: "${}$".format(f._formatSciNotation('%1.10e' % x))
        plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(g))
    else: # otherwise, y-axis tick mark labels are written as normally
        plt.tick_params(axis='both', which='major', labelsize=8)
    # for figure G only, rotate the y-axis tick mark labels
    if rotation is True:
        plt.yticks(rotation=45)

    # when using matplotlib's annotate function, when xycoord='axes fraction', this means that 0,0 is lower left of
    # axes and 1,1 is upper right of the axes (aka of the subplot or figure).
    # annotate figure label
    ax.annotate(fig_label1,
                xy=(-.15, 1.2), xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=15)

    # add vertical dashed line to clearly separate the data points belonging to a certain number
    # of bacterial species present
    for x in [1.5, 2.5, 3.5, 4.5]:
        ax.axvline(x, color='k', ls='--', lw=.5)
    ax.set_xlabel('Number of species present', fontsize=7)
    ax.set_ylabel(y_label, fontsize=7)
    # pre-determine the y axis tick marks
    if plot_index == 1:
        ax.set_yticks([2, 3, 4])
    if plot_index == 3:
        ax.set_yticks([10, 11])
    if plot_index == 5:
        ax.set_yticks([40, 45, 50, 55])
    if plot_index == 7:
        ax.set_yticks([2.5e5, 5.0e5, 7.5e5])

    # create pie chart markers for each data point. Pass the specific colors for each data point from colors nest list.
    # also indicate the ratio that each color will make up in the whole pie graph. Plot these pie chart markers where
    # the original data point would normally be
    for i in range(1, 6):
        generate_pie_chart_markers(xs=x_vals[i],
                      ys=fitness_trait_mean[i],
                      ratios=[1],
                      sizes=40,
                      colors=col[i-1],
                      ax=ax)
    for i in range(6, 16):
        generate_pie_chart_markers(xs=x_vals[i],
                      ys=fitness_trait_mean[i],
                      ratios=[.50, .50],
                      sizes=40,
                      colors=col[i-1],
                      ax=ax)
    for i in range(16, 26):
        generate_pie_chart_markers(xs=x_vals[i],
                      ys=fitness_trait_mean[i],
                      ratios=[.3333333, .3333333, .3333333],
                      sizes=40,
                      colors=col[i - 1],
                      ax=ax)
    for i in range(26, 31):
        generate_pie_chart_markers(xs=x_vals[i],
                      ys=fitness_trait_mean[i],
                      ratios=[.25, .25, .25, .25],
                      sizes=40,
                      colors=col[i - 1],
                      ax=ax)
    generate_pie_chart_markers(xs=x_vals[31],
                  ys=fitness_trait_mean[31],
                  ratios=[.20, .20, .20, .20, .20],
                  sizes=40,
                  colors=col[30],
                  ax=ax)

    # add the 5 colored pie for measured traits legend above figure A
    colors = ['prygb']
    if plot_index == 1:
        for index, c in enumerate(colors):
            # for x in range(0, 5):
            pie_annotation_box = offset_pie_image(c, ax)
            # print(ab)
            ax.add_artist(pie_annotation_box)
    # add 'measured traits' to the measured traits legend above figure A
    if plot_index == 1:
        ax.annotate("measured traits", xy=(.45, 1.065), xycoords='axes fraction', fontsize=7)

    # make legend denoting the bacterial species that is associated with each color (this is the topmost legend)
    if plot_index == 1:
        circ1 = patches.Circle((0, 0), 1, facecolor='orangered')
        circ2 = patches.Circle((0, 0), 1, facecolor='gold')
        circ3 = patches.Circle((0, 0), 1, facecolor='darkorchid')
        circ4 = patches.Circle((0, 0), 1, facecolor='royalblue')
        circ5 = patches.Circle((0, 0), 1, facecolor='limegreen')

        plt.legend((circ1, circ2, circ3, circ4, circ5),
                   ('$\it{Lactobacillus\,plantarum}$',
                    '$\it{Lactobacillus\,brevis}$',
                    '$\it{Acetobacter\,pasteurianus}$',
                    '$\it{Lactobacillus\,tropicalis}$',
                    '$\it{Lactobacillus\,orientalis}$'),
                   ncol=5,
                   handler_map={patches.Circle: Circle(),},
                   fontsize=7,
                   loc='upper left',
                   bbox_to_anchor=(0.1, 1.48),
                   columnspacing=0.5)


    ### PREDICTIONS ###

    # FIND THE SINGLE SPECIES PREDICTIONS:
    # predict the traits of higher diversity bacterial combinations using the avg of single species diversity
    # for example, follow: f(10001) = 1/2 (f(10000) + f(00001))
    # find the average of all of the single element combinations and compare it to the actual value
    single_value_prediction = []
    single_value_prediction_err = []
    k = 0
    first_bins = fitness_trait_mean[1:6] # use experimental measurements in which a single species was present
    first_bins_err = fitness_trait_SE[1:6] # use experimental measurements SE in which a single species was present

    for bin in binary_IDS[6:]:
        tot = 0
        for bnum in bin:
            tot += int(bnum)
        if tot > 1:
            pred = 0
            err_pred = 0
            j = 0
            for bnum in bin:
                if int(bnum) == 1:
                    pred += first_bins[j]
                    err_pred += first_bins_err[j]**2
                j += 1
            pred = pred/tot
            err_pred = err_pred/tot

        single_value_prediction.append(pred)
        single_value_prediction_err.append((err_pred+fitness_trait_SE[k+6]**2)**.5)
        k += 1

    # calculate the difference b/w prediction and raw measurement
    SV_pred = np.subtract(np.array(single_value_prediction),fitness_trait_mean[6:])
    # calculate the standard error
    SV_error = np.array(single_value_prediction_err)

    # make plot for single species prediction
    ax = plt.subplot(4, 2, plot_index+1)
    pred_graph = ax.errorbar(x_vals[6:], SV_pred, [1.96 * x for x in SV_error], fmt='D',
                          capsize=2, lw=.8, ms=4, mfc='white', mew=1, color='gray',
                          label='single species', zorder=2, barsabove=True)
    ax.set_xticks([2, 3, 4, 5]) # set the x tick marks to always be 2, 3, 4, 5

    # pre-determine the y axis tick marks
    if plot_index == 1: # for figure B
        ax.set_yticks([-2, 0, 2])
    if plot_index == 3: # for figure D
        ax.set_yticks([-0.5, 0, 0.5, 1.0])
    if plot_index == 5: # for figure F
        ax.set_yticks([-5, 0, 5, 10])
    if plot_index == 7: # for figure H
        # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        # a = float("{:.1e}".format(-500000))
        # a = plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{-500000:,.1f}'))  # 2 decimal places
        ax.set_yticks([-5e5, -2.5e5, 0])
    ax.yaxis.set_tick_params(labelsize=7)
    axes = plt.gca()

    # get the minimum y coordinate on the y-axis. This will be useful for plotting the pie charts along the lowest
    # y coordinate along the x-axis
    y_min, y_max = axes.get_ylim()

    # create pie chart markers for each data point. Pass the specific colors for each data point from colors nest list.
    # also indicate the ratio that each color will make up in the whole pie graph. Plot these pie chart markers at the
    # bottom of the graph along the x-axis
    for i in range(6, 16): # bacterial combinations 6-15
        generate_pie_chart_markers(xs=x_vals[i],
                      ys=y_min-.4,
                      ratios=[.50, .50],
                      sizes=20,
                      colors=col[i-1],
                      ax=ax)
    for i in range(16, 26): # bacterial combinations 16-25
        generate_pie_chart_markers(xs=x_vals[i],
                      ys=y_min-.4,
                      ratios=[.3333333, .3333333, .3333333],
                      sizes=20,
                      colors=col[i - 1],
                      ax=ax)
    for i in range(26, 31): # bacterial combinations 26-30
        generate_pie_chart_markers(xs=x_vals[i],
                      ys=y_min-.4,
                      ratios=[.25, .25, .25, .25],
                      sizes=20,
                      colors=col[i - 1],
                      ax=ax)
    generate_pie_chart_markers(xs=x_vals[31], # bacterial combination 31 (germ free)
                  ys=y_min-.4,
                  ratios=[.20, .20, .20, .20, .20],
                  sizes=20,
                  colors=col[30],
                  ax=ax)


    # FIND THE SPECIES PAIRS PREDICTIONS:
    # for example, f(10011) = (1/3)(f(10010) + f(10001) + f(00011))
    # calculate the prediction which is essentially equal to the sum of the two-paired combinations
    two_value_prediction = []
    two_value_prediction_err = []
    k = 0
    first_bins = fitness_trait_mean[6:16]
    first_bins_err = fitness_trait_SE[6:16]
    for bin in binary_IDS[16:]:
        tot = 0
        for bnum in bin:
            tot += int(bnum)
        if tot > 2: # if the bacterial combination has 3 or more bacterial species
            pred = 0
            err_pred = 0
            j = 0
            indices = []
            i = 0
            for bnum in bin: # for every bit in a binary combo...
                if int(bnum) == 1: # if bit equals 1, append its index to indices
                    indices.append(i)
                i += 1
            # print("indices:")
            # print(indices)
            # use the indices of the position of each bit that equals 1 in a binary combo
            # find the combinations of the pairs of bits that make up the given indices
            # for example, if the indices was [0, 1, 2], the inner pairs are [(0, 1), (0, 2), (1, 2)]
            # this means that the binary code is 11000, 10100, and 01100
            inner_pairs = list(itertools.combinations(indices, 2))
            # print('inner pairs')
            # print(inner_pairs)

            # determine what binary ID the inner pair index is referring to
            for pair in inner_pairs:
                j = 0
                for bin in binary_IDS[6:16]: # loop through the binary codes that have two 1's (or 2 bacterial species present)
                    l = 0
                    check = 0

                    for bnum in bin: # for each bit in the binary code
                        # check the first number in the inner pair
                        # if the first number in the inner pair is equal to l (which represents the position of the bit
                        # in the binary ID) AND that specific position is equal to 1, add 1 to check
                        # for example, for the inner pair (1,3), 1 is the first number in the pair that is considered.
                        # 1 corresponds to the 2nd position in the binary code. For example, if the code is 01001, the
                        # 2nd position corresponds to 1. Also, determine if the bit in the second position in the
                        # binary code is equal to 1. Since the binary code is 01001, it does equal 1
                        if l == pair[0] and int(bnum) == 1:
                            check += 1
                        # using the second number (aka the index in the binary code) in the inner pair, determine the
                        # bit in this position of the binary code is equal to 1. If it does, add 1 to check
                        if l == pair[1] and int(bnum) == 1:
                            check += 1
                        l += 1

                    # once check = 2, this means that there are two 1's in the binary code. We have reached a pair.
                    # use the binary code at index j
                    if check == 2:
                        pred = pred + first_bins[j]
                        err_pred += first_bins_err[j]**2
                    j = j + 1

            pred = pred * 2 / (tot*(tot-1))
            err_pred = err_pred / tot
        two_value_prediction.append(pred)
        two_value_prediction_err.append((err_pred + fitness_trait_SE[k + 16] ** 2) ** .5)

        k = k + 1
    # print(two_value_prediction)

    TV_pred = np.subtract(np.array(two_value_prediction), fitness_trait_mean[16:])
    TV_error = np.array(two_value_prediction_err)

    ax = plt.subplot(4, 2, plot_index + 1)

    # adjust the x values so that the data points from the two-species pair predictions don't overlap with the single
    # species pair predictions
    jitter = np.repeat(.05, 32-16)

    # plot the two-species pair predictions
    pred_graph = ax.errorbar(np.add(x_vals[16:], jitter), TV_pred, [1.96 * x for x in TV_error], fmt='x',
                             capsize=2, lw=.5, ms=4, mfc='black', mew=.8, color='black',
                             label='species pairs', zorder=3)

    if sci_notation is True: # scientific notation on y-axis of graph H only
        plt.tick_params(axis='both', which='major', labelsize=8)
        f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
        g = lambda x, pos: "${}$".format(f._formatSciNotation('%1.10e' % x))
        plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(g))
    else:
        plt.tick_params(axis='both', which='major', labelsize=8)
    if rotation is True: # rotation of tick labels on y-axis of graph H only
        plt.yticks(rotation=45)

    # add figure label (B, D, F, G)
    ax.annotate(fig_label2,
                xy=(-.17, 1.2), xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=15)

    # add vertical dashed line to clearly separate the data points belonging to a certain number
    # of bacterial species present
    for x in [2.5, 3.5, 4.5]:
        ax.axvline(x, color='k', ls='--', lw=.5)
    ax.axhline(y=0, color='r', ls='--', lw=.5)
    ax.set_xlabel('Number of species present', fontsize=7)
    ax.set_ylabel(y_label, fontsize=7)

    # plot legend above figure 2B
    if plot_index == 1:
        ax.legend(title='prediction error: ', title_fontsize=7, loc='upper left', ncol=2, bbox_to_anchor=(0.18, 1.26),
                  fontsize=7, frameon=False, edgecolor='black', handletextpad=1, columnspacing=2)

    confidence_interval_capture(SV_pred, SV_error, TV_pred, TV_error, fig_label2)

def confidence_interval_capture(SV_pred, SV_error, TV_pred, TV_error, fig_label2):
    # CALCULATE HOW MANY PREDICTIVE MEASUREMENTS ARE WITHIN 95% CONFIDENCE INTERVAL

    single_species = 0 # the number of predictions captured for single species model
    pairwise_species = 0 # the number of predictions captured for pairwise species model
    SV_error = [1.96 * x for x in SV_error] # SE for single species model
    TV_error = [1.96 * x for x in TV_error] # SE for pairwise species model

    # determine number of predictions captured within 95% CI for single species model
    for prediction, err in zip(SV_pred[10:], SV_error[10:]):
        if (prediction - err) < 0 and (prediction + err) > 0:
            single_species += 1

    # determine number of predictions captured within 95% CI for pairwise model
    for prediction, err in zip(TV_pred, TV_error):
        if (prediction - err) < 0 and (prediction + err) > 0:
            pairwise_species += 1

    print('For graph ' + fig_label2 + ':')
    print('The number predictions using the single species model that were in the 95% confidence interval was '
          + str(single_species) + ' out of 16')
    print('The number predictions using the pairwise model that were in the 95% confidence interval was '
          + str(pairwise_species) + ' out of 16 ')

    # perform Fisher's Exact Test
    Fisher_exact_test(single_species, pairwise_species)

def Fisher_exact_test(single_species, pairwise_species):
    # the paper used Fisher's exact test in order to determine if there was a significant association
    # between the number of predictions captured in the 95% confidence interval between the single species prediction
    # and the pairwise species predictions
    # note, the data table used as to be a 2x2 contingency table

    print("Fisher's: ",
          + scipy.stats.fisher_exact([[single_species, pairwise_species],
                                   [16-single_species, 16-pairwise_species]])
          [1])

def get_pie_png(name):
    # open pie png images for all 32 pie charts
    path_to_file = "/Users/Sheila/Documents/BF550/Homework/Project1/pies/{}.png".format(name)
    image = plt.imread(path_to_file)
    return image

def offset_pie_image(name, ax):
    # using the pie png images, add each png image to an annotation box
    # used this function to create the 5 colored pie chart in the legend for measured traits above figure A

    image = get_pie_png(name)
    im = OffsetImage(image, zoom=.01)
    im.image.axes = ax

    pie_AB = AnnotationBbox(im, xy=(2.4, 4.859), xybox=(0, 10), xycoords='data',
                        frameon=False, boxcoords="offset points", pad=0)
    return pie_AB

def pie_colors():
    # indicate the color(s) used for each pie chart data point

    # 100% pie graphs (1 treatment combinations)
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
    byg = ('dodgerblue', 'limegreen', 'gold')
    pgb = ('darkorchid', 'royalblue', 'limegreen')

    # make 25%/25%/25%/25% graphs (4 treatment combinations)
    pryb = ['darkorchid', 'royalblue', 'gold', 'orangered']
    pryg = ['darkorchid', 'limegreen', 'gold', 'orangered']
    bryg = ['royalblue', 'limegreen', 'gold', 'orangered']
    prgb = ['darkorchid', 'royalblue', 'limegreen', 'orangered']
    pygb = ['darkorchid', 'royalblue', 'limegreen', 'gold']

    # make 20%/20%/20%/20%/20% graph (5 treatment combination)
    prygb = ['darkorchid', 'royalblue', 'limegreen', 'gold', 'orangered']

    colors = [r, y, p, b, g,
              ry, rp, rb, rg, yp, yb, yg, pb, pg, bg,
              pry, bry, gry, prb, prg, brg, pyb, pyg, byg, pgb,
              pryb, pryg, bryg, prgb, pygb, prygb]
    # print(len(colors))
    # print(colors)
    return colors

def generate_pie_chart_markers(xs, ys, ratios, sizes, colors, ax=None):
    # create pie chart markers for each data point

    a = np.pi/2
    pie_chart_markers = []

    # calculate the points of each of the pie sections
    for color, ratio in zip(colors, ratios):
        b = 2 * np.pi * ratio + a
        x = [0] + np.cos(np.linspace(a, b, 10)).tolist() + [0]
        y = [0] + np.sin(np.linspace(a, b, 10)).tolist() + [0]
        xy = np.column_stack([x, y])
        a = b
        # print("s:")
        # print(np.abs(xy).max()**2*np.array(sizes))
        pie_chart_markers.append({'s': np.abs(xy).max()**2*np.array(sizes),
                                  'marker': xy,
                                  'facecolor': color})

    # for each pie chart marker, scatter it in order to create the point
    for marker in pie_chart_markers:
        ax.scatter(xs, ys, **marker)

    return ax

class Circle(HandlerPatch):
    # create patches of circles for the top most legend for each of the five bacterial species
    # note, the function name HAS to be create_artists (it comes with a pre-defined set of parameters)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        r = 5
        x = r + width//2
        y = height//2

        # create the circle patch
        p = patches.Circle(xy=(x, y), radius=r)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)

        return [p]

def organize_data(fitness_traits):
    # the purpose of this function is to prepare and organize the data in Excel file for later usage

    # make the fitness_traits data into a numpy array
    fitness_traits_array = np.array(fitness_traits)
    # print(fitness_traits_array)
    # convert mean daily fecundity, mean time to death, mean development time, and mean bacterial CFUs into arrays
    mean_fecundity = fitness_traits_array[1:, 2].astype(np.float32)
    mean_death = fitness_traits_array[1:, 5].astype(np.float32)
    mean_development = fitness_traits_array[1:, 8].astype(np.float32)
    mean_CFU = fitness_traits_array[1:, 11].astype(np.float32)

    # convert standard errors for mean fecundity, mean death, mean development, and mean CFUs into separate arrays
    fecundity_SE = fitness_traits_array[1:, 3].astype(np.float32)
    death_SE = fitness_traits_array[1:, 6].astype(np.float32)
    development_SE = fitness_traits_array[1:, 9].astype(np.float32)
    CFU_SE = fitness_traits_array[1:, 13].astype(np.float32)

    # create a list of binary IDs. Binary IDs are representative of bacterial combination. The order of the binary code
    # represents the presence of the bacterial species in the following order -LP, LB, AP, AT, and AO. 1 means that the
    # bacterial species was present in the given combination, 0 means it was absent
    num_treatments = 32
    binary_IDS = []
    for i in range(num_treatments + 1):
        binary_IDS.append(fitness_traits[i][0])
    del binary_IDS[0]
    binary_IDS = np.array(binary_IDS)  # convert binary_IDS to array
    # print(binary_IDS)

    # compute the number of bacterial species present in each of the 32 bacterial combinations. Determine this by using
    # the binary code. This is essentially the N number of species
    N_bac_species = []
    for id in binary_IDS:
        sum = 0
        for num in id:
            sum += int(num)
        N_bac_species.append(sum)
    # print(N_bac_species)
    N_bac_species = np.array(N_bac_species)  # convert to array

    # pass along the arrays that were organized in this function to perform single species and species pair predictions
    # also make plots
    plt.figure(figsize=(9, 10)) #9, 15
    plt.rc('ytick', labelsize=7)

    fecundity = predictions(mean_fecundity, fecundity_SE, binary_IDS, 1, 'Mean fecundity (progeny/day/female)', 'A', 'B')
    development = predictions(mean_development, development_SE, binary_IDS, 3, 'Development time (days)','C', 'D')
    death = predictions(mean_death, death_SE, binary_IDS, 5, 'Time to death (days)', 'E', 'F')
    CFU = predictions(mean_CFU, CFU_SE, binary_IDS, 7, 'Bacterial load (CFUs)', 'G', 'H', True, True)

    plt.show()

def main():
    file = open('SummaryDataTable.csv', 'r')
    fitness_traits = []
    for line in file:
        fitness_traits.append(line.strip().split(","))

    organize_data(fitness_traits)

main()
