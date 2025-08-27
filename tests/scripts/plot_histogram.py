import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#from sklearn.linear_model import LinearRegression
import argparse

parser = argparse.ArgumentParser(description="histogram")
parser.add_argument('input_file')
parser.add_argument('out_file')
parser.add_argument('--title', default = "Histogram", help="title")
parser.add_argument('--bins', default = "50", help="number of bins")
parser.add_argument('--column', default = "0", help="get values from this column")
parser.add_argument('--x_label', default ="", help = "x-axis label")
parser.add_argument('--y_label', default = "Count", help = "y-axis label")
parser.add_argument('--log_y', action='store_true', help="Log y-axis")
args = parser.parse_args()
print(args.column)


f= open(args.input_file)
numbers = []

for line in f:
    if (line[0] != "#"):
        l = line.split("\t")
        numbers.append(float(l[int(args.column)]))

print(args.bins)

#Color points by how correct they are
hist, bins = np.histogram(numbers, bins=int(args.bins))
print(hist)
print(bins)

fig_width = 12
fig_height = 10
plt.figure(figsize = (fig_width, fig_height))

panel_width = .8
panel_height = .8

panel1 = plt.axes([0.1, 0.1, panel_width, panel_height])
if (args.log_y): 
    plt.yscale('log')

# Make a legend for correctness
for i in range(int(args.bins)):

    rectangle = mplpatches.Rectangle([bins[i], 0],(bins[i+1]-bins[i]),hist[i],linewidth=0)
    panel1.add_patch(rectangle)


#panel1.tick_params(axis='both', which='both',\
#                    bottom='on',labelbottom='on',\
#                    left='on', labelleft='on',\
#                    right='off', labelright='off',\
#                    top='off', labeltop='off')   
panel1.set_xlim(min(bins), max(bins))
bin_ticks = range(0, int(args.bins), (int(args.bins)//6))
panel1.set_xticks([bins[bin_ticks[0]], bins[bin_ticks[1]], bins[bin_ticks[2]], bins[bin_ticks[3]], bins[bin_ticks[4]], bins[bin_ticks[5]]])

panel1.set_ylim(min(hist) * 0.9, max(hist) * 1.1)
panel1.set_ylabel(args.y_label)
panel1.set_xlabel(args.x_label)
mean = sum(numbers)/len(numbers)
median = np.median(numbers)
panel1.set_title (args.title)
panel1.text(min(bins),max(hist)*1.05,"Mean: {}\nMedian: {}".format(mean, median))

#panel1.set_yticks(range(0, score_hist[-1], score_hist[-1]//8))
#panel1.set_yticklabels([x for x in range(offset, score_hist[-1], score_hist[-1]//8)])

#plt.show()
plt.savefig(args.out_file, dpi=400)
