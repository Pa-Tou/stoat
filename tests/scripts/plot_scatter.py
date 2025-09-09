import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from sklearn.linear_model import LinearRegression
import argparse
import gzip

parser = argparse.ArgumentParser(description="scatter plot")
parser.add_argument('input_file')
parser.add_argument('out_file')
parser.add_argument('--title', default = "Title", help="title")
parser.add_argument('--x_label', default ="", help = "x-axis label")
parser.add_argument('--y_label', default = "", help = "y-axis label")
parser.add_argument('--x_col', default = 0, type=int, help = "x column")
parser.add_argument('--y_col', default = 1, type=int, help = "y column")
parser.add_argument('--color_col', default = -1, type=int, help = "column for coloring")
parser.add_argument('--log_y', action='store_true', help = "log the y axis")

args = parser.parse_args()


if (args.input_file[-2:] == "gz"):
    f = gzip.open(args.input_file, 'rt')
else:
    f= open(args.input_file)
label_to_points = dict()
x_label = args.x_label
y_label = args.y_label
color_label = ""
big_count=0

for line in f:
    l = line.split("\t")
    if line[0] == "#":
        if x_label == "":
            x_label = l[args.x_col]
        if y_label == "":
            y_label = l[args.y_col]
        if args.color_col != -1:
            color_label = l[args.color_col] 
    else:
        if args.color_col != -1:
            label = l[args.color_col]
        else:
            label = ""
        if not label in label_to_points:
            label_to_points[label] = [[], []]
        label_to_points[label][0].append(l[args.x_col])
        label_to_points[label][1].append(l[args.y_col])

f.close()


fig_width = 12
fig_height = 10
plt.figure(figsize = (fig_width, fig_height))

panel_width = .8
panel_height = .8

panel1 = plt.axes([0.1, 0.1, panel_width, panel_height])
if args.log_y:
    plt.yscale('log')

colors = ["blue", "orange", "red", "yellow"]

color_i = 0
for label, points in label_to_points.items():
    count=len(points[0])
    plt.scatter(points[0], points[1], label=label+": "+str(count), alpha=0.7, color = colors[color_i])
    color_i += 1
if args.color_col != -1:
    panel1.legend(title = color_label+": count")

#panel1.tick_params(axis='both', which='both',\
#                    bottom='on',labelbottom='on',\
#                    left='on', labelleft='on',\
#                    right='off', labelright='off',\
#                    top='off', labeltop='off')   

panel1.set_xlabel(x_label)
panel1.set_ylabel(y_label)
panel1.set_title(args.title)
#panel1.set_yticks(range(0, score_hist[-1], score_hist[-1]//8))
#panel1.set_yticklabels([x for x in range(offset, score_hist[-1], score_hist[-1]//8)])

#plt.show()
plt.savefig(args.out_file, dpi=400)
