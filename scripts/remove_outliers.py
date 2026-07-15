import argparse
import math

def print_non_outliers(infile_name):

    # Get a sorted list of values
    values = []
    with open(infile_name, "r") as infile:
        for line in infile:
            l = line.split()
            if l[0] == "SAMPLE":
                continue
            values.append(float(l[1]))

    values.sort()

    quartile_size = math.floor(len(values)/4)
    if quartile_size <=1 :
        return
    q1 = values[quartile_size-1]
    q3 = values[quartile_size * 3 - 1]
    iqr = q3 - q1
    lower_cutoff = q1 - iqr * 1.5
    upper_cutoff = q3 + iqr * 1.5 

    with open(infile_name, "r") as infile:
        for line in infile:
            l = line.split()
            val = float(l[1])
            if val > lower_cutoff and val < upper_cutoff:
                print(line.strip())

def main():
    parser = argparse.ArgumentParser(description="filter outliers")
    parser.add_argument('input_file')
    args = parser.parse_args()

    print_non_outliers(args.input_file)




if __name__ == "__main__":
    main()

