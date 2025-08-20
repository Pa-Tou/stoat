import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import re
import os
from cyvcf2 import VCF # type: ignore

def split_snarl(input_str):
    # Split the std::string and filter empty elements, then convert to integers
    return [str(num) for num in re.split(r'[><]', input_str) if num]

def parse_sv_rows(file_path):
    sv_indices_set = set()
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip().split('\t')
    type_index = header.index('TYPE')
    
    for i, line in enumerate(lines[1:], start=1):  # start=1 to get correct line index (excluding header)
        fields = line.strip().split('\t')
        type_field = fields[type_index]
        
        # Handle multiple values separated by commas
        values = []
        for part in type_field.split(','):
            for subpart in part.split('/'):  # in case of values like 508/508
                try:
                    values.append(int(subpart))
                except ValueError:
                    pass
        
        # Check absolute differences between all pairs
        for j in range(len(values)):
            for k in range(j + 1, len(values)):
                if abs(values[j] - values[k]) >= 50:
                    sv_indices_set.add(i)
                    break  # no need to check more pairs for this row
            else:
                continue
            break
    
    return sv_indices_set

# Function to process the frequency file and get result list with differences
def process_file(freq_file, threshold=0.2):
    df = pd.read_csv(freq_file, sep='\t')
    freq_path_list = []
    true_labels = []
    list_diff = []
 
    # Iterate in pairs and calculate the differences
    for i in range(0, len(df) - 1, 2):
        row_1 = df.iloc[i]
        row_2 = df.iloc[i + 1]
        diff = abs(row_1['freq'] - row_2['freq'])

        # If group freq diff > threshold, then true label is 0, else 1
        if diff > threshold:
            true_labels.append(0)
        else :
            true_labels.append(1)

        freq_path_list.append(f"{int(row_1['start_node'])}_{int(row_1['next_node'])}")
        list_diff.append(float(diff))

    return freq_path_list, true_labels, list_diff

def check_valid_snarl(start_node_1, next_node_1, start_node_2, next_node_2, snarl_list):
    """
    Check if a valid snarl exists in the snarl_list.

    A valid snarl is defined as:
    - At least one path contains both start_node_1 and next_node_1.
    - At least one path contains both start_node_2 and next_node_2.
    - start_node_1 and start_node_2 must be the same str.
    - One path can not containt the 4 node.
    """

    if start_node_1 != start_node_2:  # First condition: start nodes must match
        raise ValueError("start_node_1 != start_node_2")

    contains_first_pair = False
    contains_second_pair = False

    for path in snarl_list:
        # Extract numbers from the path as integers
        decomposed_path = {int(num) for num in re.findall(r'\d+', path)}
        
        # Check if the current path contains both nodes of the first pair
        if int(start_node_1) in decomposed_path and int(next_node_1) in decomposed_path:
            contains_first_pair = True

        # Check if the current path contains both nodes of the second pair
        if int(start_node_2) in decomposed_path and int(next_node_2) in decomposed_path:
            contains_second_pair = True

        # If both conditions are satisfied, we can exit early
        if contains_first_pair and contains_second_pair:
            return True

    # Return true only if both pairs are satisfied
    return contains_first_pair and contains_second_pair

def match_snarl(freq_path_list, true_labels, list_diff, p_value_file, paths_file, save_sv_snarl, type_):
    # Read both files
    p_value_df = pd.read_csv(p_value_file, sep='\t')
    paths_df = pd.read_csv(paths_file, sep='\t')[['SNARL', 'PATHS']]

    # Merge on SNARL to align PATHS with p-values
    p_value_df = p_value_df.merge(paths_df, on='SNARL', how='left')

    # Pre-split SNARL into numeric parts for vectorized matching
    snarl_split = p_value_df['SNARL'].str.split('_', expand=True).astype(int)
    snarl_start = snarl_split[0].values
    snarl_end = snarl_split[1].values

    # Outputs
    predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8 = [], [], []
    cleaned_true_labels, clean_list_diff, pvalue_list, num_sample, snarl_name = [], [], [], [], []

    step = 2
    for i in range(0, len(freq_path_list) - 1, step):
        start_node_1, next_node_1 = map(int, freq_path_list[i].split('_'))
        start_node_2, next_node_2 = map(int, freq_path_list[i + 1].split('_'))

        # Vectorized match filtering
        match_mask = ((snarl_end <= start_node_1) & (snarl_start >= next_node_1)) | \
                     ((snarl_start <= start_node_1) & (snarl_end >= next_node_1))
        matched_row = p_value_df[match_mask]

        if not matched_row.empty:
            if save_sv_snarl is not None:
                matched_row = matched_row[matched_row['SNARL'].isin(save_sv_snarl)]

            for _, matched in matched_row.iterrows():
                if check_valid_snarl(start_node_1, next_node_1, start_node_2, next_node_2, matched['PATHS'].split(',')):
                    # Choose p-value field
                    if type_ == 'binary':
                        p_val = matched['P_FISHER']
                        group_paths = matched.get('GROUP_PATHS', '')
                        allele_num = sum(
                            int(part.split(':')[1]) for part in group_paths.split(',') if ':' in part) if group_paths else 200
                    elif type_ == 'quantitative':
                        p_val = matched['P']
                        group_paths = matched.get('ALLELE_PATHS', '')
                        allele_num = sum(int(x) for x in group_paths.split(',')) if group_paths else 200
                    else:
                        raise ValueError("type_ must be 'binary' or 'quantitative'")

                    snarl_name.append(matched['SNARL'])
                    predicted_labels_10_2.append(0 if p_val < 0.01 else 1)
                    predicted_labels_10_5.append(0 if p_val < 1e-5 else 1)
                    predicted_labels_10_8.append(0 if p_val < 1e-8 else 1)
                    cleaned_true_labels.append(true_labels[i])
                    clean_list_diff.append(list_diff[i])
                    pvalue_list.append(p_val)
                    num_sample.append(allele_num)

    return (
        predicted_labels_10_2,
        predicted_labels_10_5,
        predicted_labels_10_8,
        cleaned_true_labels,
        clean_list_diff,
        pvalue_list,
        num_sample,
        snarl_name
    )

def conf_mat_maker(p_val, predicted_labels, true_labels, output):
        
    # Calculate confusion matrix for p-value < p_val
    print(f"\nMetrics for p-value < {p_val}:")

    # Inverse because I want the X axis to be the truth labels and Y axis to be the predicted labels
    cm = confusion_matrix(predicted_labels, true_labels)
    print(f"Confusion stoat_vcf::EdgeBySampleMatrix for p-value < {p_val}:\n{cm}")
    prec = precision_score(predicted_labels, true_labels)
    recall = recall_score(predicted_labels, true_labels)
    f1 = f1_score(predicted_labels, true_labels)
    print(f"Precision: {prec:.3f}")
    print(f"Recall: {recall:.3f}")
    print(f"F1 Score: {f1:.3f}")

    # Plot confusion matrix for p-value < p_val
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Positive', 'Negative'], 
                yticklabels=['Positive', 'Negative'],
                annot_kws={"size": 30})
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  
    plt.title(f'Confusion stoat_vcf::EdgeBySampleMatrix for p-value < {p_val}', fontsize=18)  # Increase title font size
    plt.xlabel('Truth Labels', fontsize=20)  # Increase x-label font size
    plt.ylabel('Predicted Labels', fontsize=20)  # Increase y-label font size
    plt.savefig(output + f'_{p_val}.png', format='png', dpi=300)

def print_confusion_matrix(predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, true_labels, output):
    
    p_val_10_2 = 0.01
    p_val_10_5 = 0.00001
    p_val_10_8 = 0.00000001

    conf_mat_maker(p_val_10_2, predicted_labels_10_2, true_labels, output)
    conf_mat_maker(p_val_10_5, predicted_labels_10_5, true_labels, output)
    conf_mat_maker(p_val_10_8, predicted_labels_10_8, true_labels, output)

def p_value_distribution(test_predicted_labels, cleaned_true_labels, list_diff, p_value, num_sample, snarl_name, output):
    
    false_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 1]

    print("len(false_positive_indices) : " , len(false_positive_indices))

    true_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 0]
    
    print("len(true_positive_indices) : " , len(true_positive_indices))

    false_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 0]

    print("len(false_negative_indices) : " , len(false_negative_indices))

    diff_false_positive = [list_diff[i] for i in false_positive_indices]
    pvalue_false_positive = [p_value[i] for i in false_positive_indices]
    minsample_false_positive = [num_sample[i] for i in false_positive_indices]
    snarl_name_false_positive = [snarl_name[i] for i in false_positive_indices]

    diff_true_positives = [list_diff[i] for i in true_positive_indices]
    pvalue_true_positives = [p_value[i] for i in true_positive_indices]
    minsample_true_positives = [num_sample[i] for i in true_positive_indices]
    snarl_name_true_positives = [snarl_name[i] for i in true_positive_indices]

    diff_false_negative = [list_diff[i] for i in false_negative_indices]
    pvalue_false_negative = [p_value[i] for i in false_negative_indices]
    minsample_false_negative = [num_sample[i] for i in false_negative_indices]
    snarl_name_false_negative = [snarl_name[i] for i in false_negative_indices]

    # Create a DataFrame for easy plotting
    data = {
        'P-Value': pvalue_false_positive + pvalue_true_positives + pvalue_false_negative,
        'Difference': diff_false_positive + diff_true_positives + diff_false_negative,
        'Min Sample': minsample_false_positive + minsample_true_positives + minsample_false_negative,
        'Type': ['False Positives'] * len(pvalue_false_positive) + ['True Positives'] * len(pvalue_true_positives) + ['False Negatives'] * len(pvalue_false_negative),
        'Snarl': snarl_name_false_positive + snarl_name_true_positives + snarl_name_false_negative
    }

    df = pd.DataFrame(data)

    # Create the interactive scatter plot
    fig = px.scatter(
        df, 
        x='P-Value', 
        y='Difference', 
        size='Min Sample', 
        color='Type',
        hover_name=df.index,
        hover_data={
            "P-Value": True,  # Include P-Value in hover box
            "Difference": True,  # Include Difference in hover box
            "Min Sample": True,  # Include Min Sample (size is already shown)
            "Snarl": True,  # Include Snarl name 
        },
        title="Distribution of P-Values for False Positives and True Positives",
        labels={"P-Value": "P-Value", "Difference": "Simulated Effect (Difference in Probabilities)"},
        size_max=20
    )

    fig.update_layout(
        xaxis_title="P-Value",
        yaxis_title="Simulated Effect (Difference in Probabilities)",
        legend_title="Type",
        template="plotly_white",
        xaxis_title_font=dict(size=30),  # Increase font size for x-axis title
        yaxis_title_font=dict(size=30),  # Increase font size for y-axis title
        xaxis=dict(tickfont=dict(size=25)),  # Increase font size for x-axis ticks
        yaxis=dict(tickfont=dict(size=25)),  # Increase font size for y-axis ticks
    )

    # Show the interactive plot
    fig.show()

    # Optionally, save the plot as an HTML file
    fig.write_html(f'{output}_pvalue_interactive.html')

    # Save as static PNG
    fig.write_image(f'{output}_pvalue_interactive.png')

def plot_diff_distribution(test_predicted_labels:list, cleaned_true_labels:list, clean_list_diff:list, output:str, pvalue:list):

    # True label 0 : difference freq significative
    # True label 1 : No difference freq significative
    # Predict label 0 : pvalue significative
    # Predict label 1 : No pvalue significative

    # ----------------------------- True Positive -----------------------------
    true_positive_indices= [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 0]

    true_positive_diffs = [clean_list_diff[i]*100 for i in true_positive_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(true_positive_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True Positive", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_true_positive.png', format='png', dpi=300)

    # ----------------------------- False Positive -----------------------------
    false_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 1]

    false_positive_diffs = [clean_list_diff[i]*100 for i in false_positive_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(false_positive_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Positive", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_false_positive.png', format='png', dpi=300)

    # ----------------------------- False Negative -----------------------------
    false_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 0]

    false_negative_diffs = [clean_list_diff[i]*100 for i in false_negative_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(false_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Negative", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_false_negative.png', format='png', dpi=300)

    # ----------------------------- True Negative -----------------------------
    true_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 1]

    true_negative_diffs = [clean_list_diff[i]*100 for i in true_negative_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(true_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True Negatives", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_true_negatives.png', format='png', dpi=300)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse a file with start_node, next_node, group, and freq columns.")
    parser.add_argument("--freq", help="Path to allele frequence file")
    parser.add_argument("--p_value", help="Path to p_value gwas output file")
    parser.add_argument("--paths", help="Path to snarl paths list file")
    parser.add_argument("-t", "--threshold", type=float, required=False, help="Threshold to define the truth label")
    parser.add_argument("--sv", action="store_true", required=False, help="Verify truth only for SV (>= 50 pb diff)")
    parser.add_argument("--output", default="output", required=False, help="Verify truth only for SV (>= 50 pb diff)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", action='store_true', help="binary test")
    group.add_argument("-q", "--quantitative", action='store_true', help="quantitative test")

    args = parser.parse_args()

    # Ensure output directory exists
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    else:
        print(f"Output directory already exists: {output_dir}")

    if args.binary:
        type_ = 'binary'

    elif args.quantitative:
        type_ = 'quantitative'

    # THRESHOLD_FREQ = 0.0 : Case where just a difference between both group snarl is considered like Truth label
    THRESHOLD_FREQ = 0.0
    
    save_sv_snarl = parse_sv_rows(args.p_value) if args.sv else None
    freq_test_path_list, test_true_labels, test_list_diff = process_file(args.freq, THRESHOLD_FREQ)
    test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, clean_list_diff, pvalue, num_sample, snarl_name = match_snarl(freq_test_path_list, test_true_labels, test_list_diff, args.p_value, args.paths, save_sv_snarl, type_)
    test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, clean_list_diff, pvalue, num_sample, snarl_name = match_snarl(freq_test_path_list, test_true_labels, test_list_diff, args.p_value, args.paths, save_sv_snarl, type_)
    print_confusion_matrix(test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, f"{output_dir}/confusion_matrix_{THRESHOLD_FREQ}")
    assert len(cleaned_true_labels) == len(clean_list_diff)

    print("Pourcentage of paths tested : ", (len(pvalue)/len(freq_test_path_list))*100*2) # *2 because we jump 2 per 2 the paths
    # Plot distribution of p-values for false negatives and true positives
    plot_diff_distribution(test_predicted_labels_10_2, cleaned_true_labels, clean_list_diff, output_dir, "10^-2")
    p_value_distribution(test_predicted_labels_10_2, cleaned_true_labels, clean_list_diff, pvalue, num_sample, snarl_name, output_dir)

    """
    python3 tests/scripts/verify_truth.py --freq data/quantitative/pg.snarls.freq.tsv \
    --p_value output/quantitative_table_vcf.tsv --paths output/snarl_analyse.tsv -q

    python3 tests/scripts/verify_truth.py --freq data/quantitative/pg.snarls.freq.tsv \
    --p_value output/quantitative_table_covar_vcf.tsv --paths output/snarl_analyse.tsv -q

    python3 tests/scripts/verify_truth.py --freq data/binary/pg.snarls.freq.tsv \
    --p_value output/binary_table_vcf.tsv --paths output/binary_snarl_analyse.tsv -b

    python3 tests/scripts/verify_truth.py --freq data/binary/pg.snarls.freq.tsv \
    --p_value output/binary_table_covar_vcf.tsv --paths output/binary_snarl_analyse.tsv -q
    """
