import pandas as pd
import matplotlib.pyplot as plt


CORRELATION_THRESHOLD = 0.5


def create_df_of_correlative_proteins(bac_names, correlation_file):
    df = pd.read_excel(correlation_file)
    df_selected = df[['protein_id', 'pearson_correlation']]

    correlated_prot_per_bac = {}


    for bacteria in bac_names:
        filtered_df = df_selected[
            (df_selected['protein_id'].apply(lambda x: '_'.join(x.split('_')[1:-1])) == bacteria) &
            (df_selected['pearson_correlation'] >= CORRELATION_THRESHOLD)]

        correlated_prot_per_bac[bacteria] = filtered_df['protein_id'].tolist()

    df_result = pd.DataFrame({
        'bacteria': list(correlated_prot_per_bac.keys()),
        'protein_ids': list(correlated_prot_per_bac.values())
    })

    # Save the DataFrame to an Excel file
    df_result.to_excel(f'pearson_correlations_proteins_threshold_{CORRELATION_THRESHOLD}.xlsx', index=False)

    create_plot_number_of_correlated(df_result)

    return df_result


def create_plot_number_of_correlated(correlated_prot_per_bac):
    bacteria_names = list(correlated_prot_per_bac["bacteria"])
    proteins_list = list(correlated_prot_per_bac["protein_ids"])
    counts = [len(proteins_list[i]) for i, bacteria in enumerate(bacteria_names)]

    # Sort the data by counts in descending order
    sorted_indices = sorted(range(len(counts)), key=lambda k: counts[k], reverse=True)
    bacteria_names_sorted = [bacteria_names[i] for i in sorted_indices]
    counts_sorted = [counts[i] for i in sorted_indices]

    # Plotting
    plt.figure(figsize=(10, 6))
    bars = plt.bar(range(len(bacteria_names_sorted)), counts_sorted, color='skyblue')
    plt.ylabel('Number of correlated Proteins')
    plt.title(f'Number of correlated Proteins (with correlation >= {CORRELATION_THRESHOLD}) Associated with Each Bacteria')
    plt.xticks([], [])  # Remove x-ticks
    plt.tight_layout()

    # Adding text labels on each bar
    for i, count in enumerate(counts_sorted):
        # Determine vertical alignment based on a threshold
        va = 'bottom' if count <= 70 else 'top'
        plt.text(i, count, bacteria_names_sorted[i], ha='center', va=va, rotation=90)

    plt.savefig(f"correlation_plot_threshold_{CORRELATION_THRESHOLD}.png")
    plt.show()
