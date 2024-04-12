# Indigo Hartsell, 2024-04-12
# This script reads terpene data from Logseq pages and creates a graph of terpene scores and frequencies.

import os
import matplotlib.pyplot as plt
import seaborn as sns
import typing

PAGE_DIR = "S:\Logseq\Home Graph\pages"
EXCLUDE_FILES = ["templates.md", "README.md", "index.md"] # files to exclude from search

# Set the style for seaborn and matplotlib
sns.set_style("dark")
plt.style.use("dark_background")

# Adjust the background color
plt.rcParams['axes.facecolor'] = '#383838'  # Adjust the value to your desired shade
plt.rcParams['figure.facecolor'] = '#383838'

# Define a color palette with pastel colors
# mauve, tropical indigo, nyanza, cambridge blue, hooker's green
pastel_colors = ["#C2AFF0", "#9191E9", "#D6F8D6", "#7FC6A4", "#FFB994"]

extract_terps = lambda line: line.strip().split("::")[1].split(",") # extract terpenes from line

def find_terp_data(PATH):
    terp_data = {} # {terp1: score1, terp2: score2, ...}
    terp_scores = list() # [(terp1, score1), (terp2, score2), ...]
    for filename in os.listdir(PATH):
        if filename.endswith(".md") and filename not in EXCLUDE_FILES:
            filepath = os.path.join(PATH, filename)
            with open(filepath, "r", encoding="utf-8") as file:
                no_score = False
                for line in file:
                    if line.startswith("score::"):
                        str_strain_score = line.split("::")[1].strip()
                        if str_strain_score == '':
                            no_score = True
                            continue
                        strain_score = float(str_strain_score)
                    if line.startswith("terps::") and not no_score:
                        terps = [terp.strip() for terp in extract_terps(line)] # [terp1, terp2, ...]
                        for i, terp in enumerate(terps):
                            terp_score = strain_score / (len(terps) - i) # distribute score evenly among terpenes

                            terp_scores.append((terp, terp_score))
                            # add terp_score to terp_data running average
                            if terp in terp_data:
                                terp_data[terp] = (terp_data[terp] * (len(terp_scores) - 1) + terp_score) / len(terp_scores)
                            else:
                                terp_data[terp] = terp_score

    # remove '' key and '?' key
    terp_data = {k: v for k, v in terp_data.items() if k not in ["", "?"]}
    terp_scores = [(k, v) for k, v in terp_scores if k not in ["", "?"]]

    return terp_data, terp_scores

def create_graph(terp_data: dict, terp_scores: dict):
    # Extract unique terpenes and their scores
    unique_terps = list(set([terp for terp, _ in terp_scores]))
    
    # Extract scores from terp_data
    terp_scores_dict = {terp: [score for t, score in terp_scores if t == terp] for terp in unique_terps}

    # Calculate frequencies
    frequencies = {}
    for (terp, _) in terp_scores:
        frequencies[terp] = frequencies.get(terp, 0) + 1

    print(terp_data)

    # Sort the data based on mean scores
    sorted_terps = sorted(terp_data, key=terp_data.get)
    sorted_mean_scores = [terp_data[terp] for terp in sorted_terps]
    sorted_frequencies = [frequencies[terp] for terp in sorted_terps]

    # Create a figure with two rows and two columns
    fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), gridspec_kw={'width_ratios': [2, 1]})

    # Bar plot for terpene scores
    ax1.bar(sorted_terps, sorted_mean_scores, color=pastel_colors[0])
    ax1.set_xlabel('Terpene')
    ax1.set_ylabel('Mean Score')
    ax1.set_title('Terpene Scores', color=pastel_colors[2])

    # Line plot for terpene frequencies
    sns.lineplot(x=sorted_terps, y=sorted_frequencies, ax=ax2, color=pastel_colors[0])
    ax2.set_xlabel('Terpene')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Terpene Frequencies', color=pastel_colors[2])

    # Display raw terpene data
    ax3.axis('off')
    ax4.axis('off')
    raw_data = [f"{terp}: {', '.join(map(str, scores))}" for terp, scores in terp_scores_dict.items()]
    combined_text = '\n'.join(raw_data)
    ax3.text(0, 1, combined_text, va='top', ha='left', fontfamily='monospace', color=pastel_colors[4])
    ax3.set_title('Raw Terp Score Points', color=pastel_colors[2])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    terp_data, terp_scores = find_terp_data(PAGE_DIR)

    create_graph(terp_data, terp_scores)
