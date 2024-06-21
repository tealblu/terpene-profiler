# Indigo Hartsell, 2024-04-12
# This script reads terpene data from Logseq pages and creates a graph of terpene scores and frequencies.

# Standard library imports
from datetime import datetime
import os
from collections import defaultdict
from itertools import combinations

# Visualization libraries
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import seaborn as sns

# Math and data manipulation libraries
import pandas as pd
import networkx as nx
import numpy as np

PAGE_DIR = r"S:\Logseq\Home Graph\pages"
EXCLUDE_FILES = [ # files to exclude from search
    "templates.md",
    "README.md",
    "index.md",
]

# Set the style for seaborn and matplotlib
sns.set_style("dark")
plt.style.use("dark_background")

# Adjust the background color
bgColor = "#383838"
plt.rcParams["axes.facecolor"] = bgColor # Adjust the value to your desired shade
plt.rcParams["figure.facecolor"] = bgColor

# Define a color palette with pastel colors
# mauve, tropical indigo, nyanza, cambridge blue, hooker's green
pastel_colors = ["#C2AFF0", "#9191E9", "#D6F8D6", "#7FC6A4", "#f79b6a"]

extract_terps = (
    lambda line: line.strip().split("::")[1].split(",")
)  # extract terpenes from line


class Terpene: # class to store terpene data
    def __init__(self, name: str, score: float, position: int, position_limit: int, strain: str = ""):
        self.name = name
        self.score = score
        self.position = position
        self.position_limit = position_limit
        self.strain = strain

    def __repr__(self) -> str:
        return f"{self.name}: {self.score}, weight {self.position}/{self.position_limit}"

def find_terp_data(path: str) -> tuple[list[Terpene], list, list]:
    terp_data = list() # list of terpene data instances
    terp_co_occurrences = list() # list of terpene data for each strain
    dates = list() # list of dates and scores for each strain

    for filename in os.listdir(path):
        if filename.endswith(".md") and filename not in EXCLUDE_FILES:
            filepath = os.path.join(path, filename)
            with open(filepath, "r", encoding="utf-8") as file:
                for line in file:
                    strain_name = filename[:-3] # extract strain name from filename

                    # dates can be in the following forms:
                    #   date:: [[Mar 2nd, 2024]]
                    #   date:: [[Mar 3rd, 2024]] 15:00
                    #   date:: [[Mar 4th, 2024]] 15:00

                    if line.startswith("date::"): # extract date
                        date_str = line.split("::")[1].strip()
                        if date_str == "": # if no date is found, skip terp extraction
                            continue
                        date_str = date_str.replace("[[", "").replace("]]", "") # remove brackets
                        date_str = date_str.replace("st", "").replace("nd", "").replace("rd", "").replace("th", "")
                        try:
                            date = datetime.strptime(date_str, "%b %d, %Y")
                        except ValueError:
                            date = datetime.strptime(date_str, "%b %d, %Y %H:%M")

                    if line.startswith("score::"):
                        str_strain_score = line.split("::")[1].strip() # extract score
                        if str_strain_score == "": # if no score is found, skip terp extraction
                            continue
                        strain_score = float(str_strain_score)

                        dates.append((date, strain_score))

                    if line.startswith("terps::"): # extract terps if score is found
                        terps = [
                            terp.strip() for terp in extract_terps(line)
                        ]  # [terp1, terp2, ...]
                        for i, terp in enumerate(terps): # add terpene data to list if terp is not "" or "?"
                            if terp != "" and terp != "?": terp_data.append(Terpene(terp, strain_score, i + 1, len(terps), filename[:-3]))

                        # add terps to terp_co_occurrences list
                        if "?" not in terps: 
                            terp_co_occurrences.append(terps)

    print("\nUnique terpenes:")
    for terp in set(terp.name for terp in terp_data):
        print(terp)

    print("\nData points:")
    for terp in terp_data:
        print(terp)

    return terp_data, terp_co_occurrences, dates

def weighted_score_average(terp_data: list[Terpene]) -> dict[str, float]:
    # Initialize a defaultdict to store the scores, positions, and position_limits for each terpene name
    terp_scores = defaultdict(lambda: [[], [], []])

    # Loop through the terpene data and populate the terp_scores dict
    for terp in terp_data:
        terp_scores[terp.name][0].append(terp.score)
        terp_scores[terp.name][1].append(terp.position)
        terp_scores[terp.name][2].append(terp.position_limit)

    # Calculate the weighted average score for each terpene name
    weighted_avg_scores = {}
    for name, (score_list, position_list, position_limit_list) in terp_scores.items():
        weighted_sum = 0
        total_weight = 0
        for score, position, position_limit in zip(score_list, position_list, position_limit_list):
            weight = (position_limit - position + 1) / position_limit
            weighted_sum += score * weight
            total_weight += weight
        weighted_avg_scores[name] = weighted_sum / total_weight if total_weight else 0

    # Print the weighted average scores to the console
    print("\nWeighted average scores:")
    for name, score in weighted_avg_scores.items():
        print(f"{name}: {score:.2f}")

    return weighted_avg_scores


def create_graph(terp_data: list[Terpene], terp_scores: dict, terp_co_occurrences: list[tuple], dates: list):
    # Create a figure with two rows and two columns
    fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(
        nrows=2, ncols=2, figsize=(12, 8), gridspec_kw={"width_ratios": [2, 1]}
    )

    # Extract unique terpenes
    unique_terps = set(terp.name for terp in terp_data)

    # Calculate frequencies
    frequencies = {terp: sum(1 for t in terp_data if t.name == terp) for terp in unique_terps}

    # Sort the data based on mean scores
    sorted_terps = sorted(unique_terps, key=lambda x: terp_scores[x])
    sorted_frequencies = [frequencies[terp] for terp in sorted_terps]

    # Create a DataFrame for the box-and-whisker plot
    box_plot_data = [(terp.name, terp.score) for terp in terp_data]
    box_plot_df = pd.DataFrame(box_plot_data, columns=['Terpene', 'Score'])
    box_plot_df['Mean_Score'] = box_plot_df['Terpene'].map(terp_scores)
    box_plot_df = box_plot_df.sort_values(by='Mean_Score')

    # Box-and-whisker plot for terpene scores with clustered data points
    sns.violinplot(x='Terpene', y='Score', data=box_plot_df, ax=ax1, color=pastel_colors[0])
    # sns.stripplot(x='Terpene', y='Score', data=box_plot_df, ax=ax1, color=pastel_colors[1], jitter=True, size=5)
    ax1.set_xlabel("Terpene")
    ax1.set_ylabel("Score")
    ax1.set_title("Score Distribution", color=pastel_colors[2])

    # Line plot for terpene frequencies
    sns.lineplot(x=sorted_terps, y=sorted_frequencies, ax=ax2, color=pastel_colors[0])
    for i, freq in enumerate(sorted_frequencies):
        ax2.axhline(freq, color=pastel_colors[0], linestyle="--", alpha=0.1)
    ax2.set_xlabel("Terpene")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Terpene Frequencies", color=pastel_colors[2])

    # Display a list of terpenes with their mean scores
    ax3.bar(sorted_terps, [terp_scores[terp] for terp in sorted_terps], color=pastel_colors[0])
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha="right")
    ax3.set_ylim(min(terp_scores.values()) - 0.1, max(terp_scores.values()) + 0.1)
    ax3.set_ylabel("Mean Score")
    ax3.set_title("Mean Scores", color=pastel_colors[2])

    # Create a network graph of terpene co-occurrences
    G = nx.Graph()
    for strain_terps in terp_co_occurrences:
        G.add_nodes_from(strain_terps)
        for terp1, terp2 in combinations(strain_terps, 2):
            G.add_edge(terp1, terp2)

    pos = nx.spring_layout(G)
    nx.draw(G, pos, ax=ax4, node_color=pastel_colors[1], with_labels=False, edge_color=pastel_colors[0])
    node_labels = nx.draw_networkx_labels(G, pos, font_color=pastel_colors[4], font_size=8, ax=ax4, verticalalignment='bottom')
    for node, (x, y) in pos.items():
        node_labels[node].set_position((x, y + 0.08))  # Adjust the value 0.05 as needed to increase the offset
    ax4.set_title("Terpene Co-occurrence Network", color=pastel_colors[2])

    fig.set_facecolor(bgColor)
    plt.gca().set_facecolor(bgColor)

    # plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    terp_data, terp_co_occurrences, dates = find_terp_data(PAGE_DIR)
    terp_scores = weighted_score_average(terp_data)
    create_graph(terp_data, terp_scores, terp_co_occurrences, dates)
