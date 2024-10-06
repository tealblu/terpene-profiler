# Imports and initial setup
from datetime import datetime
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

PAGE_DIR = r"/mnt/storage/Logseq/Home Graph/pages/"
EXCLUDE_FILES = ["templates.md", "README.md", "index.md"]

sns.set_style("dark")
plt.style.use("dark_background")
bgColor = "#383838"
plt.rcParams["axes.facecolor"] = bgColor
plt.rcParams["figure.facecolor"] = bgColor
pastel_colors = {
    "lavender": "#C2AFF0",
    "indigo": "#4747a6",
    "mint": "#d6f8d6",
    "teal": "#519675",
    "orange": "#f79b6a"
}

def extract_terps(line):
    return line.strip().split("::")[1].split(",")

class Terpene:
    def __init__(self, name: str, score: float, position: int, position_limit: int, strain: str = ""):
        self.name = name
        self.score = score
        self.position = position
        self.position_limit = position_limit
        self.strain = strain

    def __repr__(self) -> str:
        return f"{self.name}: {self.score}, weight {self.position}/{self.position_limit}"

strain_terps = {}
strain_scores = {}

def find_terp_data(path: str) -> tuple[list[Terpene], list, list]:
    terp_data = []
    terp_co_occurrences = []
    dates = []

    for filename in os.listdir(path):
        if filename.endswith(".md") and filename not in EXCLUDE_FILES:
            filepath = os.path.join(path, filename)
            with open(filepath, "r", encoding="utf-8") as file:
                for line in file:
                    strain_name = filename[:-3]

                    if line.startswith("date::"):
                        date_str = line.split("::")[1].strip()
                        if date_str == "":
                            continue
                        date_str = date_str.replace("[[", "").replace("]]", "").replace("st", "").replace("nd", "").replace("rd", "").replace("th", "")
                        try:
                            date = datetime.strptime(date_str, "%b %d, %Y")
                        except ValueError:
                            date = datetime.strptime(date_str, "%b %d, %Y %H:%M")

                    if line.startswith("score::"):
                        str_strain_score = line.split("::")[1].strip()
                        if str_strain_score == "":
                            continue
                        strain_score = float(str_strain_score)
                        dates.append((date, strain_score))
                        if strain_name in strain_scores:
                            strain_scores[strain_name].append(strain_score)
                        else:
                            strain_scores[strain_name] = [strain_score]

                    if line.startswith("terps::"):
                        terps = [terp.strip() for terp in extract_terps(line)]
                        for i, terp in enumerate(terps):
                            if terp != "" and terp != "?":
                                terp_data.append(Terpene(terp, strain_score, i + 1, len(terps), strain_name))
                                if strain_name in strain_terps:
                                    strain_terps[strain_name].append(terp)
                                else:
                                    strain_terps[strain_name] = [terp]

    print("\nUnique terpenes:")
    for terp in set(terp.name for terp in terp_data):
        print(terp)

    print("\nData points:")
    for terp in terp_data:
        print(terp)

    return terp_data, terp_co_occurrences, dates

def weighted_score_average(terp_data: list[Terpene]) -> dict[str, float]:
    terp_scores = defaultdict(lambda: [[], [], []])
    for terp in terp_data:
        terp_scores[terp.name][0].append(terp.score)
        terp_scores[terp.name][1].append(terp.position)
        terp_scores[terp.name][2].append(terp.position_limit)

    weighted_avg_scores = {}
    for name, (score_list, position_list, position_limit_list) in terp_scores.items():
        weighted_sum = 0
        total_weight = 0
        for score, position, position_limit in zip(score_list, position_list, position_limit_list):
            weight = (position_limit - position + 1) / position_limit
            weighted_sum += score * weight
            total_weight += weight
        weighted_avg_scores[name] = weighted_sum / total_weight if total_weight else 0

    print("\nWeighted average scores:")
    for name, score in weighted_avg_scores.items():
        print(f"{name}: {score:.2f}")

    return weighted_avg_scores

def create_graph(terp_data: list[Terpene], terp_scores: dict, terp_co_occurrences: list[tuple], dates: list):
    fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), gridspec_kw={"width_ratios": [2, 1]})

    # Step 1: Determine unique terpenes and their frequencies
    unique_terps = set(terp.name for terp in terp_data)
    frequencies = {terp: sum(1 for t in terp_data if t.name == terp) for terp in unique_terps}

    # Step 2: Sort terpenes by their mean scores
    sorted_terps = sorted(unique_terps, key=lambda x: terp_scores[x])
    sorted_frequencies = [frequencies.get(terp, 0) for terp in sorted_terps]

    # Step 3: Prepare box plot data
    box_plot_data = [(terp.name, terp.score) for terp in terp_data]
    box_plot_df = pd.DataFrame(box_plot_data, columns=['Terpene', 'Score'])
    box_plot_df['Mean_Score'] = box_plot_df['Terpene'].map(terp_scores)

    # Ensure that the order is consistent
    box_plot_df['Terpene'] = pd.Categorical(box_plot_df['Terpene'], categories=sorted_terps, ordered=True)
    box_plot_df = box_plot_df.sort_values('Terpene')

    sns.violinplot(x='Terpene', y='Score', data=box_plot_df, ax=ax1, color=pastel_colors["lavender"], bw_adjust=0.5, density_norm='width')
    ax1.set_xlabel("Terpene")
    ax1.set_ylabel("Score")
    ax1.set_title("Score Distribution", color=pastel_colors["mint"])

    sns.barplot(x=sorted_terps, y=sorted_frequencies, ax=ax2, color=pastel_colors["lavender"])
    for i, freq in enumerate(sorted_frequencies):
        ax2.axhline(freq, color=pastel_colors["lavender"], linestyle="--", alpha=0.1)
    ax2.set_xlabel("Terpene")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Terpene Frequencies", color=pastel_colors["mint"])

    ax3.bar(sorted_terps, [terp_scores[terp] for terp in sorted_terps], color=pastel_colors["lavender"])
    ax3.set_xticks(range(len(sorted_terps)))  # Explicitly set the ticks
    ax3.set_xticklabels(sorted_terps, rotation=20, ha="right")
    ax3.set_ylim(min(terp_scores.values()) - 0.1, max(terp_scores.values()) + 0.1)
    ax3.set_ylabel("Mean Score")
    ax3.set_title("Mean Scores", color=pastel_colors["mint"])

    strain_avg_scores = {}
    for strain in strain_scores.keys():
        num_scores = len(strain_scores[strain])
        total_score = sum(strain_scores[strain])  # Simplified calculation
        avg_score = total_score / num_scores
        strain_avg_scores[strain] = avg_score

    sorted_strain_avg_scores = dict(sorted(strain_avg_scores.items(), key=lambda item: item[1], reverse=True))

    raw_data_strings = []
    for strain in sorted_strain_avg_scores.keys():
        if strain in strain_terps:
            raw_data_strings.append(f"{strain} ({strain_avg_scores[strain]}): {' '.join(terp for terp in strain_terps[strain])}")

    ax4.axis('off')
    ax4.set_title("Raw Data Points", color=pastel_colors["mint"])
    for i, data_str in enumerate(raw_data_strings):
        ax4.text(0.01, 0.99 - i * 0.05, data_str, va='top', ha='left', color=pastel_colors["orange"])

    fig.set_facecolor(bgColor)
    plt.gca().set_facecolor(bgColor)
    plt.show()

if __name__ == "__main__":
    terp_data, terp_co_occurrences, dates = find_terp_data(PAGE_DIR)
    terp_scores = weighted_score_average(terp_data)
    create_graph(terp_data, terp_scores, terp_co_occurrences, dates)
