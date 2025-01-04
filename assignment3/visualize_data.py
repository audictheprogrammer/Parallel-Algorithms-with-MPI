import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Charger les fichiers CSV
file_path1 = "result_bench_intranode.csv"
file_path2 = "result_bench_internode.csv"
file_path3 = "result_bench_memdomain.csv"

# Charger les données
data1 = pd.read_csv(file_path1)
data2 = pd.read_csv(file_path2)
data3 = pd.read_csv(file_path3)

# Fonction pour visualiser un fichier
def visualize_data(data, title, output_image, heatmap_image):
    # Aperçu des données
    print(f"\n{title} - Aperçu des données:\n", data.head())

    # Tracer une courbe lissée pour chaque colonne (hors index)
    plt.figure(figsize=(10, 6))
    for column in data.columns[1:]:  # Supposons que la première colonne est un index (e.g., temps)
        sorted_data = data.sort_values(by=data.columns[0])  # Trier par la première colonne
        plt.plot(sorted_data[data.columns[0]], sorted_data[column], label=column, marker='o')

    plt.xlabel(data.columns[0])
    plt.ylabel("Valeurs")
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.savefig(output_image)
    plt.close()

    # Heatmap avec Seaborn
    plt.figure(figsize=(8, 6))
    sns.heatmap(data.corr(), annot=True, cmap="coolwarm")
    plt.title(f"Correlation Matrix - {title}")
    plt.savefig(heatmap_image)
    plt.close()

# Visualiser les trois fichiers
visualize_data(data1, "Visualisation des données - Intra-node", "image_intranode.png", "heatmap_intranode.png")
visualize_data(data2, "Visualisation des données - Inter-node", "image_internode.png", "heatmap_internode.png")
visualize_data(data3, "Visualisation des données - Memory Domain", "image_memdomain.png", "heatmap_memdomain.png")
