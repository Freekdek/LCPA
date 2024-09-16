import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Get data
veer_table = "data/test/veer_table.xlsx"
df = pd.read_excel(veer_table)

# Plotting
plt.figure(figsize=(10, 6))
sns.scatterplot(x="No. of accessory genes", y="Genome size (Mb)", hue="Group", data=df, s=100)

# Add labels and title
plt.title("Genome size vs. number of accessory genes")
plt.xlabel("Number of accessory genes")
plt.ylabel("Genome size (Mb)")

# Display & save the plot
plt.savefig("data/test/veer_scatterplot.png")
plt.show()

