import numpy as np

# Example data
data = np.array([0, 0, 1, 1, 2, 2, 3, 10, 11, 11, 11, 12])

# Calculate the 90th percentile as the threshold
threshold = np.percentile(data, 90) / 4

# Filter out the data points that are lower than the threshold.
filtered_data_indices = np.where(data > threshold)[0]
filtered_data_values = data[filtered_data_indices]

print("Threshold based on 90th percentile:", threshold)
print("Data points above threshold:", filtered_data_values)
