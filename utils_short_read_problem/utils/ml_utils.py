from pprint import pprint

from sklearn.metrics import silhouette_score
from tslearn.clustering import TimeSeriesKMeans
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
import numpy as np
from tslearn.utils import to_time_series_dataset
import matplotlib.pyplot as plt
import joblib
from tslearn.metrics import cdist_dtw
from sklearn.metrics import silhouette_score
from tslearn.utils import to_time_series_dataset


class MyModel:
    def __init__(self, model=None, labels=None):
        self.model = model
        self.labels = labels
        self.silhouette_scores = None

    def save(self, file_path_name):
        # Save the model
        joblib.dump(self.model, f"{file_path_name}_model.joblib")
        joblib.dump(self.labels, f"{file_path_name}_labels.joblib")

    def load(self, file_path_name):
        self.model = joblib.load(f"{file_path_name}_model.joblib")
        self.labels = joblib.load(f"{file_path_name}_labels.joblib")


def get_number_of_cluster_statistics(data, fastq_name):

    silhouette_scores = []
    # k=2
    # model = TimeSeriesKMeans(n_clusters=k, n_jobs=10, metric="dtw", verbose=True)
    # labels = model.fit_predict(data)

    # Assuming 'labels' are the cluster labels you obtained from TimeSeriesKMeans
    # Calculate the silhouette score using the precomputed distance matrix

    # score = silhouette_score(data, labels, metric="dtw")
    model_name = fastq_name
    model_type = "dtw"
    k_range = [2,10,20,50,100, 300,400,500,600,700,800,900]

    for k in k_range:  # Needs at least 2 clusters to compute silhouette score
        """
        see pipeline_utils.get_clusters_of_sequences()
        """

        try:
            print(f"Loading model..., k={k}")
            my_model = MyModel()
            my_model.load(f"outputs/ML/{model_name}_k_{k}_{model_type}")
            print(f"Loaded., k={k}")
        except:
            print(f"File outputs/ML/{model_name}_k_{k}_{model_type} NOT FOUND. Building model...")
            print(f"Getting model, k={k}")
            model = TimeSeriesKMeans(n_clusters=k, n_jobs=11, metric="dtw", verbose=True)

            print(f"Getting labels, k={k}")
            labels = model.fit_predict(data)

            my_model = MyModel(model=model, labels=labels)
            print(f"Saving model and labels..., k={k}")
            my_model.save(f"outputs/ML/{model_name}_k_{k}_{model_type}")

        # Calculate the DTW distance matrix for all pairs
        print(f"Calculating distnace matrix, k={k}")
        if model_type == "dtw":
            try:
                silhouette_avg = joblib.load(f'outputs/ML/{model_name}_k_{k}_{model_type}_silhouette_avg.joblib')
                print("Loaded shiluet avg")
            except:
                print(f"Calculating distance matrix for silhouette_scores, k={k}")
                distance_matrix = cdist_dtw(data, verbose=5, n_jobs=12)

                print(f"Calculating silhouette_avg, k={k}")
                silhouette_avg = silhouette_score(distance_matrix, my_model.labels, metric='precomputed')
                joblib.dump(silhouette_avg,f'outputs/ML/{model_name}_k_{k}_{model_type}_silhouette_avg.joblib')

            my_silhouette_score = silhouette_avg
        else:
            """
            this one is incomplete
            """
            try:
                my_silhouette_score=joblib.load(
                    f'outputs/ML/{model_name}_k_{k}_{model_type}_silhouette_score.joblib')
                print("Loaded shiluet avg")
            except:
                print(f"Calculating silhouette_avg, k={k}")
                my_silhouette_score = silhouette_score(data, my_model.labels, metric='precomputed')
                joblib.dump(my_silhouette_score,
                            f'outputs/ML/{model_name}_k_{k}_{model_type}_silhouette_score.joblib')

            # joblib.dump(my_silhouette_score, f'outputs/ML/{model_name}_k_{k}_{model_type}_dtw_silhouette_avg.joblib')

        silhouette_scores.append(my_silhouette_score)

    plt.plot(k_range, silhouette_scores, marker='o')
    plt.title('Silhouette Method')
    plt.xlabel('Number of clusters')
    plt.ylabel('Silhouette Score')
    plt.show()
