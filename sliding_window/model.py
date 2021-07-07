import pickle
import numpy as np
from numpy import linalg as LA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler


class ElasticNet():

    def __init__(self, sliding_window, theta=0.05, max_iter=1000, lower=3, upper=7):
        # the cutoff angle at which near-collinear neighboring columns in X are filtered
        self.theta = theta
        # set design matrix and outcome
        self.X = sliding_window.frequencies.drop('subset', 1, inplace=False)
        self.y = sliding_window.frequencies['subset']
        # set hyperparameters ranges
        self.max_iter = max_iter
        self.Cs = self.set_regularization(lower, upper)
        self.l1_ratio = np.arange(start=0.05, stop=1, step=0.1)
        # filepath for saving model output
        self.model_output_path = 'results/model_output/'

    def set_regularization(self, lower, upper):
        self.Cs = 1 / np.power(10, np.arange(lower, upper, step=0.1))

    # FIXME save how many collinear neighbors were filtered for display later
    def filter_collinear_neighbors(self):

        def unit_vector(vector):
            return vector / LA.norm(vector)

        def angle_between(v0, v1):
            return np.arccos(np.clip(np.dot(unit_vector(v0), unit_vector(v1)), -1.0, 1.0))

        mask = [angle_between(self.X[i], self.X[i+1]) > self.theta
                for i in self.X.iloc[:, :-1]]
        mask = [True] + mask
        indices = np.where(mask)
        self.X = self.X.iloc[:, mask]
        return indices

    def scale(self):
        scaler = StandardScaler()
        scaler.fit(self.X)
        self.X = scaler.transform(self.X)

    def fit_elastic_net(self):
        # FIXME indices are currently not used; will be used to find neighbors
        indices = self.filter_collinear_neighbors(self.X)
        # scale data to 0 mean and unit variance
        self.scale(self.X)
        # grid search ranges for regularization parameters
        """ tune hyperparameters for weighted multinomial logistic regression with
            elastic net penalty via stratified k-fold cross validation"""
        elastic_net = LogisticRegressionCV(multi_class='multinomial',
                                           penalty='elasticnet',
                                           class_weight='balanced',
                                           solver='saga',
                                           max_iter=self.max_iter,
                                           Cs=self.Cs,
                                           l1_ratios=self.l1_ratio)
        elastic_net.fit(self.X, self.y)
        with open(f'{self.model_output_path}elastic_net.pickle', 'wb') as f:
            pickle.dump(elastic_net, f)
