import pickle
import numpy as np
from numpy import linalg as LA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler


class ElasticNet():

    def __init__(self, sliding_window, theta=0.1, max_iter=1000, upper=1e-3, lower=1e-7, n_cv_values=20):
        # FIXME upper must be larger than lower, make a check
        # the cutoff angle at which near-collinear neighboring columns in X are filtered
        self.theta = theta
        # set design matrix and outcome
        self.y = sliding_window.frequencies['subset']
        self.X = sliding_window.frequencies.drop(
            ['subset', 'header'], axis=1, inplace=False)
        # set hyperparameters ranges
        self.max_iter = max_iter
        self.n_cv_values = n_cv_values
        self.Cs = self.make_cv_grid(upper, lower)
        self.l1_ratio = np.arange(start=0.05, stop=1, step=0.1)
        # filepath for saving model output
        self.model_output_path = 'results/model_output/'

    def make_cv_grid(self, upper, lower):

        def get_exponent(number):
            base10 = np.log10(abs(number))
            return abs(np.floor(base10))

        lower = get_exponent(lower)
        upper = get_exponent(upper)
        step = abs(upper - lower) / self.n_cv_values
        return 1 / np.power(10, np.arange(start=upper, stop=lower, step=step))

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

    def scale_X(self):
        scaler = StandardScaler()
        scaler.fit(self.X)
        self.X = scaler.transform(self.X)

    def fit_elastic_net(self):
        # FIXME indices are currently not used; will be used to find neighbors
        indices = self.filter_collinear_neighbors()
        # scale data to 0 mean and unit variance
        self.scale_X()
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
