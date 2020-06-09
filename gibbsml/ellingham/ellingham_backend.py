from .fingerprint import Fingerprint
from catlearn.regression.gaussian_process import GaussianProcess
import os


class Ellingham():

    def __init__(self, USER_API_KEY, id_mo, id_m1='', id_m2='',
                 p_O2=0.0, p_CO=0.0, p_CO2=0.0):

        self.user_api_key = USER_API_KEY
        train_fp = Fingerprint(USER_API_KEY=self.user_api_key)

        path_mod = os.path.dirname(__file__)
        train_fp.load_set(path_mod + '/trainingset_ellingham_08June2020.json')

        test_fp = Fingerprint(USER_API_KEY=self.user_api_key)

        # Get training data.
        train_x = train_fp.get_features_values()
        train_y = train_fp.get_target_features_values(
                                        target_feature='dS0_expt'
                                        )
        # Build GP model.
        kernel = [
          {'type': 'gaussian', 'width': 20., 'scaling': 1.5},
          {'type': 'linear', 'scaling': 1., 'constant': 1.},
          ]
        # Train the GP model.
        gp = GaussianProcess(kernel_list=kernel, regularization=1e-2,
                             regularization_bounds=(1e-2, 1e-1),
                             train_fp=train_x, train_target=train_y,
                             optimize_hyperparameters=False,
                             scale_data=True)
        gp.optimize_hyperparameters(global_opt=False,
                                    algomin='TNC',
                                    eval_jac=True)

        # Get test data.
        test_fp.extract_mp_features(id_mo=id_mo,
                                    id_m1=id_m1,
                                    id_m2=id_m2,
                                    selected_features='ellingham')
        test_x = test_fp.get_features_values()

        # Get predictions.
        prediction = gp.predict(test_fp=test_x, uncertainty=False)
        pred = prediction['prediction'][0][0]

        # Store data.
        label = list(test_fp.fp.keys())[0]
        self.balanced_reaction = test_fp.fp[label]['balanced reaction']
        self.dH0 = test_fp.get_feature_value('formation energy (kJ/mol)')[0][0]
        self.dS0 = pred

    def get_dG0(self, T):
        dG0 = self.dH0 + T * self.dS0
        return dG0

    def get_dH0(self):
        return self.dH0

    def get_dS0(self):
        return self.dS0

    def get_balanced_reaction(self):
        return self.balanced_reaction

