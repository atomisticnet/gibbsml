import os
from catlearn.regression.gaussian_process import GaussianProcess

from .fingerprint import Fingerprint

__author__ = "Jose A. Garrido Torres, Alexander Urban"
__email__ = "aurban@atomistic.net"
__date__ = "2021-03-18"
__version__ = "1.0"


class Ellingham():

    def __init__(self, USER_API_KEY, id_mo, id_m1='', id_m2='',
                 id_oxygen='mp-12957', p_O2=0.0, p_CO=0.0, p_CO2=0.0):
        """
        Args:
          USER_API_KEY: Personal API key for the Materials Project (MP)
            database (https://materialsproject.org)
          id_mo: Chemical formula or MP ID of the metal oxide
          id_m1, id_m2: MP IDs of the metal species; will be determined
            automatically if not provided
          id_oxygen: MP ID of the O2 entry to be used as reference
          p_O2, p_CO, p_CO2: O2, CO, and CO2 partial pressures

        """

        self.user_api_key = USER_API_KEY
        train_fp = Fingerprint(USER_API_KEY=self.user_api_key)

        path_mod = os.path.dirname(__file__)
        train_fp.load_set(path_mod + '/trainingset_ellingham_08June2020.json')

        test_fp = Fingerprint(USER_API_KEY=self.user_api_key)

        # Get training data.
        train_x = train_fp.get_features_values()
        train_y = train_fp.get_target_features_values(
            target_feature='dS0_expt')

        # Build GP model.
        kernel = [
            {'type': 'gaussian', 'width': 1., 'scaling': 1.},
            {'type': 'linear', 'scaling': 1., 'constant': 1.}]

        # Train the GP model.
        gp = GaussianProcess(kernel_list=kernel, regularization=1e-3,
                             regularization_bounds=(1e-5, 1e-1),
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
                                    id_oxygen=id_oxygen,
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
