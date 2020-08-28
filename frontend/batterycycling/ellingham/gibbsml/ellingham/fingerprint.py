from pymatgen.ext.matproj import MPRester
from ase.io import read
import os
import pandas as pd
from mendeleev import element
import numpy as np
from functools import reduce
from math import gcd
import json


class Fingerprint:

    def __init__(self, USER_API_KEY):
        """
        Parameters
        ----------
        USER_API_KEY: str
            User's API key generated from the Materials Project database.
            See https://materialsproject.org/open
        """

        self.user_api_key = USER_API_KEY
        self.fp = {}  # Dictionary containing the set of features.
        self.implemented_features = ['stoichiometry', 'electronegativity',
                                     'mass', 'volume', 'density',
                                     'bulk modulus', 'shear modulus',
                                     'poisson ratio', 'anisotropy',
                                     'spacegroup', 'ionic character']
        self.selected_features = None
        self.label = None

    def extract_mp_features(self, id_mo, id_m1='', id_m2='',
                            id_oxygen='mp-12957',
                            selected_features='all', label=None,
                            mo_energy_correction=True):
        """
        Generates a feature set for an oxidation for a given metal
        oxide (AxByOz) from the elements (A and B).

        Parameters
        ----------
        id_mo: str
            Materials Project mp-id for the metal oxide or chemical formula
            of the metal oxide, e.g. 'Al2SiO5' or 'mp-4753'.
        id_m1: str
            (optional) Materials Project mp-id for the metal A, e.g. 'mp-134'.
        id_m2: str
            (optional) Materials Project mp-id for the metal B, e.g. 'mp-149'.
        id_oxygen: str
            Materials project mp-id for oxygen in the gas phase.
        selected_features: list or str
            (option 1): list
                List of selected features to be considered to
                generate the fingerprint. Implemented are: 'stoichiometry',
                'electronegativity', 'mass', 'volume', 'density',
                'bulk modulus', 'shear modulus', 'poisson ratio',
                'anisotropy', 'spacegroup' and 'ionic character'.
            (option 2): str
                'all': Include all implemented features (see option 1).
                'ellingham': Recommended features for building models for
                             predicting Ellingham diagrams. Includes only
                             the following features:
                             'stoichiometry', 'electronegativity',
                             'density', 'bulk modulus', 'ionic character'.

        'label': str
            Defines the label tag for the fingerprint. The user can chose a
            name for the fingerprint for the data entry, e.g. 'Al2SiO5-PBEU'.
        'mo_energy_correction': bool
            If True the algorithm only selects Material Project entries which
            in which energy corrections are available. See:
            https://materialsproject.org/docs/calculations#Total_Energy_Adjustments

        """

        # Implemented features.
        if selected_features == 'all':
            self.selected_features = self.implemented_features
        if selected_features == 'ellingham':
            self.selected_features = ['stoichiometry', 'electronegativity',
                                      'density', 'bulk modulus',
                                      'ionic character']
        else:
            self.selected_features = selected_features

        # Get ID for the oxidized state.
        print("Getting information for " + id_mo)
        if "mp" not in id_mo:
            id_mo = self._find_id(id_mo)
        with MPRester(self.user_api_key) as m:
            data_mo = m.get_data(id_mo)[0]

        # Set oxygen energy and MO energy corrections.
        with MPRester(self.user_api_key) as m:
            e_adjus = m.get_entries(id_mo)[0].__dict__['energy_adjustments']
            e_mo_corr = 0
            for e_ad in e_adjus:
                e_mo_corr += e_ad.value
        if mo_energy_correction is False:
            e_mo_corr = 0.0
        data_o2 = m.get_data(id_oxygen)[0]
        e_o2 = 2 * data_o2['energy'] / data_o2['unit_cell_formula']['O']

        # Recognize whether is a unary or binary oxide.
        n_elements = data_mo['nelements']
        binary_oxide = False
        if n_elements == 3:
            binary_oxide = True
        msg = "Only unary and binary oxides are implemented."
        assert n_elements <= 3, NotImplementedError(msg)

        elements = data_mo['elements']
        element_no_ox = np.array(elements)[~np.isin(elements, 'O')]

        if binary_oxide is True:
            element_m1, element_m2 = element_no_ox[0], element_no_ox[1]
        else:
            element_m1, element_m2 = element_no_ox[0], element_no_ox[0]

        # Get info for M1 and M2.
        if "mp" not in id_m1:
            id_m1 = self._find_id(element_m1)
        if "mp" not in id_m2:
            id_m2 = self._find_id(element_m2)
        data_m1 = m.get_data(id_m1)[0]
        data_m2 = m.get_data(id_m2)[0]

        # Get formula for the metals and oxides.
        formula_mo = data_mo['pretty_formula']
        formula_m1 = data_m1['pretty_formula']
        formula_m2 = data_m2['pretty_formula']

        # Set a label for the compound if not specified.
        self.label = label
        if self.label is None:
            self.label = formula_mo + '_' + id_mo

        # Initialize dict with features and training features.
        self.fp.update({self.label: {}})  # A dictionary for each entry.
        self.fp[self.label]['target_features'] = {}
        self.fp[self.label]['features'] = {}

        # Read files and convert to ASE Atoms objects.
        atoms = []
        for i in ['m1', 'm2', 'mo']:
            f_atoms = open('tmp_Atoms.cif', 'w')
            f_atoms.write(eval("data_" + i)['cif'])
            f_atoms.close()
            atoms.append(read('tmp_Atoms.cif'))
            os.remove('tmp_Atoms.cif')

        atoms_m1, atoms_m2, atoms_mo = atoms

        # Get formula unit for M1, M2 and MO.
        n_m1, fu_m1 = 2 * (len(atoms_m1),)
        n_m2, fu_m2 = 2 * (len(atoms_m2),)
        n_mo, fu_mo = len(atoms_mo), self._get_atoms_per_unit_formula(atoms_mo)
        n_m1_in_mo = self._get_number_of_atoms_element(atoms_mo,
                                                       symbol=element_m1)
        n_m1_in_mo /= fu_mo
        n_m2_in_mo = self._get_number_of_atoms_element(atoms_mo,
                                                       symbol=element_m2)
        n_m2_in_mo /= fu_mo
        n_ox_in_mo = self._get_number_of_atoms_element(atoms_mo, symbol='O')
        n_ox_in_mo /= fu_mo

        # Update dictionary with extra info.
        self.fp[self.label]['raw data'] = {}
        self.fp[self.label]['raw data']['data mo'] = data_mo
        self.fp[self.label]['raw data']['data m1'] = data_m1
        self.fp[self.label]['raw data']['data m2'] = data_m2

        # Formation energy.
        # Get formation energies (a M1 + b M2 + O2 --> c M1xM2yOz).
        e_m1, e_m2 = data_m1['energy'], data_m2['energy']
        e_mo = data_mo['energy']
        x, y, z = n_m1_in_mo, n_m2_in_mo, n_ox_in_mo
        a = (2 / z) * x
        b = (2 / z) * y
        c = 2 / z
        if not binary_oxide:
            a /= 2
            b /= 2
        dH0 = c * (e_mo + e_mo_corr)/fu_mo
        dH0 -= a * e_m1/fu_m1
        dH0 -= b * e_m2/fu_m2
        dH0 -= e_o2
        dH0 *= 96.485  # eV to kJ/mol.
        self.add_feature(description='formation energy (kJ/mol)', value=dH0)

        balanced_reaction = None

        if not binary_oxide:
            balanced_reaction = str(2 * a) + " " + formula_m1 + " + O2"
            balanced_reaction += " --> "
            balanced_reaction += str(c) + " " + formula_mo

        if binary_oxide:
            balanced_reaction = str(a) + " " + formula_m1 + " + "
            balanced_reaction += str(b) + " " + formula_m2 + " + O2"
            balanced_reaction += " --> "
            balanced_reaction += str(c) + " " + formula_mo

        self.fp[self.label]['balanced reaction'] = balanced_reaction

        # Stoichiometry.
        if 'stoichiometry' in self.selected_features:

            # Ratio metal / oxygen:
            ratio_o_m1 = n_m1_in_mo / n_ox_in_mo
            ratio_o_m2 = n_m2_in_mo / n_ox_in_mo
            # Average oxidation state of the metals:
            av_ox_state = (n_ox_in_mo * 2) / (n_m1_in_mo + n_m2_in_mo)
            # Atomic number:
            element_m1, element_m2 = atoms_m1[0].symbol, atoms_m2[0].symbol
            z_m1 = element(element_m1).atomic_number
            z_m2 = element(element_m2).atomic_number

            self.add_feature(description='ratio metal oxygen (mean)',
                             value=np.mean([ratio_o_m1, ratio_o_m2]))
            self.add_feature(description='ratio metal oxygen (var)',
                             value=np.var([ratio_o_m1, ratio_o_m2]))
            self.add_feature(description='average oxidation state',
                             value=av_ox_state)
            self.add_feature(description='atomic number (mean)',
                             value=np.mean([z_m1, z_m2]))
            self.add_feature(description='atomic number (var)',
                             value=np.var([z_m1, z_m2]))

        # Electronegativity (Pauling).
        if 'electronegativity' in self.selected_features:
            elecneg_m1 = element(element_m1).en_pauling
            elecneg_m2 = element(element_m2).en_pauling
            self.add_feature(description='pauling electronegativity (mean)',
                              value=np.mean([elecneg_m1, elecneg_m2]))
            self.add_feature(description='pauling electronegativity (var)',
                              value=np.var([elecneg_m1, elecneg_m2]))

        # Percentage of ionic character (Pauling).
        if 'ionic character' in self.selected_features:
            elnegdif_m1 = element('O').en_pauling - element(element_m1).en_pauling
            elnegdif_m2 = element('O').en_pauling - element(element_m2).en_pauling
            pio_m1 = 100 * (1 - np.exp(-(1/2 * elnegdif_m1)**2))
            pio_m2 = 100 * (1 - np.exp(-(1/2 * elnegdif_m2)**2))
            self.add_feature(description='% ionic character (mean)',
                              value=np.mean([pio_m1, pio_m2]))
            self.add_feature(description='% ionic character (var)',
                              value=np.var([pio_m1, pio_m2]))

        # Volume.
        if 'volume' in self.selected_features:
            V_m1 = atoms_m1.get_volume()
            V_per_fu_m1 = V_m1 / fu_m1
            V_m2 = atoms_m2.get_volume()
            V_per_fu_m2 = V_m2 / fu_m2
            V_mo = atoms_mo.get_volume()
            V_per_fu_mo = V_mo / fu_mo
            self.add_feature(description='volume per formula unit (mean)',
                              value=np.mean([V_per_fu_m1, V_per_fu_m2]))
            self.add_feature(description='volume per formula unit (var)',
                              value=np.var([V_per_fu_m1, V_per_fu_m2]))
            self.add_feature(description='volume MO per formula unit',
                              value=V_per_fu_mo)
            diff_V_per_fu_m1_mo = V_per_fu_mo - V_per_fu_m1
            diff_V_per_fu_m2_mo = V_per_fu_mo - V_per_fu_m2
            self.add_feature(description='difference volume (MO-M) (mean)',
                              value=np.mean([diff_V_per_fu_m1_mo, diff_V_per_fu_m2_mo]))
            self.add_feature(description='difference volume (MO-M) (var)',
                              value=np.var([diff_V_per_fu_m1_mo, diff_V_per_fu_m2_mo]))

        # Mass.
        if 'mass' in self.selected_features:
            mass_m1 = np.average(atoms_m1.get_masses())
            mass_m2 = np.average(atoms_m2.get_masses())
            mass_mo = np.average(atoms_mo.get_masses())
            mass_per_fu_m1 = mass_m1 / fu_m1
            mass_per_fu_m2 = mass_m2 / fu_m2
            mass_per_fu_mo = mass_mo / fu_mo
            self.add_feature(description='mass per formula unit (mean)',
                              value=np.mean([mass_per_fu_m1, mass_per_fu_m2]))
            self.add_feature(description='mass per formula unit (var)',
                              value=np.var([mass_per_fu_m1, mass_per_fu_m2]))
            self.add_feature(description='mass MO per formula unit',
                              value=mass_per_fu_mo)
            diff_mass_per_fu_m1_mo = mass_per_fu_mo - mass_per_fu_m1
            diff_mass_per_fu_m2_mo = mass_per_fu_mo - mass_per_fu_m2
            self.add_feature(description='difference mass (MO-M) (mean)',
                              value=np.mean([diff_mass_per_fu_m1_mo, diff_mass_per_fu_m2_mo]))
            self.add_feature(description='difference mass (MO-M) (var)',
                              value=np.var([diff_mass_per_fu_m1_mo, diff_mass_per_fu_m2_mo]))

        # Density.
        if 'density' in self.selected_features:
            dens_m1 = data_m1['density']
            dens_m2 = data_m2['density']
            dens_mo = data_mo['density']
            self.add_feature(description='density (mean)',
                              value=np.mean([dens_m1, dens_m2]))
            self.add_feature(description='density (var)',
                              value=np.var([dens_m1, dens_m2]))
            self.add_feature(description='density MO', value=dens_mo)
            diff_dens_m1_mo = dens_mo - dens_m1
            diff_dens_m2_mo = dens_mo - dens_m2
            self.add_feature(description='difference density (MO-M) (mean)',
                              value=np.mean([diff_dens_m1_mo, diff_dens_m2_mo]))
            self.add_feature(description='difference density (MO-M) (var)',
                              value=np.var([diff_dens_m1_mo, diff_dens_m2_mo]))

        # Bulk modulus.
        if 'bulk modulus' in self.selected_features:
            elas_m1, elas_m2, elas_mo = data_m1['elasticity'], data_m2['elasticity'], data_mo['elasticity']
            Kv_m1, Kv_m2, Kv_mo = elas_m1['K_Voigt'], elas_m2['K_Voigt'], elas_mo['K_Voigt']
            self.add_feature(description='bulk modulus (mean)',
                              value=np.mean([Kv_m1, Kv_m2]))
            self.add_feature(description='bulk modulus (var)',
                              value=np.var([Kv_m1, Kv_m2]))
            self.add_feature(description='bulk modulus MO', value=Kv_mo)
            diff_Kv_m1_mo = Kv_mo - Kv_m1
            diff_Kv_m2_mo = Kv_mo - Kv_m2
            self.add_feature(description='difference bulk modulus (MO-M) (mean)',
                       value=np.mean([diff_Kv_m1_mo, diff_Kv_m2_mo]))
            self.add_feature(description='difference bulk modulus (MO-M) (var)',
                       value=np.var([diff_Kv_m1_mo, diff_Kv_m2_mo]))

        # Shear modulus.
        if 'shear modulus' in self.selected_features:
            elas_m1, elas_m2, elas_mo = data_m1['elasticity'], data_m2['elasticity'], data_mo['elasticity']
            Gv_m1, Gv_m2, Gv_mo = elas_m1['G_Voigt'], elas_m2['G_Voigt'], elas_mo['G_Voigt']
            self.add_feature(description='shear modulus (mean)',
                              value=np.mean([Gv_m1, Gv_m2]))
            self.add_feature(description='shear modulus (var)',
                              value=np.var([Gv_m1, Gv_m2]))
            self.add_feature(description='shear modulus MO', value=Gv_mo)
            diff_Gv_m1_mo = Gv_mo - Gv_m1
            diff_Gv_m2_mo = Gv_mo - Gv_m2
            self.add_feature(description='difference shear modulus (MO-M) (mean)',
                              value=np.mean([diff_Gv_m1_mo, diff_Gv_m2_mo]))
            self.add_feature(description='difference shear modulus (MO-M) (var)',
                              value=np.var([diff_Gv_m1_mo, diff_Gv_m2_mo]))

        # Poissons Ratio.
        if 'poisson ratio' in self.selected_features:
            elas_m1, elas_m2, elas_mo = data_m1['elasticity'], data_m2['elasticity'], data_mo['elasticity']
            pois_m1, pois_m2, pois_mo = elas_m1['poisson_ratio'], elas_m2['poisson_ratio'], elas_mo['poisson_ratio']
            self.add_feature(description='poisson ratio (mean)',
                              value=np.mean([pois_m1, pois_m2]))
            self.add_feature(description='poisson ratio (var)',
                              value=np.var([pois_m1, pois_m2]))
            self.add_feature(description='poisson ratio MO', value=pois_mo)
            diff_pois_m1_mo = pois_mo - pois_m1
            diff_pois_m2_mo = pois_mo - pois_m2
            self.add_feature(description='difference poisson ratio (MO-M) (mean)',
                              value=np.mean([diff_pois_m1_mo, diff_pois_m2_mo]))
            self.add_feature(description='difference poisson ratio (MO-M) (var)',
                              value=np.var([diff_pois_m1_mo, diff_pois_m2_mo]))

        # Universal anisotropy.
        if 'anisotropy' in self.selected_features:
            elas_m1, elas_m2, elas_mo = data_m1['elasticity'], data_m2['elasticity'], data_mo['elasticity']
            u_ani_m1, u_ani_m2, u_ani_mo = elas_m1['universal_anisotropy'], elas_m2['universal_anisotropy'], elas_mo['universal_anisotropy']
            el_ani_m1, el_ani_m2, el_ani_mo = elas_m1['elastic_anisotropy'], elas_m2['elastic_anisotropy'], elas_mo['elastic_anisotropy']
            self.add_feature(description='universal anisotropy (mean)',
                              value=np.mean([u_ani_m1, u_ani_m2]))
            self.add_feature(description='universal anisotropy (var)',
                              value=np.var([u_ani_m1, u_ani_m2]))
            self.add_feature(description='universal anisotropy MO',
                              value=u_ani_mo)
            diff_u_ani_m1_mo = u_ani_mo - u_ani_m1
            diff_u_ani_m2_mo = u_ani_mo - u_ani_m2
            self.add_feature(description='difference universal anisotropy (MO-M) (mean)',
                              value=np.mean([diff_u_ani_m1_mo, diff_u_ani_m2_mo]))
            self.add_feature(description='difference universal anisotropy (MO-M) (var)',
                              value=np.var([diff_u_ani_m1_mo, diff_u_ani_m2_mo]))

            # Elastic anisotropy.
            self.add_feature(description='elastic anisotropy (mean)',
                              value=np.mean([el_ani_m1, el_ani_m2]))
            self.add_feature(description='elastic anisotropy (var)',
                              value=np.var([el_ani_m1, el_ani_m2]))
            self.add_feature(description='elastic anisotropy MO',
            value=el_ani_mo)
            diff_el_ani_m1_mo = el_ani_mo - el_ani_m1
            diff_el_ani_m2_mo = el_ani_mo - el_ani_m2
            self.add_feature(description='difference elastic anisotropy (MO-M) (mean)',
                              value=np.mean([diff_el_ani_m1_mo, diff_el_ani_m2_mo]))
            self.add_feature(description='difference elastic anisotropy (MO-M) (var)',
                              value=np.var([diff_el_ani_m1_mo, diff_el_ani_m2_mo]))

        # Spacegroup.
        if 'spacegroup' in self.selected_features:
            spacegroup_m1 = data_m1['spacegroup']['number']
            spacegroup_m2 = data_m2['spacegroup']['number']
            spacegroup_mo = data_mo['spacegroup']['number']
            self.add_feature(description='Spacegroup M1', value=spacegroup_m1)
            self.add_feature(description='Spacegroup M2', value=spacegroup_m2)
            self.add_feature(description='Spacegroup MO', value=spacegroup_mo)

        print("Fingerprint for " + self.label + " completed.")

    def _find_id(self, compound):
        """ Find Materials Project ID for a given compound.

            Parameters
            ----------
            compound: str
                Compound formula.  Examples: ``'Li2O'``,
                ``'AlLiO2'``, ``'Mg2SiO4'``.
            user_api: str
                Materials Project users API.

            Returns
            -------
            id_compound: str
                Materials Project compound ID.

        """

        with MPRester(self.user_api_key) as m:
            elasticity = None
            i = 0
            while elasticity is None:
                info_MO = m.get_entries(compound, inc_structure='final',
                                        property_data=['elasticity',
                                                       'e_above_hull',
                                                       'Correction'],
                                        sort_by_e_above_hull=True)[i]
                id_compound = info_MO.__dict__['entry_id']
                elasticity = m.get_data(id_compound)[0]['elasticity']
                i += 1
        return id_compound

    def add_feature(self, description, value):
        """
        Parameters
        ----------
        description: str
            Description of the property to append (e.g. 'Formation energy').
        value: float
            Numerical value for a given property.

        Returns
        -------
        Adds a feature to the fingerprint (stored in self.fp).
        """

        self.fp[self.label]['features'].update({description: value})

    def add_target_feature(self, description, value):
        """
        Parameters
        ----------
        description: str
            Description of the property to be appended (e.g. 'Reduction
            temperature').
        value: float
            Numerical value for a given property.

        Returns
        -------
        Adds a target feature to the fingerprint (stored in self.fp).
        Note: In this case the properties and values will be only used for
        training the model (commonly know as train_y).
        """
        self.fp[self.label]['target_features'].update({description: value})

    def _get_number_of_atoms_element(self, atoms, symbol='O'):
        n_atom = 0
        for atom in atoms:
            if atom.symbol == symbol:
                n_atom += 1
        return n_atom

    def _get_atoms_per_unit_formula(self, atoms):
        symbols = []
        for atom in atoms:
            symbols.append(atom.symbol)
        a = np.unique(symbols, return_counts=True)[1]
        return reduce(gcd, a)

    def _get_number_of_species(self, atoms):
        symbols = []
        for atom in atoms:
            symbols.append(atom.symbol)
        a = np.unique(symbols)
        return len(a), a

    def get_feature_value(self, feature_name):
        species = list(self.fp.keys())
        feature_value = []
        for i in species:
            feature = self.fp[i]['features'][feature_name]
            feature_value.append([feature])
        return feature_value

    def get_features_values(self):
        species = list(self.fp.keys())
        features_names = list(self.fp[species[0]]['features'].keys())
        features_values = []
        for i in species:
            features_i = []
            for j in features_names:
                features_i.append(self.fp[i]['features'][j])
            features_values.append(features_i)
        val_shape = np.shape(features_values)
        features_values = np.reshape(features_values,
                                     (val_shape[0], val_shape[1])
                                     )
        return features_values

    def get_target_features_values(self, target_feature='dS0_expt'):
        species = list(self.fp.keys())
        target_features_values = []
        for i in species:
            feature = self.fp[i]['target_features'][target_feature]
            target_features_values.append([feature])
        val_shape = np.shape(target_features_values)
        target_features_values = np.reshape(target_features_values,
                                            (val_shape[0], -1)
                                            )
        return target_features_values

    def get_labels(self):
        """
        Returns the list of species (labelled), e.g. CaO-mp-2605.
        """
        return list(self.fp.keys())

    def get_features_names(self):
        """
        Returns a list containing the names of the features, e.g. formation
        energy (kJ/mol).
        """
        species = list(self.fp.keys())
        features_names = list(self.fp[species[0]]['features'].keys())
        return features_names

    def get_target_features_names(self):
        """
        Returns the list of target features. These features are
        user-defined and must be included with the add_target_features
        function.
        """
        species = list(self.fp.keys())
        features_names = list(self.fp[species[0]]['target_features'].keys())
        return features_names

    def dump_set(self, filename='fingerprint.json'):
        """
        Parameters
        ----------
        filename: str
            Name of the file to save the generated Fingerprint class.

        Returns
        -------
        Saves the whole Fingerprint class into a json file.
        """
        self.label = 0.0  # Remove label to prevent confusion.
        fp_dict = self.__dict__
        del fp_dict['user_api_key']  # Important! Don't store API key.
        with open(filename, 'w') as fp:
            json.dump(fp_dict, fp)

    def load_set(self, filename):
        """
        Parameters
        ----------
        filename: str
            Name of the file to load (json format).

        Returns
        -------
        Load a json file containing a previously saved Fingerprint (see
        dump_set function).
        """
        with open(filename, 'r') as fp:
            self.__dict__ = json.load(fp)
