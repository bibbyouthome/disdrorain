import pandas as pd
import copy as cp
import numpy as np


class LazyProperty(object):
    def __init__(self, func):
        self._func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__

    def __get__(self, obj, klass=None):
        if obj is None:
            return None
        result = obj.__dict__[self.__name__] = self._func(obj)
        return result


class disdrorain_wspeed(object):
    def __init__(self, classpath=None, speedpath=None, datapath=None, dataframe=None, fieldsep=' ',
                 instrument_area=5000, time_interval=60, Cspeed=3.776, bspeed=0.67):
        """
        didrorain class
        # classpath = path to disdrometer class limits file
            if None by default the RD80 Valdvogel class limits are used
        # datapath = path to file
        # dataframe = name of data frame
        # fieldsep = string separating field data
        # instrument_area = catchement area of disrometer in mm^2.
            Default value is 5000 = 50cm^2
        # time_interval = time interval of disdrometer recording in seconds.
            Default value is 60 = 1 minute
        # v(D)=Cspeed*D^bspeed -- terminal velocity of drop of diameter D
        """
        if datapath is not None:
            self.data = pd.read_csv(datapath, sep=fieldsep, header=None)
        if dataframe is not None:
            self.data = dataframe.copy()
 
        #self.data.rename(columns=lambda x: 'C' + str(x + 1), inplace=True)
        self.data.index.name = 'record number'
        if classpath is None:# Thiel
            d = {
                'C1': [0.125, 0.250], 'C2': [0.250, 0.375], 'C3': [0.375, 0.500], 'C4': [0.500, 0.750],
                'C5': [0.750, 1.000], 'C6': [1.000, 1.250], 'C7': [1.250, 1.500], 'C8': [1.500, 1.750],
                'C9': [1.750, 2.000], 'C10': [2.000, 2.500], 'C11': [2.500, 3.000], 'C12': [3.000, 3.500],
                'C13': [3.500, 4.000], 'C14': [4.000, 4.500], 'C15': [4.500, 5.000], 'C16': [5.000, 5.500],
                'C17': [5.500, 6.000], 'C18': [6.000, 6.500], 'C19': [6.500, 7.000], 'C20': [7.000, 7.500], 
                'C21': [7.500, 8.000], 'C22': [8.000, 8.500]}
            self.classlimits = pd.DataFrame(data=d)
            self.classlimits.rename(index={0: 'left', 1: 'right'}, inplace=True)
        else:
            self.classlimits = pd.read_csv(classpath, sep=fieldsep, header=None)
            self.classlimits.rename(columns=lambda x: 'C' + str(x + 1), index={0: 'left', 1: 'right'}, inplace=True)
        if speedpath is None:# Thiel
            v = {
                'S1': 0.1, 'S2': 0.3, 'S3': 0.5, 'S4': 0.7, 'S5': 0.9, 'S6': 1.2, 'S7': 1.6, 'S8': 2.0,
                'S9': 2.4, 'S10': 2.8, 'S11': 3.2, 'S12': 3.8, 'S13': 4.6, 'S14': 5.4, 'S15': 6.2,
                'S16': 7.0, 'S17': 7.8, 'S18': 8.6, 'S19': 9.5, 'S20': 10.5} 
            self.speedvalues = pd.DataFrame(data=v, index=[0])
        else:
            self.speedvalues = pd.read_csv(classpath, sep=fieldsep, header=None)
            self.speedvalues.rename(columns=lambda x: 'S' + str(x + 1), inplace=True)

        colname = []
        diameters_dict = self.classlimits.to_dict()
        speeds_dict = self.speedvalues.to_dict()
        class_eff = {}
        speed_eff = {}
        for i in range(0, self.classlimits.shape[1]):
            for j in range(0, self.speedvalues.shape[1]):
                colname.append(self.classlimits.columns[i] + '_' + self.speedvalues.columns[j])
                class_eff.update({self.classlimits.columns[i] + '_' + self.speedvalues.columns[j]: diameters_dict[f"{self.classlimits.columns[i]}"] })
                speed_eff.update({self.classlimits.columns[i] + '_' + self.speedvalues.columns[j]: speeds_dict[f"{self.speedvalues.columns[j]}"] })
        self.data.columns = colname
        self.classlimits.index.name = 'class borders'
        self.classlimits_effective = pd.DataFrame(data=class_eff)
        self.classlimits_effective.index.name = 'class borders'
        self.speedvalues_effective = pd.DataFrame(data=speed_eff)
        self.instrument_area = instrument_area
        self.time_interval = time_interval
        print()
        

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def bulkvar(self):
        return self.bulk_variables()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def phasespacepar(self):
        return self.phase_space_parameters()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def spectrum(self):
        return self.drop_pdf()

# Moment Calculator method
    def moment_calculator(self, _alpha_):
        """
        - Purpose: calculate the alpha-th moment of the drop size distribution for each record
        - Return: data frame with alpha-th moment as column
        _alpha_: list of moments order (needs to be a python list)
        """

        moments_dict = dict()
        drops_per_record = self.data.sum(axis=1)
        class_span = self.classlimits.loc['right', :] - self.classlimits.loc['left', :]
        for elem in _alpha_:
            class_span_alphap1 = pow(self.classlimits.loc['right', :], elem + 1) - pow(self.classlimits.loc['left', :], elem + 1)
            moments_dict[f"M{elem}"] = self.data.dot(class_span_alphap1 / (class_span * (elem + 1))) / drops_per_record

        _df_ = pd.DataFrame(moments_dict)
        return _df_

    def moment_calculator_cloud(self, _alpha_):
        """
        - Purpose: calculate the alpha-th moment of the drop size distribution for each record
        - Return: data frame with alpha-th moment as column
        _alpha_: list of moments order (needs to be a python list)
        """

        moments_dict = dict()
        drops_per_record = self.data.sum(axis=1)
        class_span = self.classlimits_effective.loc['right', :] - self.classlimits_effective.loc['left', :]
        for elem in _alpha_:
            class_span_alphap1 = pow(self.classlimits_effective.loc['right', :], elem + 1) - pow(self.classlimits_effective.loc['left', :], elem + 1)
            moments_dict[f"M{elem}"] = self.data.dot(class_span_alphap1 / self.speedvalues_effective / (class_span * (elem + 1))) / drops_per_record

        _df_ = pd.DataFrame(moments_dict)
        return _df_

# Remove Outliers method
    def outlier_deletion(self, change_original=False):
        """
        - Purporse: eliminate outliers counts. We look for the class with
        maximun count. From this class, we move to the left (right) until a
        class of null count is found or the minimum (maximum) diameter class
        is reached: this class is the left most (right most) class. All counts
        left of the leftmost and right of the right most class are set to zero.
        E.g. : the following 20 class disdrometer count
        [60 157 121 124 68 67 74 44 14 10 18 11 0 2 0 0 1 0 0 0] is reduced to
        [60 157 121 124 68 67 74 44 14 10 18 11 0 0 0 0 0 0 0 0]
        - Return: an element of the disdrorain class (A or B).
        A) Change input disdrorain class element by eliminating outliers
        B) Create new element of disdrorain class equal ti input except for
        records were outliers are eliminatedself.
        - change_original: flag to decide if disdrorain class element output is
        (A <-> True) or (B <-> False)
        """

# the clean up function finds
# 1) the class with max count
# 2) the class span of consencutive non zero counts which include the class with max count
# All other counts will be considered outliers
        def cleanup(x):
            imax = x.argmax()
            zeroin = np.where(x == 0)[0]
            if len(np.where(zeroin > imax)[0]) > 0:
                rbi = np.where(zeroin > imax)[0].min()
                rb = zeroin[rbi]
                for i in range(rb, len(x)):
                    x[i] = 0
            if len(np.where(zeroin < imax)[0]) > 0:
                lbi = np.where(zeroin < imax)[0].max()
                lb = zeroin[lbi]
                for i in range(0, lb):
                    x[i] = 0
            return x

        # check if we keep original disdrorain class element untouched or not
        if change_original is False:
            self_no_outliers = cp.deepcopy(self)
            rainobj = cp.copy(self_no_outliers)
        else:
            rainobj = cp.copy(self)

        sum_prior = rainobj.data.sum(axis=1)  # total number of drops in each record prior to outlier removal
        _matrix_ = rainobj.data.values
        for i, row in enumerate(_matrix_):
            _matrix_[i] = cleanup(row)
        rainobj.data = pd.DataFrame(_matrix_, columns=rainobj.data.columns, index=rainobj.data.index)
        sum_after = rainobj.data.sum(axis=1)  # total number of drops in each racord after outlier removal

        # create summary data frame -- START
        _df1_ = pd.DataFrame({'ndrops_prior': sum_prior.values})
        _df2_ = pd.DataFrame({'ndrops_after': sum_after.values})
        summary = _df1_.join(_df2_)
        summary.loc[summary.ndrops_prior != summary.ndrops_after, 'delta_drops'] = summary.ndrops_after - summary.ndrops_prior
        summary.loc[summary.ndrops_prior != summary.ndrops_after, 'perc_delta_drops'] = (summary.ndrops_after - summary.ndrops_prior) / summary.ndrops_prior * 100
        summary_final = summary.loc[summary.ndrops_prior != summary.ndrops_after, :]
        summary_final.reset_index(inplace=True)
        summary_final.rename(columns={'index': 'record number'}, inplace=True)
        # create summary data frame -- STOP

        if change_original is False:
            return rainobj, summary_final
        else:
            return summary_final

# Phase Space Parameters method
    def phase_space_parameters(self):
        """
        - Purpose: calculate for each record the drop count (N), mean (mu),
        standard deviation (sigma), skewness (gamma), kurtosis (kappa),
        fifth central moment (eta), sixth central moment (omega)
        - Return: data frame with all the above quantities as columns
        """
        list_moments_order = list([1, 2, 3, 4, 5, 6])
        _df_ = self.moment_calculator(list_moments_order)
        _df_['N'] = self.data.sum(axis=1)
        _df_['mu'] = _df_.M1
        _df_['sigma'] = pow((_df_.M2 - pow(_df_.M1, 2)), 0.5)
        _df_['gamma'] = (_df_.M3 + (2 * pow(_df_.M1, 3)) - (3 * _df_.M1 * _df_.M2)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 1.5))
        _df_['kappa'] = (_df_.M4 - (3 * pow(_df_.M1, 4)) + (6 * _df_.M2 * pow(_df_.M1, 2)) - (4 * _df_.M1 * _df_.M3)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2))
        _df_['eta'] = (_df_.M5 + (4 * pow(_df_.M1, 5)) + (10 * _df_.M3 * pow(_df_.M1, 2)) - (10 * _df_.M2 * pow(_df_.M1, 3)) - (5 * _df_.M4 * _df_.M1)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2.5))
        _df_['omega'] = (_df_.M6 - (6 * _df_.M5 * _df_.M1) + (15 * _df_.M4 * pow(_df_.M1, 2)) - (20 * _df_.M3 * pow(_df_.M1, 3)) + (15 * _df_.M2 * pow(_df_.M1, 4)) - (5 * pow(_df_.M1, 6))) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 3))
        _df_.drop(columns=['M1', 'M2', 'M3', 'M4', 'M5', 'M6'], inplace=True)
        return _df_

# Bulk Variable method
    def bulk_variables(self):
        """
        - Purpose: calculate for each record the Number of drops (N), the number of drops per unit volume N_V,
        rainfall rate (R), reflectivity  (Z), liquid water content (W), rainfall rate per unit drop (R/N=r),
        reflectivity per unit drop (Z/N_V=z), liquid water content per unit drop (W/N_V=w)
        - Return: data frame with all the above quantities as columns
        """

        PI = 3.141592653589793  # approximate value of greek pi
        seconds_in_hour = 3600  # for converting rainfall rate in mm/h
        list_moments_order = list([-self.bspeed, 3 - self.bspeed, 3, 6 - self.bspeed])
        _df_ = self.moment_calculator(list_moments_order)
        # N Nv  R Z W r z w dataframe
        NNvbulk = pd.DataFrame(columns=['N', 'Nv', 'R', 'Z', 'W', 'r', 'z', 'w'])
        NNvbulk = NNvbulk.astype(
            dtype={'N': 'int64', 'Nv': 'float64', 'R': 'float64', 'Z': 'float64', 'W': 'float64', 'r': 'float64',
                   'z': 'float64', 'w': 'float64'})

        # the division by 1000000 is necessary to have catchment area in meter squared
        _df_['N'] = self.data.sum(axis=1)
        _df_['Nv'] = 1 / ((self.instrument_area / 1000000) * self.time_interval * self.Cspeed) * _df_.N * _df_[f"M{-self.bspeed}"]  # this is the concentration per unit volume
        _df_['rainfall_rate'] = (PI / 6) * (1 / (self.instrument_area * self.time_interval)) * seconds_in_hour * _df_.N * _df_.M3
        _df_['reflectivity'] = 1 / ((self.instrument_area / 1000000) * self.time_interval * self.Cspeed) * _df_.N * _df_[f"M{6-self.bspeed}"]
        _df_['lwc'] = 1 / ((self.instrument_area / 1000000) * self.time_interval * self.Cspeed) * _df_.N * _df_[f"M{3-self.bspeed}"]
        _df_.drop(columns=[f"M{-self.bspeed}", f"M{3-self.bspeed}", 'M3', f"M{6-self.bspeed}"], inplace=True)
        return _df_

# instantaneous spectrum method
    def drop_pdf(self):
        """
        - Purpose: calculate the instantaneous drop pdf
        - Return: data frame with pdf value per class (one row per record)
        """
        dropfreq = cp.deepcopy(self.data)
        dropfreq['ndrops'] = dropfreq.sum(axis=1)
        dropfreq = dropfreq.div(dropfreq['ndrops'].values, axis=0)
        dropfreq.drop(['ndrops'], axis=1, inplace=True)
        class_span = self.classlimits.loc['right', :] - self.classlimits.loc['left', :]
        droppdf = cp.deepcopy(dropfreq)

        for row in dropfreq.itertuples():
            for j in range(1, len(row)):
                droppdf.iloc[row[0], j - 1] = row[j] / class_span[j - 1]

        return droppdf
