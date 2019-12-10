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


class disdrorain(object):
    def __init__(self, classpath=None, datapath=None, dataframe=None, fieldsep=' ',
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
        self.data.rename(columns=lambda x: 'C' + str(x + 1), inplace=True)
        self.data.index.name = 'record number'
        if classpath is None:  # RD80 Valdvogel
            d = {
                'C1': [0.313, 0.405], 'C2': [0.405, 0.505], 'C3': [0.505, 0.596], 'C4': [0.596, 0.715],
                'C5': [0.715, 0.827], 'C6': [0.827, 0.999], 'C7': [0.999, 1.232], 'C8': [1.232, 1.429],
                'C9': [1.429, 1.582], 'C10': [1.582, 1.748], 'C11': [1.748, 2.077], 'C12': [2.077, 2.441],
                'C13': [2.441, 2.727], 'C14': [2.727, 3.011], 'C15': [3.011, 3.385], 'C16': [3.385, 3.704],
                'C17': [3.704, 4.127], 'C18': [4.127, 4.573], 'C19': [4.573, 5.145], 'C20': [5.145, 5.601]}
            self.classlimits = pd.DataFrame(data=d)
            self.classlimits.rename(index={0: 'left', 1: 'right'}, inplace=True)
        else:
            self.classlimits = pd.read_csv(classpath, sep=fieldsep, header=None)
            self.classlimits.rename(columns=lambda x: 'C' + str(x + 1), index={0: 'left', 1: 'right'}, inplace=True)
        self.classlimits.index.name = 'class borders'
        self.instrument_area = instrument_area
        self.time_interval = time_interval
        self.Cspeed = Cspeed
        self.bspeed = bspeed

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
