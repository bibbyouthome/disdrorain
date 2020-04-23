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
                instrument_area=5000, time_interval=60, 
                Aplawspeed=3.776, Bplawspeed=0.67,
                Aexpospeed=9.65, Bexpospeed=10.3, Cexpospeed=0.6,
                ncells=25):
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
        # v(D)=Aplawspeed*D^Bplawspeed -- terminal velocity of drop of diameter D (power law velocity)
            see article:
        # v(D)=Aexpospeed-Bexpospeed*exp(-Cexpospeed*D) -- terminal velocity of drop of diameter D ("exponenetial" law velocity)
            see article:
        """

        # default diameter classes: RD80 Valdvogel
        default_d = {'C1': [0.313, 0.405], 'C2': [0.405, 0.505], 'C3': [0.505, 0.596], 'C4': [0.596, 0.715],
            'C5': [0.715, 0.827], 'C6': [0.827, 0.999], 'C7': [0.999, 1.232], 'C8': [1.232, 1.429],
            'C9': [1.429, 1.582], 'C10': [1.582, 1.748], 'C11': [1.748, 2.077], 'C12': [2.077, 2.441],
            'C13': [2.441, 2.727], 'C14': [2.727, 3.011], 'C15': [3.011, 3.385], 'C16': [3.385, 3.704],
            'C17': [3.704, 4.127], 'C18': [4.127, 4.573], 'C19': [4.573, 5.145], 'C20': [5.145, 5.601]}
        if datapath is not None:
            self.data = pd.read_csv(datapath, sep=fieldsep, header=None)
        if dataframe is not None:
            self.data = dataframe.copy()
        self.data.rename(columns=lambda x: 'C' + str(x + 1), inplace=True)
        self.data.index.name = 'record number'
        if classpath is None: # use default diameter classes
            self.classlimits = pd.DataFrame(data=default_d)
            self.classlimits.rename(index={0: 'left', 1: 'right'}, inplace=True)
        else:
            self.classlimits = pd.read_csv(classpath, sep=fieldsep, header=None)
            self.classlimits.rename(columns=lambda x: 'C' + str(x + 1), index={0: 'left', 1: 'right'}, inplace=True)
        self.classlimits.index.name = 'class borders'
        self.instrument_area = instrument_area
        self.time_interval = time_interval
        self.Aplawspeed = Aplawspeed
        self.Bplawspeed = Bplawspeed
        self.Aexpospeed = Aexpospeed
        self.Bexpospeed = Bexpospeed
        self.Cexpospeed = Cexpospeed
        self.ncells = ncells

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def bulkvar_vplaw(self):
        return self.bulk_variables()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def bulkvar_vexpo(self):
        return self.bulk_variables(_speed_='expo')

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def psp_flux(self):
        return self.flux_phase_space_parameters()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def psp_cloud_vplaw(self):
        return self.cloud_phase_space_parameters()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def psp_cloud_vexpo(self):
        return self.cloud_phase_space_parameters(_speed_='expo')

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    def spectrum_flux(self):
        return self.flux_drop_pdf()

# Moment Calculator method
    def flux_moment_calculator(self, _alpha_):
        """
        - Purpose: calculate the alpha-th moment of the "flux" drop size distribution for each record
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

    def cloud_moment_calculator( self,_alpha_,_speed_='plaw'):
        """
        - Purpose: calculate the alpha-th moment of the "cloud" drop size distribution for each record
        - Return: data frame with alpha-th moment as column
        _alpha_: list of moments order (needs to be a python list)

        """

        zero_drop = False
        # if _speed_ is anything but 'expo' the plaw formula for drop velocity is considered
        if (_speed_=='expo'):
            if 0 not in _alpha_:
                _alpha_.append(0)
                zero_drop = True
            moments_dict = dict()
            drops_per_record = self.data.sum(axis=1)
            class_df = self.classlimits
            _array_=np.zeros(class_df.shape[1])
            for elem in _alpha_:
                for diamclass in range(0,class_df.shape[1]):
                    class_span = round(class_df.loc['right',class_df.columns[diamclass]]-class_df.loc['left',class_df.columns[diamclass]],3)
                    _delta_ = class_span/self.ncells
                    _sum_ = 0.0
                    for j in range (0,self.ncells):
                        d_effective = class_df.loc['left',class_df.columns[diamclass]] + (_delta_/2) + j*_delta_
                        _sum_ = _sum_ + pow(d_effective,elem) / (self.Aexpospeed-self.Bexpospeed*np.exp(-self.Cexpospeed*d_effective))*_delta_
                    _array_[diamclass] = _sum_ / class_span
                _series_ = pd.Series(_array_,index=class_df.columns)
                
                moments_dict[f"X{elem}"] = self.data.dot(_series_) / drops_per_record

            #
            _df_ = pd.DataFrame(moments_dict)
            columns_drop = list()
            for elem in _alpha_:
                columns_drop.append(f"X{elem}")
    
            #
            if zero_drop == True:
                columns_drop.append('M0')
                for elem in _alpha_:
                    _df_[f"M{elem}"] = _df_[f"X{elem}"]/_df_.X0
            else:
                for elem in _alpha_:
                    if elem != 0:
                        _df_[f"M{elem}"] = _df_[f"X{elem}"]/_df_.X0
                    else:
                        _df_[f"M{elem}"] = _df_[f"X{elem}"]

            _df_.drop(columns=columns_drop,inplace=True) 

        # v(D) is power law
        else:
            if 0 not in _alpha_:
                _alpha_.append(0)
                zero_drop = True
            _alpha_new = _alpha_.copy()
            i=0
            for elem in _alpha_: 
                _alpha_new[i] = _alpha_[i]-self.Bplawspeed
                i=i+1
            _df_ = self.flux_moment_calculator(_alpha_new)/self.Aplawspeed

            #
            columns_drop = list()
            for elem in _alpha_new:
                columns_drop.append(f"M{elem}")
            
            #
            if zero_drop == True:
                columns_drop.append('X-0.67')
                for elem in _alpha_new:
                    _df_[f"X{elem}"] = _df_[f"M{elem}"]/_df_['M-0.67']
            else:
                for elem in _alpha_new:
                    if elem != -0.67:
                        _df_[f"X{elem}"] = _df_[f"M{elem}"]/_df_['M-0.67']
                    else:
                        _df_[f"X{elem}"] = _df_[f"M{elem}"]

            _df_.drop(columns=columns_drop,inplace=True) 
            
            #
            columns_name=list()
            if zero_drop == True:
                for elem in _alpha_:
                    if elem != 0:
                        columns_name.append(f"M{elem}")
            else:
                for elem in _alpha_:
                    columns_name.append(f"M{elem}")
            
            _df_.columns=columns_name 
        
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
        B) Create new element of disdrorain class equal to input except for
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
            # self_no_outliers = cp.deepcopy(self)
            rainobj = cp.deepcopy(self)
        else:
            rainobj = self

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
        summary_final = cp.deepcopy(summary.loc[summary.ndrops_prior != summary.ndrops_after, :])
        summary_final.loc[:, 'delta_drops'] = summary_final.ndrops_after - summary_final.ndrops_prior
        summary_final.loc[:, 'perc_delta_drops'] = (summary_final.ndrops_after - summary_final.ndrops_prior) / summary_final.ndrops_prior * 100
        summary_final.reset_index(inplace=True)
        summary_final.rename(columns={'index': 'record number'}, inplace=True)
        # create summary data frame -- STOP

        # return rainobj

        if change_original is False:
            return rainobj, summary_final
        else:
            return summary

# Remove "too narrow" counts: we keep records where at least _nclmin_=4 classes where occupied
    def remove_narrow(self, _nclmin_=4):
        """
        - Purporse: 
        - Return:
        """


        disdroclass_notnarrow=cp.deepcopy(self)
        disdroclass_narrow=cp.deepcopy(self)
        
        disdroclass_narrow.data['nzero_ncl'] = disdroclass_narrow.data[disdroclass_narrow.data != 0].count(axis=1)
        disdroclass_notnarrow.data['nzero_ncl'] = disdroclass_notnarrow.data[disdroclass_notnarrow.data != 0].count(axis=1)
        disdroclass_notnarrow.data = disdroclass_notnarrow.data.loc[disdroclass_notnarrow.data.nzero_ncl>=_nclmin_,:].drop(columns=['nzero_ncl'])
        disdroclass_narrow.data = disdroclass_narrow.data.loc[disdroclass_narrow.data.nzero_ncl<_nclmin_,:].drop(columns=['nzero_ncl'])

        _row_ = [[disdroclass_notnarrow.data.shape[0], disdroclass_notnarrow.data.sum().sum(), round(disdroclass_notnarrow.bulkvar_vplaw.R.sum()), 'legit']]
        _row_.append([disdroclass_narrow.data.shape[0], disdroclass_narrow.data.sum().sum(), round(disdroclass_narrow.bulkvar_vplaw.R.sum()), 'narrow'])
        _summary_ = pd.DataFrame(_row_, columns=['nrecords','ndrops','rainfall_total','type'])

        return (disdroclass_notnarrow, disdroclass_narrow, _summary_)


# Phase Space Parameters method
    def flux_phase_space_parameters(self):
        """
        - Purpose: calculate for each record the following "flux" quantity:
        drop count (N), mean (mu),
        standard deviation (sigma), skewness (gamma), kurtosis (kappa),
        fifth central moment (eta), sixth central moment (omega)
        - Return: data frame with all the above quantities as columns
        """
        list_moments_order = list([1, 2, 3, 4, 5, 6])
        _df_ = self.flux_moment_calculator(list_moments_order)
        _df_['N'] = self.data.sum(axis=1)
        _df_['mu'] = _df_.M1
        _df_['sigma'] = pow((_df_.M2 - pow(_df_.M1, 2)), 0.5)
        _df_['gamma'] = (_df_.M3 + (2 * pow(_df_.M1, 3)) - (3 * _df_.M1 * _df_.M2)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 1.5))
        _df_['kappa'] = (_df_.M4 - (3 * pow(_df_.M1, 4)) + (6 * _df_.M2 * pow(_df_.M1, 2)) - (4 * _df_.M1 * _df_.M3)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2))
        _df_['eta'] = (_df_.M5 + (4 * pow(_df_.M1, 5)) + (10 * _df_.M3 * pow(_df_.M1, 2)) - (10 * _df_.M2 * pow(_df_.M1, 3)) - (5 * _df_.M4 * _df_.M1)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2.5))
        _df_['omega'] = (_df_.M6 - (6 * _df_.M5 * _df_.M1) + (15 * _df_.M4 * pow(_df_.M1, 2)) - (20 * _df_.M3 * pow(_df_.M1, 3)) + (15 * _df_.M2 * pow(_df_.M1, 4)) - (5 * pow(_df_.M1, 6))) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 3))
        _df_.drop(columns=['M1', 'M2', 'M3', 'M4', 'M5', 'M6'], inplace=True)
        return _df_

    def cloud_phase_space_parameters(self, _speed_='plaw'):
        """
        - Purpose: calculate for each record the following "cloud" quantity:
        drop count (N), mean (mu),
        standard deviation (sigma), skewness (gamma), kurtosis (kappa),
        fifth central moment (eta), sixth central moment (omega)
        - Return: data frame with all the above quantities as columns
        """
        list_moments_order = list([0, 1, 2, 3, 4, 5, 6])
        if (_speed_=='expo'):
            _df_ = self.cloud_moment_calculator(list_moments_order,_speed_='expo')
        else:
            _df_ = self.cloud_moment_calculator(list_moments_order,_speed_='plaw')
        _df_['N'] = self.data.sum(axis=1)
        _df_['Nv'] = (round(1 / ((self.instrument_area / 1000000) * self.time_interval) * _df_.N * _df_.M0)).astype(int) 
        _df_['mu'] = _df_.M1
        _df_['sigma'] = pow((_df_.M2 - pow(_df_.M1, 2)), 0.5)
        _df_['gamma'] = (_df_.M3 + (2 * pow(_df_.M1, 3)) - (3 * _df_.M1 * _df_.M2)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 1.5))
        _df_['kappa'] = (_df_.M4 - (3 * pow(_df_.M1, 4)) + (6 * _df_.M2 * pow(_df_.M1, 2)) - (4 * _df_.M1 * _df_.M3)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2))
        _df_['eta'] = (_df_.M5 + (4 * pow(_df_.M1, 5)) + (10 * _df_.M3 * pow(_df_.M1, 2)) - (10 * _df_.M2 * pow(_df_.M1, 3)) - (5 * _df_.M4 * _df_.M1)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2.5))
        _df_['omega'] = (_df_.M6 - (6 * _df_.M5 * _df_.M1) + (15 * _df_.M4 * pow(_df_.M1, 2)) - (20 * _df_.M3 * pow(_df_.M1, 3)) + (15 * _df_.M2 * pow(_df_.M1, 4)) - (5 * pow(_df_.M1, 6))) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 3))
        _df_.drop(columns=['N','M0','M1', 'M2', 'M3', 'M4', 'M5', 'M6'], inplace=True)
        return _df_

# Bulk Variable method
    def bulk_variables(self,_speed_='plaw'):
        """
        - Purpose: calculate for each record the Number of drops (N), the number of drops per unit volume N_V,
        rainfall rate (R), reflectivity  (Z), liquid water content (W), rainfall rate per unit drop (R/N=r),
        reflectivity per unit drop (Z/N_V=z), liquid water content per unit drop (W/N_V=w)
        - Return: data frame with all the above quantities as columns
        """

        PI = 3.141592653589793  # approximate value of greek pi
        seconds_in_hour = 3600  # for converting rainfall rate in mm/h

        # initialize data frame for bulk variables
        NNvbulk = pd.DataFrame(columns=['N', 'Nv', 'R', 'Z', 'W', 'r', 'z', 'w'])
        NNvbulk = NNvbulk.astype(
            dtype={'N': 'int64', 'Nv': 'float64', 'R': 'float64', 'Z': 'float64', 'W': 'float64', 'r': 'float64',
                   'z': 'float64', 'w': 'float64'})

        # number of drops trough the catchment area N
        NNvbulk['N'] = self.data.sum(axis=1)

        # Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
        R_df_=self.flux_moment_calculator([3])
        NNvbulk['R'] = (PI / 6) * (1 / (self.instrument_area * self.time_interval)) * seconds_in_hour * NNvbulk.N * R_df_.M3

        list_moments_order = list([0, 3, 6])
        if (_speed_=='expo'):
            _df_ = self.cloud_moment_calculator(list_moments_order,_speed_='expo')
        else:
            _df_ = self.cloud_moment_calculator(list_moments_order)
        
        # the division by 1000000 is necessary to have catchment area in meter squared
        NNvbulk['Nv'] =(round(1 / ((self.instrument_area / 1000000) * self.time_interval) * NNvbulk.N * _df_.M0)).astype(int)
        NNvbulk['Nv_float'] =(1 / ((self.instrument_area / 1000000) * self.time_interval) * NNvbulk.N * _df_.M0)
        NNvbulk['Z'] = NNvbulk.Nv_float * _df_['M6']
        NNvbulk['W'] = NNvbulk.Nv_float * _df_['M3']
        NNvbulk['r'] = NNvbulk['R'] / NNvbulk['N']
        NNvbulk['z'] = NNvbulk['Z'] / NNvbulk['Nv_float']
        NNvbulk['w'] = NNvbulk['W'] / NNvbulk['Nv_float']

        NNvbulk.drop(columns=['Nv_float'],inplace=True)
        return NNvbulk

# instantaneous spectrum method
    def flux_drop_pdf(self):
        """
        - Purpose: calculate the instantaneous drop pdf
        - Return: data frame with pdf value per class (one row per record)
        """
        dropfreq = cp.deepcopy(self.data)
        dropfreq['ndrops'] = dropfreq.sum(axis=1)
        dropfreq = dropfreq.div(dropfreq['ndrops'].values, axis=0)
        dropfreq.drop(['ndrops'], axis=1, inplace=True)
        class_span = self.classlimits.loc['right', :] - self.classlimits.loc['left', :]
        droppdf = dropfreq/class_span

        return droppdf


# 
    def flux_additional_pdf_paramters(self):
        """
        - Purpose: 
        - Return: 
        """
        
        #  
        pdf_df = self.flux_drop_pdf()
        _df_=cp.deepcopy(self.flux_phase_space_parameters()[['N']])
        

        # mode variables
        idxmax = pdf_df.idxmax(axis=1) # column -- class with the maximum value of the pdf
        pdfmax = pdf_df.max(axis=1) # max value of pdf
        diammax=self.classlimits.loc['left',idxmax] + (self.classlimits.loc['right',idxmax] - self.classlimits.loc['left',idxmax]) / 2
        _df_.loc[:,'class_mode'] = idxmax.values
        _df_.loc[:,'modeclassn'] = (idxmax.str.extract(r'([0-9]+)')).values.astype(int)
        _df_.loc[:,'pdf_mode'] = pdfmax
        _df_.loc[:,'D_mode'] = diammax.values

        # span variables 
        res = pdf_df[pdf_df !=0 ].stack()
        res.rename_axis(index=['record number', 'class'],inplace=True)
        res2 = res.reset_index().drop(columns=[0])
        firstclass = res2.groupby(['record number']).head(1)
        firstclass.set_index('record number', inplace=True)
        firstclass.columns = ['first_notzero_class']
        lastclass = res2.groupby(['record number']).tail(1)
        lastclass.set_index('record number',inplace=True)
        lastclass.columns = ['last_notzero_class']

        _df_ = _df_.join(firstclass)
        _df_ = _df_.join(lastclass)
        _df_['fclassn'] = _df_.first_notzero_class.str.extract(r'([0-9]+)').astype(int)
        _df_['lclassn'] = _df_.last_notzero_class.str.extract(r'([0-9]+)').astype(int)
        _df_['D_span'] = self.classlimits.loc['right',_df_.last_notzero_class].values - self.classlimits.loc['left',_df_.first_notzero_class].values

        #_df_.drop(columns=['first_notzero_class','last_notzero_class'], inplace=True)
        #_df_.columns= ['N','class_mode','pdf_mode','D_mode','first_notzero_class','last_notzero_class','D_span']


        # gradient
        _df_.loc[:,'pln'] = (_df_.modeclassn-4)#.astype(str)
        _df_.loc[:,'prn'] = (_df_.modeclassn+4)#.astype(str)
        _df_.loc[_df_.modeclassn-4<_df_.fclassn, 'pln'] = (_df_.fclassn)#.astype(str)
        _df_.loc[_df_.modeclassn+4>_df_.lclassn, 'prn'] = (_df_.lclassn)#.astype(str)

        _df_.loc[:,'pl'] = 'C' + (_df_.modeclassn-4).astype(str)
        _df_.loc[:,'pr'] = 'C' + (_df_.modeclassn+4).astype(str)
        _df_.loc[_df_.modeclassn-4<_df_.fclassn, 'pl'] = 'C' + (_df_.fclassn).astype(str)
        _df_.loc[_df_.modeclassn+4>_df_.lclassn, 'pr'] = 'C' + (_df_.lclassn).astype(str)
        _df_.reset_index(inplace=True)

        pdf_df_matrix = pdf_df.values

        rn=_df_.index.values
        pdfpl=_df_.pln.values-1
        pdfpr=_df_.prn.values-1
        pdfmode=_df_.modeclassn.values-1

        _df_.loc[:,'delta_pdf_left']=pdf_df_matrix[rn,pdfpl] - pdf_df_matrix[rn,pdfmode]
        _df_.loc[:,'delta_pdf_right']=pdf_df_matrix[rn,pdfpr] - pdf_df_matrix[rn,pdfmode]
        _df_.loc[:,'delta_D_left'] = self.classlimits.loc['left',_df_.class_mode].values + \
            (self.classlimits.loc['right',_df_.class_mode].values - self.classlimits.loc['left',_df_.class_mode].values)/2 - \
                self.classlimits.loc['left',_df_.pl].values
        _df_.loc[:,'delta_D_right'] = self.classlimits.loc['right',_df_.pr].values - \
            (self.classlimits.loc['left',_df_.class_mode].values + \
                (self.classlimits.loc['right',_df_.class_mode].values - self.classlimits.loc['left',_df_.class_mode].values)/2)
        _df_.loc[:,'grad_left'] = _df_.delta_pdf_left / _df_.delta_D_left
        _df_.loc[:,'grad_right'] = _df_.delta_pdf_right / _df_.delta_D_right

        _df_.set_index('record number', inplace=True)
        _df_.drop(columns=['class_mode','first_notzero_class','last_notzero_class','pl','pr','pln','prn','delta_pdf_left',\
            'delta_pdf_right','delta_D_left','delta_D_right'], inplace=True)
        _df_.columns= ['N','class_mode','pdf_mode','D_mode','first_notzero_class','last_notzero_class','D_span','grad_left','grad_right']
                

        return _df_


# renormalized spectrum
    def flux_renormalized_spectrum(self, _bin_=0.2, Dr_left=-5, Dr_right=25):

        _spectrum_ = self.flux_drop_pdf()
        pdf = _spectrum_.values

        cll = cp.deepcopy(self.classlimits)


        _N_ = self.flux_phase_space_parameters()[['N']].values
        _mu_ = self.flux_phase_space_parameters()[['mu']].values
        _sigma_ = self.flux_phase_space_parameters()[['sigma']].values

        _b_ = np.append(cll.loc[['left']].values, cll.loc[['right']].values[0, cll.shape[1]-1])
        _borders_= np.tile(_b_,(_spectrum_.shape[0],1))
        _borders_renorm=(_borders_[:,:]-_mu_[:])/_sigma_[:]


        _nbins_ = np.floor((Dr_right - Dr_left)/_bin_).astype(int)  
        res=np.zeros((_spectrum_.shape[0],_nbins_))

        # print(res.shape)
        for _row_ in range (0,_spectrum_.shape[0]):
            # print(_row_)
            for i in range (0,pdf.shape[1]):
                lb=(_borders_renorm[_row_][i]+5)/_bin_
                lb_int=lb.astype(int)
                lb_rem=lb-lb_int
                rb=(_borders_renorm[_row_][i+1]+5)/_bin_
                rb_int=rb.astype(int)
                rb_rem=rb-rb_int
                # print(i,lb,lb_int,lb_rem,rb,rb_int,rb_rem)
                if (pdf[_row_][i]>0):
                    # print(_row_,i,lb,lb_int,lb_rem,rb,rb_int,rb_rem)
                    for k in range (lb_int+1,rb_int):
                        res[_row_][k] = pdf[_row_][i]*_sigma_[_row_][0]
                        # print(i,k,res[_row_][k])
                    res[_row_][lb_int] = res[_row_][lb_int] + (1-lb_rem)*pdf[_row_][i]*_sigma_[_row_][0]
                    res[_row_][rb_int] = res[_row_][rb_int] + rb_rem*pdf[_row_][i]*_sigma_[_row_][0]
            
            res[_row_][:]=res[_row_][:]*_N_[_row_][0]

        xr = np.arange(Dr_left, Dr_right, _bin_)
        renpdf = res.sum(axis=0)/_N_.sum()

        #print(xr.shape)
        #print(renpdf).shape
        _pd_ = pd.DataFrame({'Dr': xr, 'pdf': renpdf})
        _pd_ = _pd_.loc[_pd_.pdf>0,:]
        return _pd_


