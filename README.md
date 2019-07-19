# <font color=9ACD32>disdrorain</font>

Python package for analyzing disdrometer data. To use it just import __disdrorain.py__ in your python code, e.g. :

```python
import disdrorain as dr
```



### <font color=9ACD32>disdrorain class</font>:

An instance of the disdrorain class is created by specifying a

- __datapath__:  system path to disdrometer data (CSV format, no header). (default=None).  Alternatively, you can specify the name of an existing data frame via the dataframe option. If datapath is used do not use dataframe.
- __dataframe__: name of an existing data frame. (default=None).  Alternatively, you can specify a system path to an existing file via the datapath option. If datafile is used do not use datapath.
- __classpath__: system path to disdrometer cell limits file (CSV format, no header, 2columns: 1st=left border, 2nd=right border). (default=None). If no classpath is specified, the Joss-Valdvogel RD 80 disdrometer cell limits are adopted.

- __instrument_area__: the catchment are of the disdrometer in mm^2^  (default=5000: typical value of a Joss-Valdvogel impact disdrometer).

- __time_interval__: the time integration of the disdrometer in seconds (default=60: typical value of a Joss-Valdvogel impact disdrometer).

  __Cspeed__ (default=3.776)  and __bspeed__ (default=0.67):  coefficients for the rain drop terminal velocity formula $v(D)=CD^{b}$ in m/s.  This coefficients are obtained from Eq.(8) of [AtlasUlbrich]([https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281977%29016%3C1322%3APAAIRM%3E2.0.CO%3B2](https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450(1977)016<1322%3APAAIRM>2.0.CO%3B2).

- __fieldsep__ (default= ' '): The field separator for the disdrometer data set.

Examples:

```python
import disdrorain as dr

# The myrain instance of the disdrorain class is created by using '/apath/disdrodata.csv' as disdrometer data. Default values (Joss-Valdvogel RD80 disdrometer) are used for diameter class limits 
myrain = dr.disdrorain(datapath='/apath/disdrodata.csv')

# The myrain instance of the disdrorain class is created by using '/apath/disdrodata.csv' as disdrometer data. Diameter class limits are specified in the file '/apath/classlimits.csv'
myrain = dr.disdrorain(datapath='/apath/disdrodata.csv',classpath='/apath/classlimits.csv')

# The myrain instance of the disdrorain class is created by using the data frame _df_  as disdrometer data. Default values (Joss-Valdvogel RD80 disdrometer) are used for diameter class limits 
myrain = dr.disdrorain(dataframe=_df_)

# print first 5 rows of the myrain disdrometer data
myrain.data.head(5)
```

### <font color=9ACD32>Class Attributes</font>

- .data = panda data frame with disdrometer counts.
- .classlimits = panda data frame with diameter class limits.
- .instrument_area = catchment area of disdrometer (mm^2^).
- .time_interval = integration time of disdrometer (seconds).
- .Cspeep & .bspeed = coefficients of the drop terminal velocity formula $v(D)=CD^{b}$ in m/s.

### <font color=9ACD32>Class Methods</font>

- moment_calculator: calculate ground drop size distribution moment (any order).
- outlier_deletion: eliminate outliers count (see disdrometer.py for details). 
- phase_space_parameters: calculate phase space parameters of the ground spectrum for each record (integration time interval): Number of drops, mean, standard deviation, skewness, kurtosis, fifth central moment, sixth central moment). Note these parameters refers to the drop size distribution at the ground. See [this article](https://journals.ametsoc.org/doi/citedby/10.1175/JAMC-D-13-050.1) for reference.
- bulk_variables: calculate rainfall but variable for each record: Number of drops @ ground, Number of drops @ cloud, Liquid Water Content, Rainfall Rate, Reflectivity.
- drop_pdf: ground drop size distribution. Namely $p_{j}=n_{j}/\Delta_{j}$, where $n_{j}$ and $\Delta_{j}$ are, respectively, the drop count and the size of the $j-$th diameter class.
