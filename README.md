# Data-Preprocessing

Data pre-processing is an important step in the data mining process. The phrase "garbage in, garbage out" is particularly applicable to data mining and machine learning projects. Data-gathering methods are often loosely controlled, resulting in out-of-range values, impossible data combinations (e.g., Sex: Male, Pregnant: Yes), missing values. Analyzing data that has not been carefully screened for such problems can produce misleading results. Thus, the representation and quality of data is first and foremost before running an analysis.

If there is much irrelevant and redundant information present or noisy and unreliable data, then knowledge discovery during the training phase is more difficult. Data preparation and filtering steps can take considerable amount of processing time. Data pre-processing includes cleaning, normalization, transformation, feature extraction and selection. The product of data pre-processing is the final training set.

The data-preprocessing routines involve standardization(stndze),graphical summary(gs),skewness,kurtosis,creating dummy variables,box cox transformation,etc.


## Feature Selection

Not having the correct and complete data is often the most cited reason for analytics project failures, regardless of Big or Small data. To mitigate the problem, data-driven companies are giving importance to preparing and curating the data, and make it ready for analysis. It is a well-established fact that typically 60-70% of time in any analytics project is spent on data capture and preparation, and hence robust data management tools are important to drive efficiency and time savings. In a Predictive Modeling environment, data preparation is closely associated with the Pre-modeling phase. 

### The different capabilities offered by our Data-Presprocessing package is given below :  

[Standardization](https://github.com/serendio-labs/data-preprocessing-r/wiki/Standardization) - Standardize the raw feature vectors from the training data.

[Deviations](https://github.com/serendio-labs/data-preprocessing-r/wiki/Deviations) - Calculate the deviation of a particular value from the average.

[Indicator Variables](https://github.com/serendio-labs/data-preprocessing-r/wiki/Indicator-Variable) - Create Indicator variables representing the training data.

[Skewness](https://github.com/serendio-labs/data-preprocessing-r/wiki/Skewness) - Compute the skewness of a sample within a training set.

[Kurtosis](https://github.com/serendio-labs/data-preprocessing-r/wiki/Kurtosis) - Compute the kurtosis of a sample within a training set.

[Box-cox Transformation](https://github.com/serendio-labs-stage/diskoveror-datapreprocessing-R/wiki/Box-Cox-Transformation) - Transform the training vectors using Box-cox.

[Poisson Transformation](https://github.com/serendio-labs/data-preprocessing-r/wiki/Poisson-Transformation) - Transform the training vectors using Poisson.

[Proportional Transformation](https://github.com/serendio-labs/data-preprocessing-r/wiki/Proportional-Transformation) - Transform the training vectors with Proportional transformation.

[Graphical Summary](https://github.com/serendio-labs/data-preprocessing-r/wiki/Graphical-summary) - Get a pictorial representation of the training data.

[Anderson Normality](https://github.com/serendio-labs/data-preprocessing-r/wiki/Anderson-Darling-Normality)-test for normality

## Getting started

Download or pull the data-preprocessing package https://github.com/serendio-labs/data-preprocessing-r.git into the appropriate location, then refer to each of the above links to work with the respective utility.



