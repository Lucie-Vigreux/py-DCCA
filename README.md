# py-DCCA
This python code was developped during my internship at LARIS (http://laris.univ-angers.fr/fr/index.html).

This repository allows detrended cross correlation analysis with python.
It is a robust method, introduced by Podobnik and Stanley (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.084102), to quantify the intercorrelation between two non-stationary time signals.
It has been tested on ARFIMA (AutoRegressive Fractional Integral Moving Average) signals to obtain the results of the article by Podobnik and Stanley. The DCCA coefficient dominates the Pearson coefficient for nonstationary series.

Interpretation of the coefficient :  CDCCA(s) = −1 for perfectly anti-correlated series, CDCCA(s) = 1 for perfectly correlated series and CDCCA(s) = 0 for uncorrelated processes.

This package use the numpy package (http://numpy.scipy.org)
![py-DCCA Design](DCCA COEFF.png)
