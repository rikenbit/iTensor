[![DOI](https://zenodo.org/badge/402303422.svg)](https://zenodo.org/badge/latestdoi/402303422)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/iTensor)](
https://cran.r-project.org/package=iTensor)
[![Downloads](https://cranlogs.r-pkg.org/badges/iTensor)](https://CRAN.R-project.org/package=iTensor)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/iTensor?color=orange)](https://CRAN.R-project.org/package=iTensor)
[![:name status badge](https://rikenbit.r-universe.dev/badges/:name)](https://rikenbit.r-universe.dev)
[![:registry status badge](https://rikenbit.r-universe.dev/badges/:registry)](https://rikenbit.r-universe.dev)
[![:total status badge](https://rikenbit.r-universe.dev/badges/:total)](https://rikenbit.r-universe.dev)
[![iTensor status badge](https://rikenbit.r-universe.dev/badges/iTensor)](https://rikenbit.r-universe.dev)
![GitHub Actions](https://github.com/rikenbit/iTensor/actions/workflows/build_test_push.yml/badge.svg)
[![status](https://joss.theoj.org/papers/cffccfbad7fc8bda9395e429c3c6f35d/status.svg)](https://joss.theoj.org/papers/cffccfbad7fc8bda9395e429c3c6f35d)

# iTensor
ICA-based Matrix/Tensor Decomposition

Installation
======
~~~~
git clone https://github.com/rikenbit/iTensor/
R CMD INSTALL iTensor
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/iTensor")
~~~~

References
======
- **ICA**
  - *InfoMax*
    - Bell, A. J. et al., An information-maximization approach to blind separation and blind deconvolution. Neural computation, 7(6), 1129-1159, 1995
    - Amari, S. et al., A new learning algorithm for blind signal separation. NIPS 1995, 1995
  - *ExtInfoMax*
    - Lee, T. W., et al., Independent component analysis using an extended infomax algorithm for mixed subgaussian and supergaussian sources. Neural computation, 11(2), 417-441, 1999
  - *FastICA*
    - Hyvarinen, A. Fast and robust fixed-point algorithms for independent component analysis. IEEE transactions on Neural Networks, 10(3), 626-634, 1999
  - *JADE*
    - Cardoso, J. F. et al., Blind beamforming for non-gaussian signals, IEE Proceedings F, 140(6), 362-370, 1993
  - *AuxICA1/2*
    - Ono, N. et al., Auxiliary-Function-Based Independent Component Analysis for Super-Gaussian Sources, Lecture Notes in Computer Science, 6365, 165-172, 2010
  - *IPCA*
    - Yao, F. et al., Independent Principal Component Analysis for biologically meaningful dimension reduction of large biological data sets, BMC Bioinformatics, 13(24), 2012
  - *SIMBEC*
    - Cruces, S. et al., Criteria for the simultaneous blind extraction of arbitrary groups of sources, International Conference on ICA and BSS, 740-745, 2001
  - *AMUSE*
    - Tong, L. et al., Indeterminacy and identifiability of blind identification, IEEE Transactions on Circuits and Systems, 38(5), 499-509, 1991
  - *SOBI*
    - Belouchrani, A. et al., A blind source separation technique using second-order statistics, IEEE Transactions on Signal Processing, 45(2), 434-444, 1997
  - *FOBI*
    - Cardoso, J.-F. et al., Source separation using higher order moments, International Conference on Acoustics, Speech, and Signal Processing, 4, 2109-2112, 1989
  - *ProDenICA*
    - Hastie, T. et al.,  Independent Components Analysis through Product Density Estimation, NIPS 2002, 2002
  - *RICA*
    - Le, Q. et al., ICA with Reconstruction Cost for Efficient Overcomplete Feature Learning, NIPS 2011, 2011
- **GroupICA**
  - Calhourn V. D. et al, A review of group ICA for fMRI data and ICA for joint inference of imaging, genetic, and ERP data. Neuroimage. 45(1 Suppl), S163-72, 2009
  - Pfister, N. et al., groupICA: Independent component analysis for grouped data. arXiv, 2018
- **MICA**
  - Akaho, S. et al., MICA: Multimodal independent component analysis. IJCNN'99, 2, 927-932, 1999
- **MultilinearICA**
  - Vasilescu, M. A. O. et al., Multilinear Independent Component Analysis, IEEE CVPR 2005, 2005
- **CorrIndex**
  - Sobhani, E. et al., CorrIndex: a permutation invariant performance index, Signal Processing, 195, 108457, 2022

## Contributing

If you have suggestions for how `iTensor` could be improved, or want to report a bug, open an issue! We'd love all and any contributions.

For more, check out the [Contributing Guide](CONTRIBUTING.md).

## Authors
- Koki Tsuyuzaki
