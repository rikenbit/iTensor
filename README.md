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
    - Bell, Anthony J., and Terrence J. Sejnowski. "An information-maximization approach to blind separation and blind deconvolution." Neural computation 7.6 (1995): 1129-1159.
    - Amari, Shun-ichi, Andrzej Cichocki, and Howard Hua Yang. "A new learning algorithm for blind signal separation." Advances in neural information processing systems. Morgan Kaufmann Publishers, 1996.
  - *ExtInfoMax*: Lee, Te-Won, Mark Girolami, and Terrence J. Sejnowski. "Independent component analysis using an extended infomax algorithm for mixed subgaussian and supergaussian sources." Neural computation 11.2 (1999): 417-441.
  - *FastICA*: Hyvarinen, Aapo. Fast and robust fixed-point algorithms for independent component analysis. IEEE transactions on Neural Networks 10.3 (1999): 626-634.
- **MICA**: Akaho, Shotaro, Yasuhiko Kiuchi, and Shinji Umeyama. MICA: Multimodal independent component analysis. IJCNN'99. International Joint Conference on Neural Networks. Proceedings (Cat. No. 99CH36339). Vol. 2. IEEE, 1999.
- **MultilinearICA**: M. Alex O. Vasilescu, and Demetri Terzopoulos, Multilinear Independent Component Analysis, IEEE CVPR2005, 2005.

## Authors
- Koki Tsuyuzaki
