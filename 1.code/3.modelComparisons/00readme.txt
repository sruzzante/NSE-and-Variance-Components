Author: Sacha Ruzzante
sachawruzzante@gmail.com
Last Update: 2025-06-18

Technical Note: The Nash Sutcliffe Efficiency conceals inferior representations of hydrologic variability in seasonal regimes

These scripts calculate all the goodness of fit statistics for existing hydrologic models, which are used to produce Figure 4.

The models come from the following publications:

Arsenault, R., Martel, J.-L., Brunet, F., Brissette, F., & Mai, J. (2023). Continuous streamflow prediction in ungauged basins: Long short-term memory neural networks clearly outperform traditional hydrological models. Hydrology and Earth System Sciences, 27(1), 139–157. https://doi.org/10.5194/hess-27-139-2023

Kraft, B., Schirmer, M., Aeberhard, W. H., Zappa, M., Seneviratne, S. I., & Gudmundsson, L. (2025). CH-RUN: A deep-learning-based spatially contiguous runoff reconstruction for Switzerland. Hydrology and Earth System Sciences, 29(4), 1061–1082. https://doi.org/10.5194/hess-29-1061-2025

Kratzert, F., Gauch, M., Klotz, D., & Nearing, G. (2024). HESS Opinions: Never train a Long Short-Term Memory (LSTM) network on a single basin. Hydrology and Earth System Sciences, 28(17), 4187–4201. https://doi.org/10.5194/hess-28-4187-2024

Kratzert, F. (2019). CAMELS benchmark models [Dataset]. HydroShare. https://doi.org/10.4211/hs.474ecc37e7db45baa425cdb4fc1b61e1

Seibert, J., Vis, M. J. P., Lewis, E., & van Meerveld, H. j. (2018). Upper and lower benchmarks in hydrological modelling. Hydrological Processes, 32(8), 1120–1125. https://doi.org/10.1002/hyp.11476

Mizukami, N., Rakovec, O., Newman, A. J., Clark, M. P., Wood, A. W., Gupta, H. V., & Kumar, R. (2019). On the choice of calibration metrics for “high-flow” estimation using hydrologic models. Hydrology and Earth System Sciences, 23(6), 2601–2614. https://doi.org/10.5194/hess-23-2601-2019

Newman, A. J., Mizukami, N., Clark, M. P., Wood, A. W., Nijssen, B., & Nearing, G. (2017). Benchmarking of a Physically Based Hydrologic Model. https://doi.org/10.1175/JHM-D-16-0284.1

The model for Brazil was developed for this publication - see the folder /lstm-camelsbr/
The model for Arsenault et al (2023) was rerun using the provided code at https://osf.io/3s2pq/
