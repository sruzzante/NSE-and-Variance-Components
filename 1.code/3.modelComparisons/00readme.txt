Author: Sacha Ruzzante
sachawruzzante@gmail.com
Last Update: 2025-08-06

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

Yang, Y., Feng, D., Beck, H. E., Hu, W., Abbas, A., Sengupta, A., Delle Monache, L., Hartman, R., Lin, P., Shen, C., & Pan, M. (2025). Global Daily Discharge Estimation Based on Grid Long Short-Term Memory (LSTM) Model and River Routing. Water Resources Research, 61(6), e2024WR039764. https://doi.org/10.1029/2024WR039764

Klingler, C., Schulz, K., and Herrnegger, M.: LamaH-CE: LArge-SaMple DAta for Hydrology and Environmental Sciences for Central Europe, Earth Syst. Sci. Data, 13, 4529–4565, https://doi.org/10.5194/essd-13-4529-2021, 2021.

Song, Y., Bindas, T., Shen, C., Ji, H., Knoben, W. J. M., Lonzarich, L., Clark, M. P., Liu, J., van Werkhoven, K., Lamont, S., Denno, M., Pan, M., Yang, Y., Rapp, J., Kumar, M., Rahmani, F., Thébault, C., Adkins, R., Halgren, J., … Lawson, K. (2025). High-Resolution National-Scale Water Modeling Is Enhanced by Multiscale Differentiable Physics-Informed Machine Learning. Water Resources Research, 61(4), e2024WR038928. https://doi.org/10.1029/2024WR038928

Chagas, V. B. P., Chaffe, P. L. B., Addor, N., Fan, F. M., Fleischmann, A. S., Paiva, R. C. D., & Siqueira, V. A. (2020). CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil. Earth System Science Data, 12(3), 2075–2096. https://doi.org/10.5194/essd-12-2075-2020

Nearing, G., Cohen, D., Dube, V., Gauch, M., Gilon, O., Harrigan, S., Hassidim, A., Klotz, D., Kratzert, F., Metzger, A., Nevo, S., Pappenberger, F., Prudhomme, C., Shalev, G., Shenzis, S., Tekalign, T. Y., Weitzner, D., & Matias, Y. (2024). Global prediction of extreme floods in ungauged watersheds. Nature, 627(8004), 559–563. https://doi.org/10.1038/s41586-024-07145-1

Regan, R. S., Juracek, K. E., Hay, L. E., Markstrom, S. L., Viger, R. J., Driscoll, J. M., LaFontaine, J. H., & Norton, P. A. (2019). The U. S. Geological Survey National Hydrologic Model infrastructure: Rationale, description, and application of a watershed-scale model for the conterminous United States. Environmental Modelling & Software, 111, 192–203. https://doi.org/10.1016/j.envsoft.2018.09.023

Schnorbus, M. (2018). VIC Glacier (VIC-GL)—Description of VIC model changes and upgrades, VIC Generation 2 Deployment Report (p. 40). Pacific Climate Impacts Consortium, University of Victoria. https://dspace.library.uvic.ca/server/api/core/bitstreams/850d5363-e326-4b13-bddc-f5cd6227c928/content


The model for Brazil was developed for this publication - see the folder /lstm-camelsbr/
The model for Arsenault et al (2023) was rerun using the provided code at https://osf.io/3s2pq/
For the Yang et al (2025) model you need to source all the observation data yourself, then run part 1 to compile it.
