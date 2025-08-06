Author: Sacha Ruzzante
sachawruzzante@gmail.com
Last Update: 2025-06-18


Technical Note: High Nash Sutcliffe Efficiencies conceal poor simulations of interannual variance in tropical, alpine, and polar catchments

These scripts calculate goodness of fit statistics (NSE, KGE, and KGE components) for 19 large-sample hydrology datasets. The scripts also decompose the time series into interannual, seasonal, and irregular components, using three techniques: classical decomposition, a variation of STL, and a decomposition based on a fast fourier transform, as described in the manuscript. 

All scripts follow a similar pattern. The code to read in data from each source differs.

The data are from the following sources:

Alvarez-Garreton, C., Mendoza, P. A., Boisier, J. P., Addor, N., Galleguillos, M., Zambrano-Bigiarini, M., Lara, A., Puelma, C., Cortes, G., Garreaud, R., McPhee, J., & Ayala, A. (2018). The CAMELS-CL dataset: Catchment attributes and meteorology for large sample studies – Chile dataset. Hydrology and Earth System Sciences, 22(11), 5817–5846. https://doi.org/10.5194/hess-22-5817-2018
Arsenault, R., Martel, J.-L., Brunet, F., Brissette, F., & Mai, J. (2023). Continuous streamflow prediction in ungauged basins: Long short-term memory neural networks clearly outperform traditional hydrological models. Hydrology and Earth System Sciences, 27(1), 139–157. https://doi.org/10.5194/hess-27-139-2023
Casado Rodríguez, J. (2023). CAMELS-ES: Catchment Attributes and Meteorology for Large-Sample Studies – Spain (Version 1.0.2) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.8428374
Chagas, V. B. P., Chaffe, P. L. B., Addor, N., Fan, F. M., Fleischmann, A. S., Paiva, R. C. D., & Siqueira, V. A. (2020). CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil. Earth System Science Data, 12(3), 2075–2096. https://doi.org/10.5194/essd-12-2075-2020
Coxon, G., Addor, N., Bloomfield, J. P., Freer, J., Fry, M., Hannaford, J., Howden, N. J. K., Lane, R., Lewis, M., Robinson, E. L., Wagener, T., & Woods, R. (2020). CAMELS-GB: Hydrometeorological time series and landscape attributes for 671 catchments in Great Britain. Earth System Science Data, 12(4), 2459–2483. https://doi.org/10.5194/essd-12-2459-2020
Delaigue, O., Guimarães, G. M., Brigode, P., Génot, B., Perrin, C., Soubeyroux, J.-M., Janet, B., Addor, N., & Andréassian, V. (2024). CAMELS-FR dataset: A large-sample hydroclimatic dataset for France to explore hydrological diversity and support model benchmarking. Earth System Science Data Discussions, 1–27. https://doi.org/10.5194/essd-2024-415
Dolich, A., Maharjan, A., Mälicke, M., Manoj J, A., & Loritz, R. (2025). Caravan-DE: Caravan extension Germany - German dataset for large-sample hydrology [Dataset]. Zenodo. https://zenodo.org/records/14755229
Efrat, M. (2025). Caravan extension Israel—Israel dataset for large-sample hydrology [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.15003600
Fowler, K. J. A., Zhang, Z., & Hou, X. (2024). CAMELS-AUS v2: Updated hydrometeorological timeseries and landscape attributes for an enlarged set of catchments in Australia. Earth System Science Data Discussions, 1–21. https://doi.org/10.5194/essd-2024-263
Global Runoff Data Centre, The (2025) 56068 Koblenz, Germany. https://grdc.bafg.de/
Helgason, H. B., & Nijssen, B. (2024). LamaH-Ice: LArge-SaMple DAta for Hydrology and Environmental Sciences for Iceland. Earth System Science Data, 16(6), 2741–2771. https://doi.org/10.5194/essd-16-2741-2024
Höge, M., Kauzlaric, M., Siber, R., Schönenberger, U., Horton, P., Schwanbeck, J., Floriancic, M. G., Viviroli, D., Wilhelm, S., Sikorska-Senoner, A. E., Addor, N., Brunner, M., Pool, S., Zappa, M., & Fenicia, F. (2023). CAMELS-CH: Hydro-meteorological time series and landscape attributes for 331 catchments in hydrologic Switzerland. Earth System Science Data, 15(12), 5755–5784. https://doi.org/10.5194/essd-15-5755-2023
Klingler, C., Schulz, K., & Herrnegger, M. (2021). LamaH-CE: LArge-SaMple DAta for Hydrology and Environmental Sciences for Central Europe. Earth System Science Data, 13(9), 4529–4565. https://doi.org/10.5194/essd-13-4529-2021
Lammers, R. B., & Shiklomanov, A. I. (2000). R-ArcticNet, A Regional Hydrographic Data Network for the Pan-Arctic Region. [Dataset]. https://www.r-arcticnet.sr.unh.edu/v4.0/AllData/index.html
Liu, J., Koch, J., Stisen, S., Troldborg, L., Højberg, A. L., Thodsen, H., Hansen, M. F. T., & Schneider, R. J. M. (2024). CAMELS-DK: Hydrometeorological Time Series and Landscape Attributes for 3330 Catchments in Denmark. Earth System Science Data Discussions, 1–30. https://doi.org/10.5194/essd-2024-292
Loritz, R., Dolich, A., Acuña Espinoza, E., Ebeling, P., Guse, B., Götte, J., Hassler, S. K., Hauffe, C., Heidbüchel, I., Kiesel, J., Mälicke, M., Müller-Thomy, H., Stölzle, M., & Tarasova, L. (2024). CAMELS-DE: Hydro-meteorological time series and attributes for 1582 catchments in Germany. Earth System Science Data, 16(12), 5625–5642. https://doi.org/10.5194/essd-16-5625-2024
Mangukiya, N. K., Kumar, K. B., Dey, P., Sharma, S., Bejagam, V., Mujumdar, P. P., & Sharma, A. (2025). CAMELS-IND: Hydrometeorological time series and catchment attributes for 228 catchments in Peninsular India. Earth System Science Data, 17(2), 461–491. https://doi.org/10.5194/essd-17-461-2025
Newman, A. J., Clark, M. P., Sampson, K., Wood, A., Hay, L. E., Bock, A., Viger, R. J., Blodgett, D., Brekke, L., Arnold, J. R., Hopson, T., & Duan, Q. (2015). Development of a large-sample watershed-scale hydrometeorological data set for the contiguous USA: Data set characteristics and assessment of regional variability in hydrologic model performance. Hydrology and Earth System Sciences, 19(1), 209–223. https://doi.org/10.5194/hess-19-209-2015
Nijzink, J., Loritz, R., Gourdol, L., Zoccatelli, D., Iffly, J. F., & Pfister, L. (2024). CAMELS-LUX: Highly Resolved Hydro-Meteorological and Atmospheric Data for Physiographically Characterized Catchments around Luxembourg [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.13846620

