# NSE-and-Variance-Components
Author: Sacha Ruzzante
Last Updated: 2025-08-06

This code was written to produce the analyses and plots in the manuscript:

Ruzzante, S. W., Knoben, W. J. M., Wagener, T., Gleeson, T., and Schnorbus, M. Technical Note: High Nash Sutcliffe Efficiencies conceal poor simulations of interannual variance in tropical, alpine, and polar catchments. Submitted to Hydrology and Earth System Sciences, 2025

Streamflow time series can be decomposed into interannual, seasonal, and irregular components, with regionally varying contributions of each component. Seasonal variance dominates in many tropical, alpine, and polar regions, while irregular variance dominates in most other regions. Year-to-year variations in streamflow are known to strongly influence human and ecological systems, though we find that historical interannual variance is usually only a small fraction of the total variance. We show that hydrologic models often simulate one component well while failing to simulate the others, a fact that is hidden by popular performance metrics such as the Nash-Sutcliffe Efficiency (NSE) and the Kling-Gupta Efficiency (KGE) which aggregate performance to a single number. We analyse 18 regional and global hydrologic models and found that in highly seasonal catchments where the NSE and KGE are consistently the highest, the interannual performance is almost always worse. The NSE of the interannual component is lower in highly seasonal catchments, and year-to-year changes in ecologically relevant hydrologic signatures are more poorly simulated. This is concerning because it indicates that these hydrologic models may struggle to predict long-term responses to climate change, especially in tropical, alpine, and polar regions, which are some of the most vulnerable regimes regarding climate change. 

The repository is organized as follows:

[/1.code/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/1.code) contains all R codes

[/1.code/1.climatologicalBenchmarkGOF/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/1.code/1.climatologicalBenchmarkGOF/) contains the scripts to analyse streamflow data and calculate the variance components and benchmark NSEs for each catchment

[/1.code/3.modelComparison/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/1.code/3.modelComparisons) contains the scripts to calculate goodness-of-fit statistics for all 18 models analysed. 

[/1.code/4.Plotting/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/1.code/4.Plotting) contains the scripts to make all four figures in the manuscript 

[/1.code/5.utils/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/1.code/5.utils) contains utility functions used throughout other scripts

[/2.data/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/2.data) contains data produced by this project 

[/2.data/GOF/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/2.data/GOF) contains goodness-of-fit stats for the climatological benchmark model,  (Figure 2)

[/2.data/highLowBenchmarkGOF/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/2.data/highLowBenchmarkGOF) contains goodness-of-fit statistics for each of the 18 models.

[/2.data/NHM_Q/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/2.data/NHM_Q) contains simulation data from the National Hydrologic Model

[/2.data/varComponents/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/2.data/varComponents) contains the variance components for 17,245 gauges (Figure 2)

[/2.data/worldclim/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/2.data/worldclim/indices) includes climate indices calculated from WorldClim for each gauges

[/3.figures/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/3.figures) is where figures will be saved

[/4.lstm-camelsbr/](https://github.com/sruzzante/NSE-and-Variance-Components/tree/main/4.lstm-camelsbr) is where the LSTM model for the camels-br data is saved. This was created using the NeuralHydrology package for Python.