MODEL-BASED BAYESIAN DIRECTION OF ARRIVAL ANALYSIS FOR SOUND SOURCES USING A SPHERICAL MICROPHONE ARRAY

This algorithm estimates the directions of arrival (DoA) and relative levels of an arbitrary number of sound sources from sound data recorded on Rensselaer Polytechnic Institute’s (RPI) 16-channel spherical microphone array. This program is built upon a framework of Bayesian probability to estimate DoAs accurately. Nested and slice sampling methods were leveraged to execute this framework. Spherical harmonics and beamforming were utilized to spatially filter the incoming sound data as well as construct models. 

NOTE: Functions and other data files not included. Please contact me if you would like access to these files and/or the associated research papers/Thesis. crlandschoot@gmail.com

NOTE: This program was constructed for a specific use in the context of research and will take several changes to function within another ecosystem. For example, the main algorithm only works with RPI's 16-channel spherical microphone array. Its sampling constants are hardcoded and would need to be changed to use another sound source. The purpose of this repo is to provide the framework of the algorithm as a demonstration.

This work is the intellectual property of Christopher Landschoot, Ning Xiang, and Rensselaer Polytechnic Institute.



ABSTRACT - http://dx.doi.org/10.1121/1.5138126

In many room acoustics and noise control applications, it is often challenging to determine the directions of arrival (DoAs) of incoming sound sources. This work seeks to solve this problem reliably by beamforming, or spatially filtering, incoming sound data with a spherical microphone array via a probabilistic method. When estimating the DoA, the signal under consideration may contain one or multiple concurrent sound sources originating from different directions. This leads to a two-tiered challenge of first identifying the correct number of sources, followed by determining the directional information of each source. To this end, a probabilistic method of model-based Bayesian analysis is leveraged. This entails generating analytic models of the experimental data, individually defined by a specific number of sound sources and their locations in physical space, and evaluating each model to fit the measured data. Through this process, the number of sources is first estimated, and then the DoA information of those sources is extracted from the model that is the most concise to fit the experimental data. This paper will present the analytic models, the Bayesian formulation, and preliminary results to demonstrate the potential usefulness of this model-based Bayesian analysis for complex noise environments with potentially multiple concurrent sources.
