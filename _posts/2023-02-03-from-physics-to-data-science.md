---
layout: post
title: From Physics to data science
---

![SuperCDMS data scheme]({{ 'assets/img/SuperCDMS_data.png' | relative_url}})

Big data is the backbone of modern particle physics. Did you know that [CERN](https://home.cern/science/computing/storage) produces millions of gigabytes of data in an average experimental run? While CERN smashes highly energetic particle beams to produce a large number of collision events (and thus a large amount of data), dark matter searches are typically what are called "rare event searches" in the particle physics community. Rare event searches, as the name suggests, seek to find extremely rare physical phenomena in the universe, such as a dark matter particle scattering off a crystal in a deep underground laboratory. By definition, you would naturally expect less data from a rare event search experiment compared to a collider experiment such as the [LHC](https://home.cern/science/accelerators/large-hadron-collider). Yet, my experiment, [SuperCDMS](https://supercdms.slac.stanford.edu/) e.g., is expected to produce ~100 TB of data in a year. To provide a sense of how massive this is, the training data for [GPT-3](https://arxiv.org/pdf/2005.14165.pdf) was 570 GB post filtering of the acquired 45 TB of data from 2016-2019.

So what does our data look like? Our raw data acquired from the detectors looks are digitized waveforms sampled at 625 kHz rate, i.e., time series data. It could look like this:

![Example pulse 1]({{ 'assets/img/pulse1.png' | relative_url}})

or, this:

![Example pulse 2]({{ 'assets/img/pulse2.png' | relative_url}})

or even this:

![Example pulse 3]({{ 'assets/img/pulse3.png' | relative_url}})

Can you spot the pulse in the last one? Regression, PCA and even classifiers with CNNs are approaches used to study data of this kind. We also often process these data sets to decompose it into physically interpretible and human analyzable components. Processing maps the raw pulses you saw above into a feature space. Reduced features (called reduced quantities in SuperCDMS) are, thus, another form of data we study. This could be the location of events clustered in a channel like this ([image related to this paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.022002)):

![Triangle plot]({{ 'assets/img/Triangle_plot.png' | relative_url}})

or the correlation between various detector channels like this:

![Correlation matrix between channels]({{ 'assets/img/NxMHybridCorrelationcomparison.png' | relative_url}})

or the event energy distribution ([image related to this paper for example](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.062001)) like this:

![Energy spectrum example]({{ 'assets/img/Energy_spectrum.png' | relative_url}})

Big data calls for extensive and creative data analysis techniques, making every particle physicist also a data scientist.

Want to know how machine learning can help detect dark matter? Check out [this project page]({{'/dark_matter' | relative_url}})