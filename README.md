# Code for expected risk estimates from IVF and egg donor data.

That's only the code for the main calculations of the risk reduction, as the full code was many different Rmd files printing tables in order to be able to recreate the different plots without
direct access to the data.

utils.R defines the functions used for the calculations, while egg_donor.R and ivf_data.R are where the calculatios are done for each dataset.

Finally, stop_model.R was an attempt to look at a model where after each embryo transfer there was a probability of stopping, and looking if that alone could explain the difference between the group that had untransferred embryos and the one without. While it seems to produce similar data, we didn't include it in the paper in the end.