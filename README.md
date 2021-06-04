# optimize_ego
Fit neural data to a cosine function.
## 1. modelMe
Use position, spike time, and angular information to search for a set of optimal parameters which describe the set of circumstances under which egocentric modulation of the unit is most likely. This function is modeled after the reference-heading model introducted in Pablo Jercog et al.'s work in 2019.
Inputs:
> * P:   X and Y coordinates of the moving agent (in centimeters). Format in the form [t x(LED1) y(LED1) x(LED2) y(LED2)] where t is a vector of timestamps (in seconds).
> * ST:  Timestamps for spike events (in seconds).
> * Z:   Angular variable (in degrees). This is usually the agent's head direction.
Outputs:
> * out: Struct with information about ratemaps computed from the date, model-predicted ratemaps, and corresponding measures.
>> * out.model.Rxyh: 10x10x10 model-predicted spatial ratemap which is conditioned on the angular variable, Z. 
>> * out.model.error: 1x1 scalar describing the mean squared error (MSE) between the true and model-predicted conditional ratemaps computed from a cosine function with best-fit parameters.
>> * out.model.fitParams.g: 1x1 scalar describing the optimal value of the fit parameter, g, which can be interpreted as the strength of modulation or height of the cosine curve.
>> * out.model.fitParams.thetaP: 1x1 scalar describing the optimal value of the fit parameter, thetaP, which is the egocentric bearing (or reference-heading) in degrees. This translates into the x-coordinate of the peak of the model-predicted cosine curve.
