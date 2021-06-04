# optimize_ego
Fit neural data to a cosine function.
## 1. modelMe
Use position, spike time, and angular information to search for a set of optimal parameters which describe the set of circumstances under which egocentric modulation of the unit is most likely. This function is modeled after the reference-heading model introducted in Pablo Jercog et al.'s work in 2019.

> **Inputs:**
> * **P:** X and Y coordinates of the moving agent (in centimeters). Format in the form [t x(LED1) y(LED1) x(LED2) y(LED2)] where t is a vector of timestamps (in seconds).
> * **ST:** Timestamps for spike events (in seconds).
> * **Z:** Angular variable (in degrees). This is usually the agent's head direction.

> **Outputs:**
> * **out:** Struct with information about ratemaps computed from the date, model-predicted ratemaps, and corresponding measures.
>> * **out.model.Rxyh:** 10x10x10 matrix; model-predicted spatial ratemap which is conditioned on the angular variable, Z. 
>> * **out.model.error:** 1x1 scalar; mean squared error (MSE) between the true and model-predicted conditional ratemaps computed from a cosine function with best-fit parameters.
>> * **out.model.fitParams.g:** 1x1 scalar; optimal value of the fit parameter, g, which can be interpreted as the strength of modulation or height of the cosine curve.
>> * **out.model.fitParams.thetaP:** 1x1 scalar; optimal value of the fit parameter, thetaP, which is the egocentric bearing (or reference-heading) in degrees. This translates into the x-coordinate of the peak of the model-predicted cosine curve.
>> * **out.model.fitParams.xref:** 1x1 scalar; optimal value of the fit parameter, xref, which describes the x-coordinate of the model-predicted egocentric reference point. This output is given in units of spatial bins, where the recording arena is confined to bins 1 through 10, though this parameter is not constrained and can be infinitley distant from the arena.
>> * **out.model.fitParams.yref:** 1x1 scalar; same as above, but the corresponding y-coordinate.
>> * **out.data.Rxyh:** 10x10x10 matrix; spatial ratemap calculated for the real data which is conditioned on the angular variable, Z and normalized by the average firing rate in each spatial bin.
>> * **out.data.rxyh:** 10x10x10 matrix; same as above, but without normalization.
>> * **out.data.rxy:** 10x10 matrix; classic spatial ratemap calculated for the real data.
>> * **out.data.RxyhN:** 10x10x10 matrix; same as *out.data.Rxyh*, but with NaNs replacing spatial+angular bins which do not pass the sampling and firing rate criteria set by the reference-heading model.
>> * **out.data.rxyhN:** 10x10x10 matrix; same as *out.data.rxyh* with NaN replacement (as above).
>> * **out.data.rxyN:** 10x10x10 matrix; same as *out.data.rxy* with NaN replacement (as above).
>> * **out.data.occ.count:** 10x10x10 matrix; spatial+angular occupancy expressed in *counts*.
>> * **out.data.occ.time:** 10x10x10 matrix; spatial+angular occupancy expressed in *seconds*.
>> * **out.measures.VE.place:** 1x1 scalar; proxy for the variance explained (VE) by place tuning.
>> * **out.measures.VE.RH:** 1x1 scalar; proxy for the variance explained (VE) by the reference-heading (RH) model.
>> * **out.measures.MVL.HD:** 10x10 matrix; mean vector length of the tuning curve for the true values of the angular variable (*e.g.-* head direction; HD) in each of 100 spatial bins. This will be used to generate a vector field with the *plotMe* function.
>> * **out.measures.MVL.RH:** 10x10 matrix; same as above, but for [RH] model-predicted tuning curves.
>> * **out.measures.TS.HD:** 1x1 scalar; tuning strength (TS) of the unit to the angular variable (*e.g.-* head direction; HD). This is taken as the linear average of *out.measures.MVL.HD*.
>> * **out.measures.TS.RH:** 1x1 scalar; same as above, but for *out.measures.MVL.RH*.
>> * **out.measures.mu.HD:** 10x10 matrix; preferred local head direction orientation. This is calculated circularly using the true angular tuning curves in each of 100 spatial bins.
>> * **out.measures.mu.RH:** 10x10 matrix; same as above, but for [RH] model-predicted tuning curves.