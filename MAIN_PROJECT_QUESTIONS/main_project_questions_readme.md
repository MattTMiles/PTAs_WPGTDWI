Matt Miles CW investigatory project

-----------------------

WHERE ARE WE AT?

-----------------------

Okay, a lot has been done in this project now to set up the infrastructure to investigate everything. Now, in these directories I want to do more of a target approach to some very discrete questions.

-----------------------

WHAT ARE THE QUESTIONS?

-----------------------

Question 1: How does detection & localisation change with pulsar distance constraints?

How to do this?

- Define Baseline 

- - Ring of 30 pulsars
- - One CW at moderate frequency ~ log10=-8 (so evolution happens but isn't insane)
- - Phases: WN only to find CW SNR -> Add RN -> Add a GWB
- - SNR as a control parameter (5, 10, 20 ?)


Test different scenarios as a function of SNR and distance information:

SNR at: 5, 10, 20
Distance priors at: 
- Real: Current priors 
- Tier 1: Loosened. This shows what the CW can do alone, and quantifies degeneracy breaking from the CW phase evolution. Also helpful for future PTA puslars. (3x sigma? 10 for extreme?)
- Tier 2: Tightened for various precision increases (1.6, 3, 5?). This shows the increases that we'll get over improved measurements. GAIA alone will give a factor of 1.6.
- Tier 3: Exact, dream scenario.

Across these, measure: 
- Detection statistic
- - Either BF or delta_logL? Or using Dave's optimiser maybe we can get something better?
- Sky localisation area
- - 90% credible area of the localisation pars. Compare to whether the sky position is actually in there. Credible level at the injected precision is probably good actually.
- Posterior sigma(dist) / inj_dist
- - I think this is pretty much the same as "how well did the posterior get the distance".
- - Another option might be prior sigma(dist) / inj_dist vs posterior sigma(dist) / inj_dist, because that'll tell us how much the CW improved the distance measurement. I think definitely do this as well.


Ring PTA to real PTA relationship needs to be respected I think. 
I think the way to do this is to draw from the distribution of fractional uncertainty of the pulsars as they currently stand. In this way, when we take the pulsars from our larger sample we retain that distribution. This feels "fair"?


Scenarios?
- Scenario A:
- - WN only

- Scenario B:
- - WN + RN

- Scenario C:
- - WN + RN + GWB

-----------------------

Question 2: How do pulsar distance constraints improve as you add CW information?

How to do this?

Setup?
- Base this off the Q1 work, but output pulsar distance posteriors.

- Ring PTA again

Signals?
- N phase connected CWs
- Noise matching Scenarios A, B, C

- CW freq structure:
- - Distribution of freqs
- - Use the minimum frequency separation of the sources?

- Geometry:
- - Distribution of cosmu
- - Measure as mean lever arm over sources?

- SNRs:
- - Same as Q1 ~ 5, 10, 20

- Distance priors:
- - Same as Q1

Across these, measure:
- A gain factor of prior vs posterior fractional uncertainty for each pulsar.
- - The distribution of this is informative as a function of: S/N, N_CW, freq structure, geometry
- Bias for each pulsar in the ring array? Does that make sense?

- Joint PTA CW model vs factorised per-CW model
- - Does factorised approach give systematically larger gains in distance inference? I think yes.


-----------------------

Question 3: What changes when a phase-connected pulsar-term model cannot be used?

*** Okay do we actually want to ask this? ***

-----------------------

Question 4: How does precise distance knowledge for some pulsars help infer the distances of others?


-----------------------

Question 5: What is the best way to sample the CW nodes effectively?

- Contained in the CW_node_sampling directory, this has some good progress.


-----------------------
